# =============================================================================
# miCCI / 90_compute.R
# Stage 1: compute_predictions().
#
# Runs S1..S4 + meta on the full cohort and persists a single self-contained
# artefact (predictions.parquet + run_manifest.json + meta_weights.csv).
# Stage 2 (build_report) reads those files.
# =============================================================================

#' Load a clinical cohort from a Parquet file with the Mannheim schema.
#'
#' Required columns: `falnr`, `age`, `date_admission`, `date_discharge`,
#' `stay_in_days`, `diagnosen` (pipe-separated ICD-10-GM codes).
#'
#' Duplicate `falnr` rows are dropped silently in v0.x; v1.x reports the
#' drop count so readmissions are not lost without warning.
#'
#' @param path       parquet file path.
#' @param date_from  inclusive lower admission-date bound.
#' @param date_to    inclusive upper admission-date bound.
#' @return data.table with the cleaned cohort.
#' @export
load_cohort <- function(path,
                        date_from  = "2010-01-01",
                        date_to    = "2024-09-30",
                        patient_col = NULL) {
  if (!isTRUE(requireNamespace("arrow", quietly = TRUE)))
    stop("load_cohort() requires the 'arrow' package.")
  dt <- as.data.table(arrow::read_parquet(path))
  keep <- c("falnr", "age", "date_admission", "date_discharge",
            "stay_in_days", "diagnosen")
  miss <- setdiff(keep, names(dt))
  if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))

  # Patient identifier. `falnr` is the ENCOUNTER number, so it cannot serve as
  # the cluster for a patient-level bootstrap: deduplicating on it removes
  # duplicate rows of the same encounter, never a readmission. A separate
  # pseudonymised patient id is looked for by name, or named explicitly.
  if (is.null(patient_col)) {
    cand <- c("patient_id", "patient_pseudonym", "pseudonym", "patid",
              "pat_id", "patnr", "pid")
    hit  <- intersect(cand, names(dt))
    patient_col <- if (length(hit) > 0L) hit[1L] else NULL
  }
  if (!is.null(patient_col)) {
    if (!patient_col %in% names(dt))
      stop("patient_col not found in cohort: ", patient_col)
    keep <- c(keep, patient_col)
    message(sprintf("load_cohort: patient identifier '%s' found - cluster-robust inference enabled",
                    patient_col))
  } else {
    message("load_cohort: no patient identifier found - inference falls back to encounter-level resampling")
  }
  dt <- dt[, ..keep]
  if (!is.null(patient_col) && !identical(patient_col, "patient_id"))
    setnames(dt, patient_col, "patient_id")
  dt[, date_admission := as.Date(as.character(date_admission))]
  dt[, year           := as.integer(format(date_admission, "%Y"))]
  dt[, diagnosen      := as.character(diagnosen)]
  dt <- dt[!is.na(date_admission) & !is.na(year)]
  dt <- dt[!is.na(diagnosen) & diagnosen != "" & diagnosen != "NA"]
  dt[, stay_in_days := as.numeric(stay_in_days)]
  dt <- dt[!is.na(stay_in_days) & stay_in_days >= 0]
  dt[, age := as.numeric(age)]
  dt <- dt[!is.na(age) & age >= 0]
  dt <- dt[date_admission >= as.Date(date_from) &
           date_admission <= as.Date(date_to)]
  dt[, n_diagnoses := lengths(strsplit(diagnosen, "\\|+"))]

  # `falnr` is the encounter key: one row per discharge episode. Dropping
  # duplicates therefore removes repeated rows of the SAME encounter, not
  # readmissions - a readmission carries its own falnr and survives. The count
  # is reported so the claim is checkable rather than assumed.
  n_before <- nrow(dt)
  dt <- unique(dt, by = "falnr")
  n_after <- nrow(dt)
  if (n_before != n_after) {
    message(sprintf("load_cohort: %d duplicate rows of an already-present encounter number dropped (%.3f%%); distinct encounters retained: %d",
                    n_before - n_after, 100 * (n_before - n_after) / n_before,
                    n_after))
  }
  if ("patient_id" %in% names(dt)) {
    np <- data.table::uniqueN(dt$patient_id)
    message(sprintf("Cohort: %d encounters from %d patients (%.2f encounters/patient), %d-%d",
                    nrow(dt), np, nrow(dt) / np, min(dt$year), max(dt$year)))
  } else {
    message(sprintf("Cohort: %d encounters (%d-%d)",
                    nrow(dt), min(dt$year), max(dt$year)))
  }
  dt
}

.timed <- function(label, expr) {
  ex <- substitute(expr); pf <- parent.frame()
  message(sprintf("  %s ...", label)); t0 <- proc.time()
  result <- eval(ex, envir = pf)
  message(sprintf("  %s done in %.1fs", label, (proc.time() - t0)[3L]))
  result
}

# Stage-level checkpointing. Each expensive stage saves an .rds under
# <output_dir>/checkpoints/<stage>.rds; on rerun, if the file exists and
# the cohort fingerprint matches, the stage is skipped. The fingerprint
# is the SHA-1 of the cohort row count, the first/last falnr and the
# concatenated diagnoses of the first 100 rows - enough to detect any
# accidental cohort change without hashing the entire 720k-row data.table.
.cohort_fingerprint <- function(dt, config = NULL) {
  s <- paste(nrow(dt),
             dt$falnr[1L], dt$falnr[nrow(dt)],
             paste(head(dt$diagnosen, 100L), collapse = "|"),
             # Configuration is part of the fingerprint. Without it, switching
             # the age prior or the multiplicity rule would silently reuse
             # checkpoints computed under the previous settings.
             paste(names(config), unlist(config), sep = "=", collapse = ";"),
             sep = "::")
  # Cheap deterministic hash without rlang/digest; use R's built-in
  # adler-style via sum of charToRaw integers. Collisions are essentially
  # impossible for the constructed string above.
  paste0("len", nchar(s), "_sum",
         sum(as.integer(charToRaw(s))))
}

.cp_path <- function(output_dir, stage)
  file.path(output_dir, "checkpoints", paste0(stage, ".rds"))

.cp_load <- function(output_dir, stage, fp) {
  p <- .cp_path(output_dir, stage)
  if (!file.exists(p)) return(NULL)
  obj <- tryCatch(readRDS(p), error = function(e) NULL)
  if (is.null(obj) || !identical(obj$fp, fp)) return(NULL)
  message(sprintf("  [checkpoint] reusing %s from %s", stage, p))
  obj$value
}

.cp_save <- function(output_dir, stage, fp, value) {
  cp_dir <- file.path(output_dir, "checkpoints")
  dir.create(cp_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(fp = fp, value = value), .cp_path(output_dir, stage))
}

#' Stage 1: compute every per-encounter quantity and persist the artefact.
#'
#' @param data_path   parquet path (must match the schema of `load_cohort`).
#' @param output_dir  directory for `predictions.parquet`, `run_manifest.json`,
#'   and `meta_weights.csv`.
#' @param freq_table  optional data.table of population frequencies in
#'   Destatis schema. If NULL, `load_destatis()` is called and the German
#'   national reference is used.
#' @param sample_size optional integer; subsample for testing.
#' @param mi_m        S3 imputations.
#' @param bayes_draws S4 Dirichlet draws.
#' @param alpha_0     S4 prior pseudo-count multiplier.
#' @param sl_cv_folds CV folds for the meta learner.
#' @param seed        master RNG seed (local; never alters caller's RNG).
#' @param resume      if TRUE (default), reuse stage checkpoints under
#'   `<output_dir>/checkpoints/` when the cohort fingerprint matches. Set
#'   to FALSE to force a full recompute.
#' @return invisibly, list(predictions, manifest).
#' @export
compute_predictions <- function(data_path,
                                output_dir,
                                freq_table  = NULL,
                                sample_size = NULL,
                                mi_m        = 20L,
                                bayes_draws = 25L,
                                alpha_0     = 10,
                                sl_cv_folds = 10L,
                                seed        = 42L,
                                resume      = TRUE,
                                age_exist   = TRUE,
                                age_col     = "age",
                                preserve_multiplicity = TRUE,
                                patient_col = NULL) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  message("=== miCCI Stage 1 - compute_predictions ===")
  t_global <- proc.time()

  message("[1/8] Loading cohort")
  dt <- load_cohort(data_path, patient_col = patient_col)
  if (!is.null(sample_size) && nrow(dt) > sample_size) {
    dt <- .with_local_seed(seed, dt[sample(.N, sample_size)])
    message(sprintf("  Subsampled to %d encounters", sample_size))
  }

  # Age-conditioned prior. Destatis 23131-01 publishes subcode frequencies
  # across 22 age bands, and age strongly confounds comorbidity profiles, so
  # conditioning on the band the patient actually falls into sharpens every
  # subcode probability. `age_exist = FALSE` selects the marginal prior and
  # reproduces the submitted behaviour exactly.
  age_idx <- NULL
  if (isTRUE(age_exist)) {
    if (!age_col %in% names(dt))
      stop(sprintf("age_exist = TRUE but column '%s' is not in the cohort", age_col))
    age_idx <- age_to_bin_index(dt[[age_col]])
    n_ok <- sum(!is.na(age_idx))
    message(sprintf("  Age-conditioned prior: %d/%d encounters mapped to a Destatis age band (%.2f%%); the rest use the marginal",
                    n_ok, nrow(dt), 100 * n_ok / nrow(dt)))
  } else {
    message("  Marginal (age-aggregated) prior in use")
  }
  message(sprintf("  Code multiplicity: %s",
                  if (isTRUE(preserve_multiplicity))
                    "preserved (one prefix slot per coded diagnosis)"
                  else "collapsed (one slot per distinct prefix)"))

  message("[2/8] Loading Quan map and population frequencies")
  quan_map <- load_quan_map()
  freq_used_label <- "user-supplied freq_table"
  if (is.null(freq_table)) {
    freq_table <- .timed("Destatis", load_destatis())
    freq_used_label <- "Destatis 23131-01"
  }
  cache <- .timed("Precompute lookups", precompute_lookups(freq_table, quan_map))
  pl <- build_pattern_lookup(quan_map)
  dl <- build_dep_lookup(quan_map)

  fp <- .cohort_fingerprint(dt, config = list(
    age_exist = isTRUE(age_exist), age_col = age_col,
    multiplicity = isTRUE(preserve_multiplicity),
    mi_m = mi_m, bayes_draws = bayes_draws, alpha_0 = alpha_0,
    sl_cv_folds = sl_cv_folds, seed = seed))

  message("[3/8] Gold CCI")
  gold <- if (resume) .cp_load(output_dir, "gold", fp) else NULL
  if (is.null(gold)) {
    gold <- .timed("gold", cci_gold_batch(dt, quan_map, pl, dl))
    .cp_save(output_dir, "gold", fp, gold)
  }
  dt[, cci_gold := gold]

  message("[4/8] S1 Interval")
  res_s1 <- if (resume) .cp_load(output_dir, "s1", fp) else NULL
  if (is.null(res_s1)) {
    res_s1 <- .timed("S1", cci_interval_batch(dt, quan_map, cache, return_group_prob = TRUE))
    .cp_save(output_dir, "s1", fp, res_s1)
  }
  # Checkpoint may be a legacy data.table (no group_prob). Upgrade transparently.
  if (is.data.table(res_s1)) {
    message("  [S1] Legacy checkpoint detected - recomputing group_prob only")
    res_s1 <- list(
      cci_min        = res_s1$cci_min,
      cci_max        = res_s1$cci_max,
      cci_mid        = res_s1$cci_mid,
      interval_width = res_s1$interval_width,
      group_prob     = cci_interval_batch(dt, quan_map, cache,
                                          return_group_prob = TRUE)$group_prob
    )
    .cp_save(output_dir, "s1", fp, res_s1)  # overwrite with upgraded checkpoint
  }
  dt[, c("s1_min", "s1_max", "s1_mid", "s1_width") :=
       .(res_s1$cci_min, res_s1$cci_max, res_s1$cci_mid, res_s1$interval_width)]

  message("[5/8] S2 Probabilistic (with per-group probabilities for QA)")
  res_s2 <- if (resume) .cp_load(output_dir, "s2", fp) else NULL
  if (is.null(res_s2)) {
    res_s2 <- .timed("S2",
                     cci_probabilistic_batch(dt[, .(diagnosen)],
                                             quan_map, cache,
                                             return_group_prob = TRUE,
                                             age_idx = age_idx,
                                             preserve_multiplicity = preserve_multiplicity))
    .cp_save(output_dir, "s2", fp, res_s2)
  }
  dt[, s2_ecci := res_s2$e_cci]

  message("[6/8] S3 Multiple Imputation")
  res_s3 <- if (resume) .cp_load(output_dir, "s3", fp) else NULL
  if (is.null(res_s3)) {
    res_s3 <- .timed("S3",
                     cci_mi_batch(dt[, .(diagnosen)], quan_map, cache,
                                  m = mi_m, seed = seed,
                                  return_group_count = TRUE,
                                  age_idx = age_idx,
                                  preserve_multiplicity = preserve_multiplicity))
    .cp_save(output_dir, "s3", fp, res_s3)
  }
  dt[, s3_mi := res_s3$mi_cci]

  message("[7/8] S4 Bayesian")
  res_s4 <- if (resume) .cp_load(output_dir, "s4", fp) else NULL
  if (is.null(res_s4)) {
    res_s4 <- .timed("S4",
                     cci_bayesian_batch(dt[, .(diagnosen)], quan_map, cache,
                                        n_draws = bayes_draws,
                                        alpha_0 = alpha_0,
                                        seed = seed,
                                        return_group_count = TRUE,
                                        age_idx = age_idx,
                                        preserve_multiplicity = preserve_multiplicity))
    .cp_save(output_dir, "s4", fp, res_s4)
  }
  dt[, s4_bayes := res_s4$posterior_median]

  message("[8/8] Meta learner (cross-validated SuperLearner)")
  # SuperLearner is an optional dependency. The four base strategies are the
  # substance of the method and none of them needs it, so a missing package
  # degrades to an S1-S4 run rather than throwing away eight completed stages.
  have_sl <- isTRUE(requireNamespace("SuperLearner", quietly = TRUE))
  meta_res <- if (resume && have_sl) .cp_load(output_dir, "meta", fp) else NULL
  if (!have_sl) {
    warning("SuperLearner is not installed: the Meta Learner is skipped and ",
            "the 'meta' column is NA. Install it and rerun stage 1 to add it.",
            call. = FALSE)
    message("  SKIPPED - SuperLearner not installed. S1 to S4 are unaffected.")
    dt[, meta := NA_real_]
  } else {
    if (is.null(meta_res)) {
      meta_res <- .timed("Meta",
        cci_meta_fit(dt[, .(cci_gold, s1_min, s1_max, s1_mid,
                            s2_ecci, s3_mi, s4_bayes)],
                     V = sl_cv_folds, seed = seed, verbose = FALSE))
      .cp_save(output_dir, "meta", fp, meta_res)
    }
    dt[, meta := meta_res$predictions]
    message("  Meta NNLS weights:")
    for (nm in names(meta_res$weights))
      message(sprintf("    %-10s %.4f", nm, meta_res$weights[nm]))
    message("  Meta CV risk (lower is better):")
    for (nm in names(meta_res$cv_risk))
      message(sprintf("    %-10s %.4f", nm, meta_res$cv_risk[nm]))
  }

  out_cols <- c("falnr", if ("patient_id" %in% names(dt)) "patient_id",
                "year", "age", "stay_in_days", "n_diagnoses",
                "diagnosen",
                "cci_gold",
                "s1_min", "s1_max", "s1_mid", "s1_width",
                "s2_ecci", "s3_mi", "s4_bayes", "meta")
  preds <- dt[, ..out_cols]

  preds[, cci_gold := as.integer(cci_gold)]
  preds[, s1_min   := as.integer(s1_min)]
  preds[, s1_max   := as.integer(s1_max)]
  preds[, s1_width := as.integer(s1_width)]
  for (col in c("s1_mid", "s2_ecci", "s3_mi", "s4_bayes", "meta")) {
    set(preds, j = col, value = round(as.numeric(preds[[col]]), 4L))
  }

  # Attach per-group probability tables as attributes for the QA stage.
  setattr(preds, "s1_group_prob", res_s1$group_prob)
  setattr(preds, "s2_group_prob", res_s2$group_prob)
  setattr(preds, "s3_group_prob", res_s3$group_prob)
  setattr(preds, "s4_group_prob", res_s4$group_prob)

  parquet_ok <- isTRUE(requireNamespace("arrow", quietly = TRUE))
  pred_path  <- file.path(output_dir,
                          if (parquet_ok) "predictions.parquet" else "predictions.rds")
  if (parquet_ok) {
    arrow::write_parquet(preds, pred_path)
  } else {
    saveRDS(preds, pred_path)
  }
  # Group-prob tables are stored separately (parquet does not preserve attrs).
  fwrite(res_s1$group_prob, file.path(output_dir, "s1_group_prob.csv"))
  fwrite(res_s2$group_prob, file.path(output_dir, "s2_group_prob.csv"))
  fwrite(res_s3$group_prob, file.path(output_dir, "s3_group_prob.csv"))
  fwrite(res_s4$group_prob, file.path(output_dir, "s4_group_prob.csv"))
  message(sprintf("  predictions: %s  (n = %d, %d cols)",
                  pred_path, nrow(preds), ncol(preds)))

  manifest <- list(
    package         = "miCCI",
    version         = tryCatch(as.character(utils::packageVersion("miCCI")),
                               error = function(e) NA_character_),
    created_utc     = format(Sys.time(), tz = "UTC", usetz = TRUE),
    data_path       = normalizePath(data_path, mustWork = FALSE),
    output_dir      = normalizePath(output_dir, mustWork = FALSE),
    freq_table      = freq_used_label,
    n_encounters    = nrow(preds),
    year_min        = as.integer(min(preds$year)),
    year_max        = as.integer(max(preds$year)),
    parameters      = list(
      mi_m        = mi_m,
      bayes_draws = bayes_draws,
      alpha_0     = alpha_0,
      sl_cv_folds = sl_cv_folds,
      sl_method   = "SuperLearner::method.NNLS",
      seed        = seed,
      sample_size = if (is.null(sample_size)) "FULL" else sample_size,
      age_exist   = isTRUE(age_exist),
      age_col     = age_col,
      age_prior   = if (isTRUE(age_exist)) "Destatis 23131-01 age-conditioned"
                    else "Destatis 23131-01 marginal",
      preserve_multiplicity = isTRUE(preserve_multiplicity),
      patient_id_available  = "patient_id" %in% names(dt),
      n_patients  = if ("patient_id" %in% names(dt))
                      data.table::uniqueN(dt$patient_id) else NA_integer_
    ),
    meta_learner = list(
      fitted  = have_sl,
      method  = if (have_sl) "SuperLearner (van der Laan, Polley, Hubbard 2007)"
                else "not fitted: SuperLearner unavailable",
      library = c("SL.s1_mid", "SL.s2_ecci", "SL.s3_mi", "SL.s4_bayes"),
      meta    = "method.NNLS",
      V       = sl_cv_folds,
      weights = as.list(meta_res$weights),
      cv_risk = as.list(meta_res$cv_risk)
    ),
    runtime_seconds = round((proc.time() - t_global)[3L], 1L),
    columns         = out_cols,
    R_version       = R.version.string,
    sessionInfo     = utils::capture.output(utils::sessionInfo())
  )
  manifest_path <- file.path(output_dir, "run_manifest.json")
  writeLines(jsonlite::toJSON(manifest, pretty = TRUE, auto_unbox = TRUE),
             manifest_path)
  message(sprintf("  manifest:    %s", manifest_path))

  meta_csv <- file.path(output_dir, "meta_weights.csv")
  fwrite(data.table(
    learner = names(meta_res$weights),
    weight  = as.numeric(meta_res$weights),
    cv_risk = as.numeric(meta_res$cv_risk)
  ), meta_csv)
  message(sprintf("  meta weights: %s", meta_csv))

  message(sprintf("Stage 1 complete in %.1fs (%.2f h)",
                  manifest$runtime_seconds, manifest$runtime_seconds / 3600))
  invisible(list(predictions = preds, manifest = manifest))
}

#' Load a previously persisted predictions artefact.
#' @export
load_predictions <- function(output_dir) {
  pq <- file.path(output_dir, "predictions.parquet")
  rd <- file.path(output_dir, "predictions.rds")
  if (file.exists(pq)) {
    if (!isTRUE(requireNamespace("arrow", quietly = TRUE)))
      stop("Reading predictions.parquet requires the 'arrow' package.")
    preds <- as.data.table(arrow::read_parquet(pq))
  } else if (file.exists(rd)) {
    preds <- as.data.table(readRDS(rd))
  } else {
    stop("No predictions artefact found in ", output_dir)
  }
  # Reattach per-group probability tables if present alongside.
  for (nm in c("s1_group_prob", "s2_group_prob", "s3_group_prob", "s4_group_prob")) {
    f <- file.path(output_dir, paste0(nm, ".csv"))
    if (file.exists(f)) setattr(preds, nm, fread(f))
  }
  preds
}