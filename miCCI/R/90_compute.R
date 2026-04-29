# =============================================================================
# miCCI / 90_compute.R
# Stage 1: compute_predictions().
#
# Runs S1..S4 + meta on the full cohort and persists a single self-contained
# artefact (predictions.parquet + run_manifest.json + meta_weights.csv).
# Stage 2 (build_report) reads those files and never recomputes anything.
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
                        date_from = "2010-01-01",
                        date_to   = "2024-09-30") {
  if (!isTRUE(requireNamespace("arrow", quietly = TRUE)))
    stop("load_cohort() requires the 'arrow' package.")
  dt <- as.data.table(arrow::read_parquet(path))
  keep <- c("falnr", "age", "date_admission", "date_discharge",
            "stay_in_days", "diagnosen")
  miss <- setdiff(keep, names(dt))
  if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))
  dt <- dt[, ..keep]
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

  n_before <- nrow(dt)
  dt <- unique(dt, by = "falnr")
  n_after <- nrow(dt)
  if (n_before != n_after) {
    message(sprintf("load_cohort: %d duplicate falnr rows dropped (readmissions hidden by uniqueness constraint)",
                    n_before - n_after))
  }
  message(sprintf("Cohort: %d encounters (%d-%d)",
                  nrow(dt), min(dt$year), max(dt$year)))
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
.cohort_fingerprint <- function(dt) {
  s <- paste(nrow(dt),
             dt$falnr[1L], dt$falnr[nrow(dt)],
             paste(head(dt$diagnosen, 100L), collapse = "|"),
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
                                resume      = TRUE) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  message("=== miCCI Stage 1 - compute_predictions ===")
  t_global <- proc.time()

  message("[1/8] Loading cohort")
  dt <- load_cohort(data_path)
  if (!is.null(sample_size) && nrow(dt) > sample_size) {
    dt <- .with_local_seed(seed, dt[sample(.N, sample_size)])
    message(sprintf("  Subsampled to %d encounters", sample_size))
  }

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

  fp <- .cohort_fingerprint(dt)

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
    res_s1 <- .timed("S1", cci_interval_batch(dt, quan_map, cache))
    .cp_save(output_dir, "s1", fp, res_s1)
  }
  dt[, c("s1_min", "s1_max", "s1_mid", "s1_width") :=
       .(res_s1$cci_min, res_s1$cci_max, res_s1$cci_mid, res_s1$interval_width)]

  message("[5/8] S2 Probabilistic (with per-group probabilities for QA)")
  res_s2 <- if (resume) .cp_load(output_dir, "s2", fp) else NULL
  if (is.null(res_s2)) {
    res_s2 <- .timed("S2",
                     cci_probabilistic_batch(dt[, .(diagnosen)],
                                             quan_map, cache,
                                             return_group_prob = TRUE))
    .cp_save(output_dir, "s2", fp, res_s2)
  }
  dt[, s2_ecci := res_s2$e_cci]

  message("[6/8] S3 Multiple Imputation")
  res_s3 <- if (resume) .cp_load(output_dir, "s3", fp) else NULL
  if (is.null(res_s3)) {
    res_s3 <- .timed("S3",
                     cci_mi_batch(dt[, .(diagnosen)], quan_map, cache,
                                  m = mi_m, seed = seed,
                                  return_group_count = TRUE))
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
                                        return_group_count = TRUE))
    .cp_save(output_dir, "s4", fp, res_s4)
  }
  dt[, s4_bayes := res_s4$posterior_median]

  message("[8/8] Meta learner (cross-validated SuperLearner)")
  meta_res <- if (resume) .cp_load(output_dir, "meta", fp) else NULL
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

  out_cols <- c("falnr", "year", "age", "stay_in_days", "n_diagnoses",
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
  fwrite(res_s2$group_prob, file.path(output_dir, "s2_group_prob.csv"))
  fwrite(res_s3$group_prob, file.path(output_dir, "s3_group_prob.csv"))
  fwrite(res_s4$group_prob, file.path(output_dir, "s4_group_prob.csv"))
  message(sprintf("  predictions: %s  (n = %d, %d cols)",
                  pred_path, nrow(preds), ncol(preds)))

  manifest <- list(
    package         = "miCCI",
    version         = tryCatch(as.character(utils::packageVersion("miCCI")),
                               error = function(e) "1.0.0"),
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
      sample_size = if (is.null(sample_size)) "FULL" else sample_size
    ),
    meta_learner = list(
      method  = "SuperLearner (van der Laan, Polley, Hubbard 2007)",
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
  for (nm in c("s2_group_prob", "s3_group_prob", "s4_group_prob")) {
    f <- file.path(output_dir, paste0(nm, ".csv"))
    if (file.exists(f)) setattr(preds, nm, fread(f))
  }
  preds
}
