# =============================================================================
# miCCI v1.0.0 — 90_compute.R
# Stage 1 of the two-stage pipeline.
#
# compute_predictions() runs the full long-running computation (S1-S4 + S6,
# vectorised) on the entire consistent 2010-2024 cohort and writes a single
# self-contained artefact (predictions.parquet + run_manifest.json) that can
# later be consumed by build_report() (stage 2) without ever recomputing.
#
# No train/calibrate/validate split. No ML. No NNLS training. No rounding:
# all probabilistic outputs are floats.
# =============================================================================

#' Load the consistent Mannheim cohort 2010-2024.
#'
#' Identical cleaning to v0.4.3 load_mannheim_data(), but the entire interval
#' is treated as one cohort — no temporal split.
#'
#' @export
load_cohort <- function(path,
                        date_from = "2010-01-01",
                        date_to   = "2024-09-30") {
  requireNamespace("arrow", quietly = TRUE)
  dt <- as.data.table(arrow::read_parquet(path))
  keep <- c("falnr", "age", "date_admission", "date_discharge",
            "stay_in_days", "diagnosen")
  missing <- setdiff(keep, names(dt))
  if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))
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
  los_r <- frank(dt$stay_in_days, ties.method = "dense")
  dt[, los_decile := as.integer(ceiling(10 * los_r / max(los_r)))]
  dt[los_decile > 10L, los_decile := 10L]
  dt <- unique(dt, by = "falnr")
  message(sprintf("\u2705 cohort: %d encounters (%d-%d)",
                  nrow(dt), min(dt$year), max(dt$year)))
  dt
}

.timed <- function(label, expr) {
  cat(sprintf("  %s... ", label)); t0 <- proc.time()
  result <- eval(expr); cat(sprintf("%.1fs\n", (proc.time() - t0)[3]))
  result
}

#' Stage 1 — compute every per-encounter quantity and persist.
#'
#' This is the >16h job. Run once.
#'
#' @param data_path     Parquet path (Mannheim or compatible schema)
#' @param output_dir    Directory for predictions.parquet + run_manifest.json
#' @param sample_size   Optional integer; subsample for testing only (NULL=all)
#' @param mi_m          Number of MI imputations (S3); default 20
#' @param bayes_draws   Number of Bayesian Dirichlet draws (S4); default 25
#' @param sl_cv_folds   Number of CV folds for the S6 Super Learner; default 10
#' @param seed          RNG seed for S3/S4/S6 reproducibility; default 42
#'
#' @return Invisibly: list(predictions = <data.table>, manifest = <list>).
#'         Side effect: writes predictions.parquet (or .rds fallback) +
#'         run_manifest.json + a console summary.
#'
#' @export
compute_predictions <- function(data_path,
                                output_dir,
                                sample_size = NULL,
                                mi_m        = 20L,
                                bayes_draws = 25L,
                                sl_cv_folds = 10L,
                                seed        = 42L) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  cat("=== miCCI v1.0.0 — Stage 1: compute_predictions ===\n\n")
  t_global <- proc.time()

  # 1. Load
  cat("[1/7] Loading cohort\n")
  dt <- load_cohort(data_path)
  if (!is.null(sample_size) && nrow(dt) > sample_size) {
    set.seed(seed)
    dt <- dt[sample(.N, sample_size)]
    cat(sprintf("  Subsampled to %d encounters\n", sample_size))
  }
  n <- nrow(dt)

  # 2. Precompute Destatis lookups
  cat("[2/7] Precomputing Destatis + Quan lookups\n")
  quan_map    <- load_quan_map()
  dt_destatis <- .timed("Destatis", load_destatis())
  cache       <- .timed("Cache",    precompute_lookups(dt_destatis, quan_map))
  pl          <- build_pattern_lookup(quan_map)
  dl          <- build_dep_lookup(quan_map)

  # 3. Gold CCI (vectorised)
  cat("[3/7] Gold CCI (reference)\n")
  dt[, cci_gold := .timed("gold", cci_gold_batch(.SD, quan_map, pl, dl))]

  # 4. S1 — Interval
  cat("[4/7] S1 Interval\n")
  res <- .timed("S1", cci_interval_batch(dt, quan_map, cache))
  dt[, c("s1_min", "s1_max", "s1_mid", "s1_width") :=
       .(res$cci_min, res$cci_max, res$cci_mid, res$interval_width)]

  # 5. S2 — Probabilistic
  cat("[5/7] S2 Probabilistic\n")
  dt[, s2_ecci := .timed("S2", cci_probabilistic_batch(.SD, quan_map, cache))]

  # 6. S3 — Multiple Imputation
  cat("[6/7] S3 Multiple Imputation\n")
  dt[, s3_mi := .timed("S3",
                       cci_mi_batch(.SD, quan_map, cache,
                                    m = mi_m, seed = seed))]

  # 7. S4 — Bayesian
  cat("[7/7] S4 Bayesian\n")
  dt[, s4_bayes := .timed("S4",
                          cci_bayesian_batch(.SD, quan_map, cache,
                                             n_draws = bayes_draws,
                                             seed = seed))]

  # 8. S6 — Cross-validated Super Learner meta-estimator
  cat("[+]   S6 Super Learner meta-estimator (10-fold CV)\n")
  sl_res <- .timed("S6",
    cci_meta_fit(dt[, .(cci_gold, s1_min, s1_max, s1_mid,
                        s2_ecci, s3_mi, s4_bayes)],
                 V = sl_cv_folds, seed = seed, verbose = FALSE))
  dt[, s6_meta := sl_res$predictions]
  cat("  SL weights (NNLS, sum-to-1):\n")
  for (nm in names(sl_res$weights))
    cat(sprintf("    %-10s %.4f\n", nm, sl_res$weights[nm]))
  cat("  SL CV risk per learner:\n")
  for (nm in names(sl_res$cv_risk))
    cat(sprintf("    %-10s %.4f\n", nm, sl_res$cv_risk[nm]))

  # 9. Persist artefact
  cat("\n=== SAVING ARTEFACT ===\n")
  out_cols <- c("falnr", "year", "age", "stay_in_days", "n_diagnoses",
                "los_decile", "diagnosen",
                "cci_gold",
                "s1_min", "s1_max", "s1_mid", "s1_width",
                "s2_ecci", "s3_mi", "s4_bayes", "s6_meta")
  preds <- dt[, ..out_cols]

  # Output rounding: 2 decimal places for every *predicted* CCI column.
  # The gold standard stays integer by definition; S1 min/max are integer
  # sums of integer Charlson weights; everything else is a float that we
  # round only for storage so the artefact is human-readable and compact.
  # Rounding happens ONCE here, on the final output, so the internal
  # computation (S2 expectation, S3 imputation mean, S4 posterior median,
  # S6 meta) was carried out in full double precision.
  preds[, cci_gold := as.integer(cci_gold)]
  preds[, s1_min   := as.integer(s1_min)]
  preds[, s1_max   := as.integer(s1_max)]
  preds[, s1_width := as.integer(s1_width)]
  for (col in c("s1_mid", "s2_ecci", "s3_mi", "s4_bayes", "s6_meta")) {
    set(preds, j = col, value = round(as.numeric(preds[[col]]), 2))
  }

  parquet_ok <- requireNamespace("arrow", quietly = TRUE)
  pred_path  <- file.path(output_dir,
                          if (parquet_ok) "predictions.parquet" else "predictions.rds")
  if (parquet_ok) {
    arrow::write_parquet(preds, pred_path)
  } else {
    saveRDS(preds, pred_path)
  }
  cat(sprintf("  predictions: %s  (n = %d, %d cols)\n",
              pred_path, nrow(preds), ncol(preds)))

  # 10. Run manifest (machine + human readable)
  manifest <- list(
    package         = "miCCI",
    version         = "1.0.0",
    created_utc     = format(Sys.time(), tz = "UTC", usetz = TRUE),
    data_path       = normalizePath(data_path, mustWork = FALSE),
    output_dir      = normalizePath(output_dir, mustWork = FALSE),
    n_encounters    = nrow(preds),
    year_min        = as.integer(min(preds$year)),
    year_max        = as.integer(max(preds$year)),
    parameters      = list(
      mi_m        = mi_m,
      bayes_draws = bayes_draws,
      sl_cv_folds = sl_cv_folds,
      sl_method   = "SuperLearner::method.NNLS",
      seed        = seed,
      sample_size = if (is.null(sample_size)) "FULL" else sample_size
    ),
    super_learner   = list(
      method  = "SuperLearner (van der Laan, Polley, Hubbard 2007)",
      library = c("SL.s1_mid", "SL.s2_ecci", "SL.s3_mi", "SL.s4_bayes"),
      meta    = "method.NNLS",
      V       = sl_cv_folds,
      weights = as.list(sl_res$weights),
      cv_risk = as.list(sl_res$cv_risk)
    ),
    runtime_seconds = round((proc.time() - t_global)[3], 1),
    columns         = out_cols,
    R_version       = R.version.string,
    sessionInfo     = utils::capture.output(utils::sessionInfo())
  )
  manifest_path <- file.path(output_dir, "run_manifest.json")
  writeLines(jsonlite::toJSON(manifest, pretty = TRUE, auto_unbox = TRUE),
             manifest_path)
  cat(sprintf("  manifest:    %s\n", manifest_path))

  # Also drop a dedicated CSV with the SL weights — handy for the paper
  sl_csv <- file.path(output_dir, "superlearner_weights.csv")
  fwrite(data.table(
    learner = names(sl_res$weights),
    weight  = as.numeric(sl_res$weights),
    cv_risk = as.numeric(sl_res$cv_risk)
  ), sl_csv)
  cat(sprintf("  sl weights:  %s\n", sl_csv))

  cat(sprintf("\nTotal stage-1 runtime: %.1f s (%.2f h)\n",
              manifest$runtime_seconds, manifest$runtime_seconds / 3600))
  cat("=== STAGE 1 DONE ===\n")
  cat("Next step: Rscript scripts/02_run_report.R\n\n")

  invisible(list(predictions = preds, manifest = manifest))
}

#' Helper: load a previously written predictions artefact.
#' @export
load_predictions <- function(output_dir) {
  pq <- file.path(output_dir, "predictions.parquet")
  rd <- file.path(output_dir, "predictions.rds")
  if (file.exists(pq)) {
    requireNamespace("arrow", quietly = TRUE)
    return(as.data.table(arrow::read_parquet(pq)))
  }
  if (file.exists(rd)) return(as.data.table(readRDS(rd)))
  stop("No predictions artefact found in ", output_dir,
       " (looked for predictions.parquet and predictions.rds).")
}
