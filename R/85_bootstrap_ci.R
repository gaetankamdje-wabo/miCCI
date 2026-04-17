# =============================================================================
# miCCI v1.0.0 — 85_bootstrap_ci.R
# Non-parametric encounter-level bootstrap for MAE / RMSE / R^2 / Bias.
#
# Used both for the global per-strategy table and inside the stratified
# evaluation. Sequential by default; if the optional package future.apply
# is available the resampling loop runs in parallel.
# =============================================================================

#' Compute MAE / RMSE / R^2 / Bias on a single (predicted, gold) pair.
#' Returns numeric vector of length 4.
#' @keywords internal
.metrics4 <- function(pred, gold) {
  res <- pred - gold
  ss_tot <- sum((gold - mean(gold))^2)
  c(
    mae  = mean(abs(res)),
    rmse = sqrt(mean(res^2)),
    r2   = if (ss_tot > 0) 1 - sum(res^2) / ss_tot else NA_real_,
    bias = mean(res)
  )
}

#' Bootstrap 95% percentile CIs for MAE / RMSE / R^2 / Bias.
#'
#' @param pred numeric vector of predictions
#' @param gold numeric vector of gold values (same length)
#' @param B    number of bootstrap resamples (default 1000)
#' @param seed RNG seed
#' @param parallel If TRUE and future.apply is installed, resamples are
#'   parallelised across available cores.
#'
#' @return data.table with one row, columns:
#'   n, mae, mae_lo, mae_hi,
#'   rmse, rmse_lo, rmse_hi,
#'   r2, r2_lo, r2_hi,
#'   bias, bias_lo, bias_hi
#' @export
bootstrap_metrics <- function(pred, gold, B = 1000L, seed = 42L,
                              parallel = FALSE) {
  stopifnot(length(pred) == length(gold))
  ok   <- is.finite(pred) & is.finite(gold)
  pred <- as.numeric(pred[ok])
  gold <- as.numeric(gold[ok])
  n <- length(pred)
  if (n < 2L) {
    return(data.table(
      n = n,
      mae = NA_real_,  mae_lo  = NA_real_, mae_hi  = NA_real_,
      rmse = NA_real_, rmse_lo = NA_real_, rmse_hi = NA_real_,
      r2 = NA_real_,   r2_lo   = NA_real_, r2_hi   = NA_real_,
      bias = NA_real_, bias_lo = NA_real_, bias_hi = NA_real_
    ))
  }

  point <- .metrics4(pred, gold)

  one_replicate <- function(b) {
    idx <- sample.int(n, n, replace = TRUE)
    .metrics4(pred[idx], gold[idx])
  }

  set.seed(seed)
  if (parallel && requireNamespace("future.apply", quietly = TRUE)) {
    seeds <- sample.int(.Machine$integer.max, B)
    boot_mat <- future.apply::future_vapply(
      seq_len(B),
      function(b) { set.seed(seeds[b]); one_replicate(b) },
      FUN.VALUE = numeric(4),
      future.seed = TRUE
    )
  } else {
    boot_mat <- vapply(seq_len(B), one_replicate, FUN.VALUE = numeric(4))
  }
  # boot_mat is 4 x B
  q_lo <- apply(boot_mat, 1, function(x) quantile(x, 0.025, na.rm = TRUE))
  q_hi <- apply(boot_mat, 1, function(x) quantile(x, 0.975, na.rm = TRUE))

  data.table(
    n        = n,
    mae      = point["mae"],
    mae_lo   = q_lo["mae"],   mae_hi  = q_hi["mae"],
    rmse     = point["rmse"],
    rmse_lo  = q_lo["rmse"],  rmse_hi = q_hi["rmse"],
    r2       = point["r2"],
    r2_lo    = q_lo["r2"],    r2_hi   = q_hi["r2"],
    bias     = point["bias"],
    bias_lo  = q_lo["bias"],  bias_hi = q_hi["bias"]
  )
}

#' Apply bootstrap_metrics() to every strategy column of a predictions table.
#'
#' @param preds data.table with column cci_gold and the strategy columns
#' @param strategies named character vector — names = column names in preds,
#'                  values = pretty labels for the report
#' @export
bootstrap_strategies <- function(preds, strategies, B = 1000L, seed = 42L,
                                 parallel = FALSE) {
  rows <- list()
  for (i in seq_along(strategies)) {
    col <- names(strategies)[i]
    lab <- strategies[i]
    if (!col %in% names(preds)) next
    cat(sprintf("    bootstrap %-12s (B=%d)... ", lab, B))
    t0 <- proc.time()
    r <- bootstrap_metrics(preds[[col]], preds$cci_gold,
                           B = B, seed = seed, parallel = parallel)
    r[, strategy := col]
    r[, label    := lab]
    rows[[i]] <- r
    cat(sprintf("%.1fs\n", (proc.time() - t0)[3]))
  }
  out <- rbindlist(rows, use.names = TRUE)
  setcolorder(out, c("strategy", "label", "n",
                     "mae", "mae_lo", "mae_hi",
                     "rmse", "rmse_lo", "rmse_hi",
                     "r2", "r2_lo", "r2_hi",
                     "bias", "bias_lo", "bias_hi"))
  out[]
}
