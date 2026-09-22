# =============================================================================
# miCCI / 85_bootstrap_ci.R
# Bootstrap inference for MAE / RMSE / R^2 / Bias.
#
# CLUSTER-ROBUST RESAMPLING (added in this revision)
# -------------------------------------------
# Previously the bootstrap drew encounters independently. A patient
# readmitted several times contributes several encounters whose reconstruction
# errors are correlated, because they share the same coding habits, the same
# chronic diagnoses and often the same prefixes. Treating them as independent
# units understates the variance of every metric and yields confidence
# intervals that are too narrow.
#
# When a cluster identifier is supplied, resampling happens at the patient
# level: patients are drawn with replacement and every encounter of a drawn
# patient enters the replicate together. That is the standard cluster (block)
# bootstrap and it is valid under arbitrary within-patient dependence. Cluster
# sizes are preserved, so a replicate has a random number of encounters.
#
# PAIRED DIFFERENCES
# ------------------
# Comparing two strategies through their marginal intervals ignores that both
# are evaluated on the identical encounters, which induces a strong positive
# correlation between their errors. `paired_bootstrap_diff()` resamples once
# per replicate and recomputes BOTH strategies on that replicate, so the
# interval is for the difference itself. It honours the same cluster argument.
# =============================================================================

#' Compute MAE / RMSE / R^2 / Bias for one (predicted, gold) pair.
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

#' Build a reusable resampler for one cohort.
#'
#' Returns a function of no arguments that yields a vector of row indices for
#' one bootstrap replicate. With `cluster = NULL` this is the ordinary
#' encounter-level bootstrap; otherwise it is the cluster bootstrap over the
#' levels of `cluster`. The index splitting happens once, not once per
#' replicate, which is what keeps B = 1000 cluster replicates affordable on a
#' half-million-row cohort.
#' @keywords internal
.make_resampler <- function(n, cluster = NULL) {
  if (is.null(cluster)) {
    return(function() sample.int(n, n, replace = TRUE))
  }
  groups <- split(seq_len(n), cluster)
  G <- length(groups)
  function() unlist(groups[sample.int(G, G, replace = TRUE)],
                    use.names = FALSE)
}

#' Percentile bootstrap CIs for MAE / RMSE / R^2 / Bias.
#'
#' @param pred numeric predictions.
#' @param gold numeric reference values.
#' @param B number of bootstrap replicates.
#' @param seed RNG seed (local; never alters the caller's RNG state).
#' @param cluster optional vector of cluster identifiers, one per element of
#'   `pred`, typically the pseudonymised patient id. When supplied, replicates
#'   are drawn at the cluster level so within-patient correlation from
#'   readmissions is carried into the interval.
#' @param conf confidence level (default 0.95). Values below 0.95 widen the
#'   interval and are how a multiplicity correction is applied: pass
#'   `1 - alpha/K` for a Bonferroni-adjusted family of K comparisons.
#' @param parallel parallelise replicates when `future.apply` is available.
#' @export
bootstrap_metrics <- function(pred, gold, B = 1000L, seed = 42L,
                              cluster = NULL, conf = 0.95,
                              parallel = FALSE) {
  stopifnot(length(pred) == length(gold))
  ok   <- is.finite(pred) & is.finite(gold)
  if (!is.null(cluster)) {
    stopifnot(length(cluster) == length(pred))
    cluster <- cluster[ok]
  }
  pred <- as.numeric(pred[ok])
  gold <- as.numeric(gold[ok])
  n <- length(pred)
  if (n < 2L) {
    return(data.table(
      n = n, n_clusters = NA_integer_,
      mae = NA_real_,  mae_lo  = NA_real_, mae_hi  = NA_real_,
      rmse = NA_real_, rmse_lo = NA_real_, rmse_hi = NA_real_,
      r2 = NA_real_,   r2_lo   = NA_real_, r2_hi   = NA_real_,
      bias = NA_real_, bias_lo = NA_real_, bias_hi = NA_real_
    ))
  }

  point <- .metrics4(pred, gold)
  draw  <- .make_resampler(n, cluster)
  one_replicate <- function(b) {
    idx <- draw()
    .metrics4(pred[idx], gold[idx])
  }

  boot_mat <- .with_local_seed(seed, {
    if (parallel && requireNamespace("future.apply", quietly = TRUE)) {
      seeds <- sample.int(.Machine$integer.max, B)
      future.apply::future_vapply(
        seq_len(B),
        function(b) { set.seed(seeds[b]); one_replicate(b) },
        FUN.VALUE = numeric(4L),
        future.seed = TRUE
      )
    } else {
      vapply(seq_len(B), one_replicate, FUN.VALUE = numeric(4L))
    }
  })

  a <- (1 - conf) / 2
  q_lo <- apply(boot_mat, 1L, function(x) stats::quantile(x, a,     na.rm = TRUE))
  q_hi <- apply(boot_mat, 1L, function(x) stats::quantile(x, 1 - a, na.rm = TRUE))

  data.table(
    n          = n,
    n_clusters = if (is.null(cluster)) NA_integer_ else data.table::uniqueN(cluster),
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

#' Apply `bootstrap_metrics()` to every strategy column.
#'
#' @param preds data.table of predictions plus `cci_gold`.
#' @param strategies named character vector (column name = pretty label).
#' @param cluster_col optional column name in `preds` holding the cluster
#'   (patient) identifier. When present, all intervals are cluster-robust.
#' @inheritParams bootstrap_metrics
#' @export
bootstrap_strategies <- function(preds, strategies, B = 1000L, seed = 42L,
                                 cluster_col = NULL, conf = 0.95,
                                 parallel = FALSE) {
  cl <- if (!is.null(cluster_col) && cluster_col %in% names(preds))
    preds[[cluster_col]] else NULL
  rows <- list()
  for (i in seq_along(strategies)) {
    col <- names(strategies)[i]
    lab <- strategies[i]
    if (!col %in% names(preds)) next
    r <- bootstrap_metrics(preds[[col]], preds$cci_gold,
                           B = B, seed = seed, cluster = cl, conf = conf,
                           parallel = parallel)
    r[, strategy := col]
    r[, label    := lab]
    rows[[i]] <- r
  }
  out <- rbindlist(rows, use.names = TRUE)
  setcolorder(out, c("strategy", "label", "n", "n_clusters",
                     "mae", "mae_lo", "mae_hi",
                     "rmse", "rmse_lo", "rmse_hi",
                     "r2", "r2_lo", "r2_hi",
                     "bias", "bias_lo", "bias_hi"))
  out[]
}

#' Paired bootstrap of the MAE difference between two strategies.
#'
#' Both strategies are recomputed on the same replicate, so the interval is for
#' the difference and absorbs the cross-strategy correlation that evaluation on
#' identical encounters induces. Comparing two marginal intervals instead is
#' conservative in the wrong direction: it can call a real difference a tie.
#'
#' With a cluster, each replicate draws patients with replacement. The MAE
#' difference is a ratio of sums, so a replicate needs only the per-patient
#' totals of the absolute-error difference and the per-patient encounter
#' counts: `sum(D[j]) / sum(N[j])` over the drawn patients `j`. The patient
#' indices come from the same `sample.int(G, G, replace = TRUE)` call as before,
#' so the replicates are unchanged, only no longer expanded back into
#' encounters on every draw. On the full cohort that is the difference between
#' minutes and hours.
#'
#' @param pred_a,pred_b numeric predictions of the two strategies.
#' @param gold numeric reference values.
#' @param B replicates.
#' @param seed RNG seed (local).
#' @param cluster optional cluster identifiers (see `bootstrap_metrics()`).
#' @param conf confidence level of the primary interval.
#' @param conf_adj optional second, adjusted level (e.g. Bonferroni
#'   `1 - alpha / K`), read from the SAME replicates. Percentile intervals from
#'   one replicate vector are nested, so the adjusted interval always contains
#'   the primary one and an adjusted "decisive" implies a primary "decisive".
#' @return one-row data.table: observed difference (a minus b), its interval,
#'   `decisive` when the interval excludes zero, and the adjusted counterparts
#'   when `conf_adj` is given.
#' @export
paired_bootstrap_diff <- function(pred_a, pred_b, gold, B = 1000L, seed = 42L,
                                  cluster = NULL, conf = 0.95, conf_adj = NULL) {
  ok <- is.finite(pred_a) & is.finite(pred_b) & is.finite(gold)
  if (!is.null(cluster)) cluster <- cluster[ok]
  a <- as.numeric(pred_a[ok]); b <- as.numeric(pred_b[ok])
  g <- as.numeric(gold[ok])
  n <- length(g)
  na_out <- data.table(diff = NA_real_, diff_lo = NA_real_, diff_hi = NA_real_,
                       conf = conf, B = B, n = n, decisive = NA)
  if (!is.null(conf_adj))
    na_out[, `:=`(diff_lo_adj = NA_real_, diff_hi_adj = NA_real_,
                  conf_adj = conf_adj, decisive_adj = NA)]
  if (n < 2L) return(na_out)

  d   <- abs(a - g) - abs(b - g)
  obs <- mean(d)

  reps <- if (is.null(cluster)) {
    .with_local_seed(seed, vapply(seq_len(B), function(i) {
      idx <- sample.int(n, n, replace = TRUE)
      mean(d[idx])
    }, numeric(1L)))
  } else {
    # split() orders groups exactly as the encounter-expanding resampler did,
    # so index j refers to the same patient in both formulations.
    D  <- vapply(split(d, cluster), sum, numeric(1L))
    Nc <- vapply(split(d, cluster), length, integer(1L))
    G  <- length(D)
    .with_local_seed(seed, vapply(seq_len(B), function(i) {
      j <- sample.int(G, G, replace = TRUE)
      sum(D[j]) / sum(Nc[j])
    }, numeric(1L)))
  }

  ci <- function(level) {
    a2 <- (1 - level) / 2
    unname(stats::quantile(reps, c(a2, 1 - a2), na.rm = TRUE))
  }
  q <- ci(conf)
  out <- data.table(diff = obs, diff_lo = q[1L], diff_hi = q[2L],
                    conf = conf, B = B, n = n,
                    decisive = (q[1L] > 0) | (q[2L] < 0))
  if (!is.null(conf_adj)) {
    qa <- ci(conf_adj)
    out[, `:=`(diff_lo_adj = qa[1L], diff_hi_adj = qa[2L],
               conf_adj = conf_adj,
               decisive_adj = (qa[1L] > 0) | (qa[2L] < 0))]
  }
  out
}
