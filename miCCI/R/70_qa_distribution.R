# =============================================================================
# miCCI / 70_qa_distribution.R
# Quality-assurance / distributional plausibility analysis.
#
# Three deliverables:
#   (a) JSON coverage   - per Charlson group, the fraction of the cohort
#                          for which the gold algorithm flags it active.
#   (b) Mass conservation - per Charlson group, the mean expected
#                          contribution to the CCI under each strategy
#                          divided by the gold contribution. A faithful
#                          strategy gives ratios near 1 across groups.
#                          IMPORTANT: this requires per-(idx, gk)
#                          probabilities for the strategy. We compute
#                          them honestly here using the strategy's own
#                          probabilistic structure - never by sharing a
#                          single cohort-level mean across groups.
#   (c) Score distribution - histogram across CCI bins, plus KS distance
#                          vs. gold and posterior coverage tables.
# =============================================================================

#' Per-encounter active gold groups (after dependency suppression).
#'
#' Wrapper around the internal helper so it has a stable export.
#' @export
gold_active_long <- function(dt, quan_map) {
  .gold_active_long_internal(dt, quan_map)
}

#' Build the JSON-coverage / mass-conservation table.
#'
#' Mass conservation is computed per Charlson group as
#'    mean_idx [ contribution_strategy(idx, gk) ]
#'  / mean_idx [ contribution_gold(idx, gk) ]
#' where the per-encounter contribution is `P(group active) * weight`
#' for the probabilistic strategies (S2 union, S3 imputation rate, S4
#' Dirichlet rate) and `1 * weight` if the gold algorithm activates the
#' group for that encounter.
#'
#' @param preds     data.table with at least `diagnosen`, `cci_gold`. If
#'   the optional per-group probabilities `s2_group_prob`, `s3_group_prob`,
#'   `s4_group_prob` (each a per-(idx, gk) data.table with column `p`)
#'   are stored as attributes of `preds`, they are used directly. Otherwise
#'   the function falls back to recomputing them from `cache` if supplied.
#' @param quan_map  output of `load_quan_map()`.
#' @param cache     optional precomputed lookup environment. If `preds`
#'   does not carry the strategy group_prob attributes and `cache` is NULL,
#'   only the gold side of the mass-conservation table is filled in.
#' @return data.table, one row per Charlson group, with columns:
#'   group, weight, n_patterns_in_json, n_gold_active, pct_gold_active,
#'   gold_mass, s2_mass, s3_mass, s4_mass, ratio_s2, ratio_s3, ratio_s4.
#' @export
qa_group_coverage <- function(preds, quan_map, cache = NULL) {
  gn <- names(quan_map)
  weights <- vapply(gn, function(g) as.numeric(quan_map[[g]]$weight), numeric(1L))
  n_pats  <- vapply(gn, function(g) length(.extract_patterns(quan_map[[g]])), integer(1L))
  N <- nrow(preds)

  active_gold <- .gold_active_long_internal(preds, quan_map)
  n_per_group <- active_gold[, .(n_active = uniqueN(idx)), by = gk]
  setnames(n_per_group, "gk", "group")

  cov <- data.table(group = gn,
                    weight = weights,
                    n_patterns_in_json = n_pats)
  cov <- merge(cov, n_per_group, by = "group", all.x = TRUE)
  cov[is.na(n_active), n_active := 0L]
  cov[, pct_gold_active := round(100 * n_active / N, 3)]
  cov[, gold_mass := n_active * weight / N]

  # Strategy mass: needs per-(idx, gk) probabilities.
  get_attr_dt <- function(attr_name) {
    a <- attr(preds, attr_name, exact = TRUE)
    if (is.null(a)) return(NULL)
    if (!is.data.table(a)) a <- as.data.table(a)
    a
  }
  s2_gp <- get_attr_dt("s2_group_prob")
  s3_gp <- get_attr_dt("s3_group_prob")
  s4_gp <- get_attr_dt("s4_group_prob")

  # Optional fallback: recompute S2 group probabilities if a cache is given.
  if (is.null(s2_gp) && !is.null(cache)) {
    res <- cci_probabilistic_batch(preds[, .(diagnosen)], quan_map, cache,
                                   return_group_prob = TRUE)
    s2_gp <- res$group_prob
  }

  per_group_mass <- function(gp_dt) {
    if (is.null(gp_dt) || nrow(gp_dt) == 0L) {
      return(setNames(rep(NA_real_, length(gn)), gn))
    }
    # mean over the cohort (treat missing (idx, gk) rows as p = 0)
    sums <- gp_dt[, .(s = sum(p)), by = gk]
    out <- setNames(rep(0, length(gn)), gn)
    out[sums$gk] <- sums$s / N
    # Multiply by weight to get mass.
    out * weights[names(out)]
  }

  s2m <- per_group_mass(s2_gp)
  s3m <- per_group_mass(s3_gp)
  s4m <- per_group_mass(s4_gp)

  cov[, s2_mass := s2m[group]]
  cov[, s3_mass := s3m[group]]
  cov[, s4_mass := s4m[group]]

  cov[, ratio_s2 := ifelse(gold_mass > 0, s2_mass / gold_mass, NA_real_)]
  cov[, ratio_s3 := ifelse(gold_mass > 0, s3_mass / gold_mass, NA_real_)]
  cov[, ratio_s4 := ifelse(gold_mass > 0, s4_mass / gold_mass, NA_real_)]

  setorder(cov, -pct_gold_active)
  cov[]
}

#' Score-distribution comparison: bin counts plus KS distance vs gold.
#' @export
qa_score_distribution <- function(preds, strategies = c("s2_ecci", "s3_mi",
                                                        "s4_bayes", "meta")) {
  bin_fn <- function(x) {
    cut(round(x), breaks = c(-Inf, 0, 1, 2, 3, 6, 10, Inf),
        labels = c("0", "1", "2", "3", "4-6", "7-10", "11+"),
        right = TRUE)
  }
  rows <- list()
  for (col in c("cci_gold", strategies)) {
    if (!col %in% names(preds)) next
    b <- bin_fn(preds[[col]])
    tab <- as.data.table(table(bin = b))
    tab[, source := col]
    tab[, pct := round(100 * N / sum(N), 3)]
    rows[[length(rows) + 1L]] <- tab
  }
  freq <- rbindlist(rows, use.names = TRUE)

  ks_rows <- list()
  for (col in strategies) {
    if (!col %in% names(preds)) next
    suppressWarnings({
      ks <- stats::ks.test(preds[[col]], preds$cci_gold)
    })
    ks_rows[[length(ks_rows) + 1L]] <- data.table(
      strategy  = col,
      ks_d      = round(unname(ks$statistic), 4),
      ks_p      = signif(ks$p.value, 3),
      mean_pred = round(mean(preds[[col]], na.rm = TRUE), 4),
      mean_gold = round(mean(preds$cci_gold, na.rm = TRUE), 4),
      mean_diff = round(mean(preds[[col]] - preds$cci_gold, na.rm = TRUE), 4)
    )
  }
  list(frequency = freq, ks = rbindlist(ks_rows))
}

#' Posterior-coverage check: share of encounters whose prediction is
#' within 0.5, 1, 2 CCI points of gold.
#' @export
qa_posterior_coverage <- function(preds,
                                  strategies = c("s2_ecci", "s3_mi",
                                                 "s4_bayes", "meta")) {
  rows <- list()
  for (col in strategies) {
    if (!col %in% names(preds)) next
    d <- abs(preds[[col]] - preds$cci_gold)
    rows[[length(rows) + 1L]] <- data.table(
      strategy        = col,
      pct_within_0_5  = round(100 * mean(d <= 0.5, na.rm = TRUE), 2),
      pct_within_1    = round(100 * mean(d <= 1.0, na.rm = TRUE), 2),
      pct_within_2    = round(100 * mean(d <= 2.0, na.rm = TRUE), 2),
      mean_abs_error  = round(mean(d, na.rm = TRUE), 4),
      max_abs_error   = round(max(d, na.rm = TRUE), 4)
    )
  }
  rbindlist(rows)
}
