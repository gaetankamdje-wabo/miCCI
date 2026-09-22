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
  s1_gp <- get_attr_dt("s1_group_prob")
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

  s1m <- per_group_mass(s1_gp)
  s2m <- per_group_mass(s2_gp)
  s3m <- per_group_mass(s3_gp)
  s4m <- per_group_mass(s4_gp)

  cov[, s1_mass := s1m[group]]
  cov[, s2_mass := s2m[group]]
  cov[, s3_mass := s3m[group]]
  cov[, s4_mass := s4m[group]]

  cov[, ratio_s1 := ifelse(gold_mass > 0, s1_mass / gold_mass, NA_real_)]
  cov[, ratio_s2 := ifelse(gold_mass > 0, s2_mass / gold_mass, NA_real_)]
  cov[, ratio_s3 := ifelse(gold_mass > 0, s3_mass / gold_mass, NA_real_)]
  cov[, ratio_s4 := ifelse(gold_mass > 0, s4_mass / gold_mass, NA_real_)]

  setorder(cov, -pct_gold_active)
  cov[]
}

#' Exact-agreement table: no rounding, no binning.
#'
#' Binning a reconstructed score forces a rounding rule onto estimators whose
#' natural support differs. S1 is half-integer, S2, S3 and the Meta Learner are
#' real-valued, S4 and the reference are integer-valued. Under R's default
#' round-half-to-even, S1's half-integers land disproportionately on even bins,
#' so a histogram built that way shows a distributional pattern that is really
#' an artefact of the rounding convention.
#'
#' This function sidesteps the question entirely and reports what the
#' estimators actually produce:
#'
#' \describe{
#'   \item{`by_value`}{per integer reference score v: how many encounters have
#'     gold exactly v, how many have the strategy exactly v, and how many have
#'     BOTH exactly v (the exact hit).}
#'   \item{`summary`}{per strategy: the share of predictions that are
#'     integer-valued at all, the share that are exact half-integers, and the
#'     overall exact-agreement rate against gold.}
#' }
#'
#' Reported alongside the binned histogram, this separates genuine
#' distributional behaviour from the rounding rule.
#'
#' @param preds predictions table including `cci_gold`.
#' @param strategies strategy column names.
#' @param max_value highest reference score reported individually; everything
#'   above is pooled into one row so the table stays readable.
#' @export
qa_exact_agreement <- function(preds,
                               strategies = c("s1_mid", "s2_ecci", "s3_mi",
                                              "s4_bayes", "meta"),
                               max_value = 10L) {
  N <- nrow(preds)
  gold <- preds$cci_gold
  vals <- 0:max_value

  by_rows <- list(); sum_rows <- list()
  for (col in strategies) {
    if (!col %in% names(preds)) next
    p <- preds[[col]]
    int_p <- .is_int_valued(p)

    for (v in vals) {
      g_is <- gold == v
      p_is <- int_p & (round(p) == v)
      by_rows[[length(by_rows) + 1L]] <- data.table(
        strategy       = col,
        cci_value      = as.character(v),
        n_gold         = sum(g_is, na.rm = TRUE),
        n_pred_exact   = sum(p_is, na.rm = TRUE),
        n_exact_match  = sum(g_is & p_is, na.rm = TRUE)
      )
    }
    g_hi <- gold > max_value
    p_hi <- int_p & (round(p) > max_value)
    by_rows[[length(by_rows) + 1L]] <- data.table(
      strategy       = col,
      cci_value      = paste0(">", max_value),
      n_gold         = sum(g_hi, na.rm = TRUE),
      n_pred_exact   = sum(p_hi, na.rm = TRUE),
      n_exact_match  = sum(g_hi & p_hi, na.rm = TRUE)
    )

    half <- abs(p - floor(p) - 0.5) < 1e-8
    sum_rows[[length(sum_rows) + 1L]] <- data.table(
      strategy            = col,
      n                   = N,
      pct_integer_valued  = round(100 * mean(int_p, na.rm = TRUE), 2),
      pct_half_integer    = round(100 * mean(half, na.rm = TRUE), 2),
      n_exact_agreement   = sum(int_p & (round(p) == gold), na.rm = TRUE),
      pct_exact_agreement = round(100 * mean(int_p & (round(p) == gold),
                                             na.rm = TRUE), 2)
    )
  }

  by_value <- rbindlist(by_rows, use.names = TRUE)
  if (nrow(by_value) > 0L) {
    by_value[, pct_of_gold := ifelse(n_gold > 0,
                                     round(100 * n_exact_match / n_gold, 2),
                                     NA_real_)]
    by_value[, pct_of_cohort := round(100 * n_exact_match / N, 3)]
  }
  list(by_value = by_value, summary = rbindlist(sum_rows, use.names = TRUE))
}

#' Score-distribution comparison: bin counts plus KS distance vs gold.
#'
#' Binning uses `.round_half_up()` rather than `base::round()`, so half-integer
#' predictions are not pushed preferentially onto even bins. Read this table
#' together with `qa_exact_agreement()`, which needs no rounding at all.
#' @export
qa_score_distribution <- function(preds, strategies = c("s2_ecci", "s3_mi",
                                                        "s4_bayes", "meta")) {
  bin_fn <- function(x) {
    cut(.round_half_up(x), breaks = c(-Inf, 0, 1, 2, 3, 6, 10, Inf),
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