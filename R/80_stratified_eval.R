# =============================================================================
# miCCI / 80_stratified_eval.R
# Per-stratum benchmark: which strategy wins where?
#
# Stratifications:
#   (1) ICD chapter (~22 strata) - dominant chapter of the encounter's codes.
#   (2) Active Charlson group   - sub-cohort: encounters where the listed
#                                  group is gold-active.
#   (3) Top-K multimorbid pattern - sorted tuple of gold-active groups.
#
# For every (stratum, strategy) we compute MAE / RMSE / R^2 / Bias plus a
# bootstrap 95% CI for MAE. The "winner" table picks, per stratum, the
# strategy with the smallest point MAE; if no other strategy's CI overlaps
# the winner's, the verdict is "use_<strategy>"; otherwise "tie".
# =============================================================================

#' ICD-10 chapter for a vector of three-character codes (vectorised via fcase).
#' @keywords internal
.icd_chapter <- function(code3) {
  if (length(code3) == 0L) return(character(0L))
  letter <- substr(code3, 1L, 1L)
  num    <- suppressWarnings(as.integer(substr(code3, 2L, 3L)))

  data.table::fcase(
    letter == "A" | letter == "B",                                 "I  Infectious",
    letter == "C",                                                 "II  Neoplasms",
    letter == "D" & !is.na(num) & num <= 48L,                      "II  Neoplasms",
    letter == "D" & !is.na(num) & num >= 49L,                      "III  Blood/Immune",
    letter == "E",                                                 "IV  Endocrine/Metab",
    letter == "F",                                                 "V  Mental",
    letter == "G",                                                 "VI  Nervous",
    letter == "H" & !is.na(num) & num <= 59L,                      "VII  Eye",
    letter == "H" & !is.na(num) & num >= 60L,                      "VIII  Ear",
    letter == "I",                                                 "IX  Circulatory",
    letter == "J",                                                 "X  Respiratory",
    letter == "K",                                                 "XI  Digestive",
    letter == "L",                                                 "XII  Skin",
    letter == "M",                                                 "XIII  Musculoskeletal",
    letter == "N",                                                 "XIV  Genitourinary",
    letter == "O",                                                 "XV  Pregnancy",
    letter == "P",                                                 "XVI  Perinatal",
    letter == "Q",                                                 "XVII  Congenital",
    letter == "R",                                                 "XVIII  Symptoms",
    letter == "S" | letter == "T",                                 "XIX  Injury",
    letter == "V" | letter == "W" | letter == "X" | letter == "Y", "XX  External",
    letter == "Z",                                                 "XXI  Health-status",
    letter == "U",                                                 "XXII  Special",
    default = "Unclassified"
  )
}

#' Dominant ICD chapter per encounter (chapter contributing the most
#' distinct three-character codes; ties broken alphabetically).
#' @keywords internal
.dominant_chapter <- function(diagnosen) {
  vapply(diagnosen, function(s) {
    codes <- unique(substr(toupper(gsub("[^A-Z0-9]", "",
                  unlist(strsplit(s, "\\|+")))), 1L, 3L))
    codes <- codes[nchar(codes) >= 3L]
    if (length(codes) == 0L) return(NA_character_)
    ch <- .icd_chapter(codes)
    tab <- sort(table(ch), decreasing = TRUE)
    names(tab)[1L]
  }, character(1L))
}

#' Per-encounter sorted active Charlson pattern (e.g. "dm_complicated|hf|kidney").
#' @keywords internal
.gold_pattern <- function(preds, quan_map) {
  active <- .gold_active_long_internal(preds, quan_map)
  if (nrow(active) == 0L) return(rep("(none)", nrow(preds)))
  pat <- active[order(idx, gk),
                .(pattern = paste(sort(unique(gk)), collapse = "|")),
                by = idx]
  out <- rep("(none)", nrow(preds))
  out[pat$idx] <- pat$pattern
  out
}

#' Per-stratum metrics with bootstrap CIs for every strategy.
#'
#' @param preds        data.table from `compute_predictions()`.
#' @param strategies   named character vector (column = pretty label).
#' @param quan_map     output of `load_quan_map()`.
#' @param B            bootstrap resamples per stratum.
#' @param min_n        minimum encounter count for a stratum to be reported.
#' @param top_patterns number of most-frequent patterns to evaluate.
#' @param parallel     if TRUE and `future.apply` is available, parallelise
#'   the bootstrap loop.
#' @return list of three data.tables: chapter, group, pattern.
#' @export
stratified_evaluation <- function(preds, strategies, quan_map,
                                  B = 500L, min_n = 50L,
                                  top_patterns = 30L,
                                  cluster_col = NULL,
                                  parallel = FALSE) {
  message("  computing chapter labels")
  preds <- copy(preds)
  preds[, .chapter := .dominant_chapter(diagnosen)]
  message("  computing gold-active groups + patterns")
  active <- .gold_active_long_internal(preds, quan_map)
  preds[, .pattern := .gold_pattern(.SD, quan_map)]

  # Row membership per stratum, carried as an attribute so the winner test can
  # rerun a paired bootstrap on the same encounters without recomputing the
  # stratification.
  membership <- list(chapter = list(), group = list(), pattern = list())

  empty_strata_dt <- function() {
    data.table(
      stratum  = character(0L),
      strategy = character(0L), label = character(0L),
      n       = integer(0L),
      mae     = numeric(0L), mae_lo  = numeric(0L), mae_hi  = numeric(0L),
      rmse    = numeric(0L), rmse_lo = numeric(0L), rmse_hi = numeric(0L),
      r2      = numeric(0L), r2_lo   = numeric(0L), r2_hi   = numeric(0L),
      bias    = numeric(0L), bias_lo = numeric(0L), bias_hi = numeric(0L)
    )
  }
  safe_set_stratum_first <- function(dt) {
    if (!is.data.table(dt) || nrow(dt) == 0L) return(dt)
    if ("stratum" %in% names(dt))
      setcolorder(dt, c("stratum", setdiff(names(dt), "stratum")))
    dt
  }
  do_strata <- function(label_col, sub_dt, slot_name) {
    if (!label_col %in% names(sub_dt)) return(empty_strata_dt())
    keys <- sub_dt[, .N, by = c(label_col)][N >= min_n][[label_col]]
    if (length(keys) == 0L) return(empty_strata_dt())
    out_rows <- list()
    for (k in keys) {
      rows <- which(sub_dt[[label_col]] == k)
      if (length(rows) < min_n) next
      sub <- sub_dt[rows]
      r <- bootstrap_strategies(sub, strategies, B = B, seed = 42L,
                                cluster_col = cluster_col,
                                parallel = parallel)
      if (!is.data.table(r) || nrow(r) == 0L) next
      r[, stratum := k]
      out_rows[[length(out_rows) + 1L]] <- r
      membership[[slot_name]][[as.character(k)]] <<- rows
    }
    if (length(out_rows) == 0L) return(empty_strata_dt())
    safe_set_stratum_first(rbindlist(out_rows, use.names = TRUE, fill = TRUE))
  }

  message("  evaluating ICD chapters")
  chap <- do_strata(".chapter", preds, "chapter")

  message("  evaluating per-Charlson-group strata")
  grp_rows <- list()
  if (nrow(active) > 0L) {
    for (g in unique(active$gk)) {
      idxs <- active[gk == g, unique(idx)]
      if (length(idxs) < min_n) next
      sub <- preds[idxs]
      r <- bootstrap_strategies(sub, strategies, B = B, seed = 42L,
                                cluster_col = cluster_col,
                                parallel = parallel)
      if (!is.data.table(r) || nrow(r) == 0L) next
      r[, stratum := g]
      grp_rows[[length(grp_rows) + 1L]] <- r
      membership$group[[as.character(g)]] <- idxs
    }
  }
  grp <- if (length(grp_rows) > 0L)
           rbindlist(grp_rows, use.names = TRUE, fill = TRUE)
         else empty_strata_dt()
  grp <- safe_set_stratum_first(grp)

  message(sprintf("  evaluating top-%d multimorbid patterns", top_patterns))
  pat_freq <- preds[, .N, by = .pattern][order(-N)]
  top <- head(pat_freq[.pattern != "(none)" & N >= min_n], top_patterns)
  pat_rows <- list()
  if (nrow(top) > 0L) {
    for (p in top$.pattern) {
      rows <- which(preds$.pattern == p)
      if (length(rows) < min_n) next
      sub <- preds[rows]
      r <- bootstrap_strategies(sub, strategies, B = B, seed = 42L,
                                cluster_col = cluster_col,
                                parallel = parallel)
      if (!is.data.table(r) || nrow(r) == 0L) next
      r[, stratum := p]
      pat_rows[[length(pat_rows) + 1L]] <- r
      membership$pattern[[as.character(p)]] <- rows
    }
  }
  pat <- if (length(pat_rows) > 0L)
           rbindlist(pat_rows, use.names = TRUE, fill = TRUE)
         else empty_strata_dt()
  pat <- safe_set_stratum_first(pat)

  setattr(chap, "membership", membership$chapter)
  setattr(grp,  "membership", membership$group)
  setattr(pat,  "membership", membership$pattern)
  list(chapter = chap, group = grp, pattern = pat)
}

#' Per-stratum winner with a paired-bootstrap verdict.
#'
#' For each stratum the lowest-MAE strategy is identified, then tested against
#' the runner-up.
#'
#' Two verdict methods exist and the one used is recorded in the `method`
#' column, so a downstream table can never be ambiguous about which was run.
#'
#' \describe{
#'   \item{`paired`}{the default whenever `preds` and `strategies` are
#'     supplied. Both strategies are recomputed on the same bootstrap
#'     replicate and the interval is for the MAE difference. This is the
#'     correct test, because the two strategies are scored on identical
#'     encounters and their errors are strongly correlated.}
#'   \item{`ci_overlap`}{the fallback when the raw predictions are not
#'     available. It asks only whether the marginal interval of the winner
#'     overlaps any competitor. It is a weaker, more conservative screen and
#'     must not be described as a paired test.}
#' }
#'
#' @param stratum_dt per-stratum metric table from `stratified_evaluation()`.
#' @param preds optional predictions table, required for the paired method.
#' @param strategies optional named strategy vector, required for the paired
#'   method.
#' @param membership optional named list of row indices per stratum. Defaults
#'   to the `membership` attribute that `stratified_evaluation()` attaches.
#' @param B paired-bootstrap replicates.
#' @param conf confidence level for the paired interval. Pass
#'   `1 - alpha / K` to apply a Bonferroni correction across the K strata
#'   tested in the same family.
#' @param cluster_col optional patient-id column in `preds` for a
#'   cluster-robust paired bootstrap.
#' @param seed RNG seed (local).
#' @export
stratum_winners <- function(stratum_dt,
                            preds = NULL, strategies = NULL,
                            membership = NULL,
                            B = 1000L, conf = 0.95,
                            cluster_col = NULL, seed = 42L) {
  empty_winners <- function() data.table(
    stratum       = character(0L),
    n_strategies  = integer(0L),
    n_encounters  = integer(0L),
    best_strategy = character(0L),
    best_mae      = numeric(0L),
    best_mae_lo   = numeric(0L),
    best_mae_hi   = numeric(0L),
    runner_up     = character(0L),
    runner_mae    = numeric(0L),
    diff          = numeric(0L),
    diff_lo       = numeric(0L),
    diff_hi       = numeric(0L),
    conf          = numeric(0L),
    method        = character(0L),
    verdict       = character(0L)
  )
  if (!is.data.table(stratum_dt) || nrow(stratum_dt) == 0L) return(empty_winners())
  if (!"stratum" %in% names(stratum_dt)) return(empty_winners())

  if (is.null(membership)) membership <- attr(stratum_dt, "membership", exact = TRUE)
  use_paired <- !is.null(preds) && !is.null(strategies) && !is.null(membership)
  col_of <- if (!is.null(strategies))
    setNames(names(strategies), unname(strategies)) else character(0L)
  cl_all <- if (use_paired && !is.null(cluster_col) && cluster_col %in% names(preds))
    preds[[cluster_col]] else NULL

  out <- stratum_dt[, {
    ord <- order(mae)
    bi  <- ord[1L]
    best_lab <- label[bi]
    if (.N >= 2L) {
      runner_lab <- label[ord[2L]]
      runner_val <- round(mae[ord[2L]], 4L)
    } else {
      runner_lab <- NA_character_; runner_val <- NA_real_
    }

    d <- NA_real_; d_lo <- NA_real_; d_hi <- NA_real_
    meth <- "ci_overlap"; dec <- NA

    if (use_paired && .N >= 2L) {
      rows <- membership[[as.character(.BY$stratum)]]
      ca <- col_of[[best_lab]]; cb <- col_of[[runner_lab]]
      if (!is.null(rows) && length(rows) > 1L &&
          !is.null(ca) && !is.null(cb)) {
        pb <- paired_bootstrap_diff(
          preds[[ca]][rows], preds[[cb]][rows], preds$cci_gold[rows],
          B = B, seed = seed, conf = conf,
          cluster = if (is.null(cl_all)) NULL else cl_all[rows])
        d <- pb$diff; d_lo <- pb$diff_lo; d_hi <- pb$diff_hi
        dec <- isTRUE(pb$decisive); meth <- "paired"
      }
    }
    if (identical(meth, "ci_overlap")) {
      dec <- if (.N >= 2L)
        !any(mae_lo[-bi] <= mae_hi[bi], na.rm = TRUE) else FALSE
    }

    .(
      n_strategies  = .N,
      n_encounters  = n[bi],
      best_strategy = best_lab,
      best_mae      = round(mae[bi], 4L),
      best_mae_lo   = round(mae_lo[bi], 4L),
      best_mae_hi   = round(mae_hi[bi], 4L),
      runner_up     = runner_lab,
      runner_mae    = runner_val,
      diff          = round(d, 5L),
      diff_lo       = round(d_lo, 5L),
      diff_hi       = round(d_hi, 5L),
      conf          = conf,
      method        = meth,
      verdict       = if (isTRUE(dec)) paste0("use_", best_lab) else "tie"
    )
  }, by = stratum]

  if (nrow(out) == 0L) return(empty_winners())
  setorder(out, best_mae)
  out[]
}
