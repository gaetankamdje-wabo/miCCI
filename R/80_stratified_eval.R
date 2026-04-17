# =============================================================================
# miCCI v1.0.0 — 80_stratified_eval.R
# Per-stratum benchmark: which strategy wins where?
#
# Three stratifications are produced:
#
#   (1) ICD chapter (22 strata)            — based on the dominant chapter
#                                              of the encounter's diagnoses.
#   (2) Active Charlson group (17 strata)  — sub-cohort: encounters in which
#                                              that single group is gold-active.
#                                              An encounter can appear in
#                                              several group strata.
#   (3) Multimorbid pattern (top 30)       — the sorted tuple of gold-active
#                                              Charlson groups, restricted to
#                                              the 30 most frequent patterns.
#
# For every (stratum, strategy) we compute MAE / RMSE / R^2 / Bias plus a
# bootstrap 95% CI for MAE. The "winner" table contains, per stratum, the
# strategy whose MAE upper-CI is the smallest, with a categorical
# recommendation:
#   * "use_<strategy>" — if its MAE 95% CI does not overlap with any other
#                        strategy's MAE 95% CI for that stratum.
#   * "tie"            — if at least one other strategy's CI overlaps.
# =============================================================================

#' ICD-10 chapter classification (22 chapters, A00-Z99 + U).
#' Returns chapter label for a single 3-char code.
#' @keywords internal
.icd_chapter <- function(code3) {
  if (length(code3) == 0) return(character(0))
  letter <- substr(code3, 1, 1)
  num    <- suppressWarnings(as.integer(substr(code3, 2, 3)))
  out <- character(length(code3))

  rules <- list(
    list(let = "A", lo =  0, hi = 99, lab = "I  Infectious"),
    list(let = "B", lo =  0, hi = 99, lab = "I  Infectious"),
    list(let = "C", lo =  0, hi = 97, lab = "II  Neoplasms"),
    list(let = "D", lo =  0, hi = 48, lab = "II  Neoplasms"),
    list(let = "D", lo = 49, hi = 89, lab = "III  Blood/Immune"),
    list(let = "E", lo =  0, hi = 90, lab = "IV  Endocrine/Metab"),
    list(let = "F", lo =  0, hi = 99, lab = "V  Mental"),
    list(let = "G", lo =  0, hi = 99, lab = "VI  Nervous"),
    list(let = "H", lo =  0, hi = 59, lab = "VII  Eye"),
    list(let = "H", lo = 60, hi = 95, lab = "VIII  Ear"),
    list(let = "I", lo =  0, hi = 99, lab = "IX  Circulatory"),
    list(let = "J", lo =  0, hi = 99, lab = "X  Respiratory"),
    list(let = "K", lo =  0, hi = 93, lab = "XI  Digestive"),
    list(let = "L", lo =  0, hi = 99, lab = "XII  Skin"),
    list(let = "M", lo =  0, hi = 99, lab = "XIII  Musculoskeletal"),
    list(let = "N", lo =  0, hi = 99, lab = "XIV  Genitourinary"),
    list(let = "O", lo =  0, hi = 99, lab = "XV  Pregnancy"),
    list(let = "P", lo =  0, hi = 96, lab = "XVI  Perinatal"),
    list(let = "Q", lo =  0, hi = 99, lab = "XVII  Congenital"),
    list(let = "R", lo =  0, hi = 99, lab = "XVIII  Symptoms"),
    list(let = "S", lo =  0, hi = 99, lab = "XIX  Injury"),
    list(let = "T", lo =  0, hi = 98, lab = "XIX  Injury"),
    list(let = "V", lo =  0, hi = 99, lab = "XX  External"),
    list(let = "W", lo =  0, hi = 99, lab = "XX  External"),
    list(let = "X", lo =  0, hi = 99, lab = "XX  External"),
    list(let = "Y", lo =  0, hi = 99, lab = "XX  External"),
    list(let = "Z", lo =  0, hi = 99, lab = "XXI  Health-status"),
    list(let = "U", lo =  0, hi = 99, lab = "XXII  Special")
  )
  out[] <- "Unclassified"
  for (r in rules) {
    sel <- !is.na(num) & letter == r$let & num >= r$lo & num <= r$hi
    out[sel] <- r$lab
  }
  out
}

#' Compute the dominant ICD chapter per encounter (the chapter contributing
#' the most distinct codes; ties broken alphabetically).
#' @keywords internal
.dominant_chapter <- function(diagnosen) {
  vapply(diagnosen, function(s) {
    codes <- unique(substr(toupper(gsub("[^A-Z0-9]", "",
                  unlist(strsplit(s, "\\|+")))), 1, 3))
    codes <- codes[nchar(codes) >= 3]
    if (length(codes) == 0) return(NA_character_)
    ch <- .icd_chapter(codes)
    tab <- sort(table(ch), decreasing = TRUE)
    names(tab)[1]
  }, character(1))
}

#' Per-encounter sorted active Charlson pattern, e.g. "dm_complicated|hf|kidney".
#' @keywords internal
.gold_pattern <- function(preds, quan_map) {
  active <- .gold_active_long(preds, quan_map)
  if (nrow(active) == 0)
    return(rep("(none)", nrow(preds)))
  pat <- active[order(idx, gk),
                .(pattern = paste(sort(unique(gk)), collapse = "|")),
                by = idx]
  out <- rep("(none)", nrow(preds))
  out[pat$idx] <- pat$pattern
  out
}

#' Compute per-stratum metrics + bootstrap CIs for every strategy.
#'
#' @param preds       data.table from compute_predictions()
#' @param strategies  named character vector — see bootstrap_strategies()
#' @param quan_map    Quan map (needed for Charlson-group strata)
#' @param B           bootstrap resamples per stratum (default 500)
#' @param min_n       minimum encounter count to include a stratum
#' @param top_patterns number of most-frequent patterns to evaluate
#'
#' @return list of three data.tables: chapter, group, pattern.
#'         Each row: stratum, n, strategy, label, mae, mae_lo, mae_hi,
#'                   rmse, rmse_lo, rmse_hi, r2, r2_lo, r2_hi,
#'                   bias, bias_lo, bias_hi
#' @export
stratified_evaluation <- function(preds, strategies, quan_map,
                                  B = 500L, min_n = 50L,
                                  top_patterns = 30L,
                                  parallel = FALSE) {
  cat("  computing chapter labels...\n")
  preds <- copy(preds)
  preds[, .chapter := .dominant_chapter(diagnosen)]
  cat("  computing gold-active groups + patterns...\n")
  active <- .gold_active_long(preds, quan_map)
  preds[, .pattern := .gold_pattern(.SD, quan_map)]

  # ---- empty-result template ------------------------------------------------
  # bootstrap_strategies() returns a data.table with this column set. When no
  # stratum produces any rows we still return an empty data.table that has
  # exactly these columns so that downstream code (setcolorder, fwrite,
  # stratum_winners) never crashes on a missing column.
  empty_strata_dt <- function() {
    data.table(
      stratum  = character(0),
      strategy = character(0), label = character(0),
      n       = integer(0),
      mae     = numeric(0), mae_lo  = numeric(0), mae_hi  = numeric(0),
      rmse    = numeric(0), rmse_lo = numeric(0), rmse_hi = numeric(0),
      r2      = numeric(0), r2_lo   = numeric(0), r2_hi   = numeric(0),
      bias    = numeric(0), bias_lo = numeric(0), bias_hi = numeric(0)
    )
  }

  safe_set_stratum_first <- function(dt) {
    if (!is.data.table(dt) || nrow(dt) == 0) return(dt)
    if ("stratum" %in% names(dt))
      setcolorder(dt, c("stratum", setdiff(names(dt), "stratum")))
    dt
  }

  do_strata <- function(label_col, sub_dt) {
    if (!label_col %in% names(sub_dt)) return(empty_strata_dt())
    keys <- sub_dt[, .N, by = c(label_col)][N >= min_n][[label_col]]
    if (length(keys) == 0) return(empty_strata_dt())
    out_rows <- list()
    for (k in keys) {
      sub <- sub_dt[get(label_col) == k]
      if (nrow(sub) < min_n) next
      r <- bootstrap_strategies(sub, strategies, B = B, seed = 42L,
                                parallel = parallel)
      if (!is.data.table(r) || nrow(r) == 0) next
      r[, stratum := k]
      out_rows[[length(out_rows) + 1L]] <- r
    }
    if (length(out_rows) == 0) return(empty_strata_dt())
    res <- rbindlist(out_rows, use.names = TRUE, fill = TRUE)
    safe_set_stratum_first(res)
  }

  # (1) chapter
  cat("  evaluating ICD chapters...\n")
  chap <- do_strata(".chapter", preds)

  # (2) per Charlson group: encounter is in stratum if group is gold-active
  cat("  evaluating per-Charlson-group strata...\n")
  grp_rows <- list()
  if (nrow(active) > 0) {
    for (g in unique(active$gk)) {
      idxs <- active[gk == g, unique(idx)]
      if (length(idxs) < min_n) next
      sub <- preds[idxs]
      r <- bootstrap_strategies(sub, strategies, B = B, seed = 42L,
                                parallel = parallel)
      if (!is.data.table(r) || nrow(r) == 0) next
      r[, stratum := g]
      grp_rows[[length(grp_rows) + 1L]] <- r
    }
  }
  grp <- if (length(grp_rows) > 0)
           rbindlist(grp_rows, use.names = TRUE, fill = TRUE)
         else empty_strata_dt()
  grp <- safe_set_stratum_first(grp)

  # (3) top multimorbid patterns
  cat(sprintf("  evaluating top-%d multimorbid patterns...\n", top_patterns))
  pat_freq <- preds[, .N, by = .pattern][order(-N)]
  top <- head(pat_freq[.pattern != "(none)" & N >= min_n], top_patterns)
  pat_rows <- list()
  if (nrow(top) > 0) {
    for (p in top$.pattern) {
      sub <- preds[.pattern == p]
      if (nrow(sub) < min_n) next
      r <- bootstrap_strategies(sub, strategies, B = B, seed = 42L,
                                parallel = parallel)
      if (!is.data.table(r) || nrow(r) == 0) next
      r[, stratum := p]
      pat_rows[[length(pat_rows) + 1L]] <- r
    }
  }
  pat <- if (length(pat_rows) > 0)
           rbindlist(pat_rows, use.names = TRUE, fill = TRUE)
         else empty_strata_dt()
  pat <- safe_set_stratum_first(pat)

  list(chapter = chap, group = grp, pattern = pat)
}

#' For each stratum, find the strategy with the lowest MAE upper-CI and
#' tag whether it is significantly better than the competitors.
#' @export
stratum_winners <- function(stratum_dt) {
  empty_winners <- function() data.table(
    stratum       = character(0),
    n_strategies  = integer(0),
    n_encounters  = integer(0),
    best_strategy = character(0),
    best_mae      = numeric(0),
    best_mae_lo   = numeric(0),
    best_mae_hi   = numeric(0),
    verdict       = character(0),
    runner_up     = character(0),
    runner_mae    = numeric(0)
  )

  if (!is.data.table(stratum_dt) || nrow(stratum_dt) == 0) return(empty_winners())
  if (!"stratum" %in% names(stratum_dt)) return(empty_winners())

  out <- stratum_dt[, {
    ord <- order(mae)                       # rank by point MAE
    best_idx <- ord[1]
    best_lab <- label[best_idx]
    best_mae_hi <- mae_hi[best_idx]
    if (.N >= 2) {
      others_lo <- mae_lo[-best_idx]
      overlap   <- any(others_lo <= best_mae_hi, na.rm = TRUE)
      runner_lab <- label[ord[2]]
      runner_val <- round(mae[ord[2]], 4)
    } else {
      overlap    <- TRUE
      runner_lab <- NA_character_
      runner_val <- NA_real_
    }
    .(
      n_strategies  = .N,
      n_encounters  = n[best_idx],
      best_strategy = best_lab,
      best_mae      = round(mae[best_idx], 4),
      best_mae_lo   = round(mae_lo[best_idx], 4),
      best_mae_hi   = round(best_mae_hi, 4),
      verdict       = ifelse(overlap, "tie", paste0("use_", best_lab)),
      runner_up     = runner_lab,
      runner_mae    = runner_val
    )
  }, by = stratum]

  if (nrow(out) == 0) return(empty_winners())
  setorder(out, best_mae)
  out[]
}
