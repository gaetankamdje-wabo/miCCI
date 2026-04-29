# =============================================================================
# miCCI / 20_cci_probabilistic.R
# S2 - Probabilistic E-CCI.
#
# For an encounter with truncated prefixes p_1, ..., p_K, the probability
# that Charlson group g is active is:
#
#   P(g active) = 1 - prod_{i=1..K} (1 - q_{i,g})
#
# where q_{i,g} = sum over children c of prefix p_i of P(c | p_i) * 1{c maps to g}.
# This is the inclusion-exclusion union under the assumption that the
# realisations of distinct truncated prefixes are independent. It replaces
# the v0.5.0 max(...) heuristic, which under-attributed mass when several
# different prefixes for the same encounter could each populate the same
# Charlson group.
#
# The expected CCI is then
#
#   E[CCI] = sum_g P(g active) * w_g
#
# with the dependency rule applied as in the gold algorithm, but in
# expectation: a child group's contribution is weighted by the probability
# that none of its parents is active, i.e.
#
#   E[CCI] = sum_g P(g active) * prod_{p in parents(g)} (1 - P(p active)) * w_g.
# =============================================================================

#' Single-encounter S2.
#' @param icd_anon character vector of truncated ICD codes.
#' @param quan_map output of `load_quan_map()`.
#' @param cache    output of `precompute_lookups()`.
#' @param age      optional numeric age for age-stratified subcode probs.
#' @return list with elements `e_cci` and `group_prob` (named numeric).
#' @export
cci_probabilistic <- function(icd_anon, quan_map, cache, age = NULL) {
  codes <- unique(icd3(icd_anon))
  gn <- get(".group_names", envir = cache)
  dm <- get(".dep_map",     envir = cache)
  wm <- get(".wt_map",      envir = cache)

  # Per-prefix per-group hit probability q_{i,g}.
  # Accumulate the union (1 - prod(1 - q)) across prefixes for each group.
  log_one_minus <- setNames(numeric(length(gn)), gn)  # log(prod(1 - q_i))
  for (code in codes) {
    pc <- get_prefix_cache(code, cache)
    if (is.null(pc)) next
    probs <- if (!is.null(age) && !is.na(age))
      get_subcode_probs_fast(code, cache, age = age)
    else
      pc$child_freqs[, .(code_nodot, prob)]
    if (nrow(probs) == 0L) next

    cg <- pc$child_groups
    for (gk in gn) {
      q <- sum(probs$prob[vapply(probs$code_nodot,
                                 function(cn) gk %in% (cg[[cn]] %||% character(0L)),
                                 logical(1L))],
               na.rm = TRUE)
      if (q > 0) log_one_minus[gk] <- log_one_minus[gk] + log1p(-min(q, 1 - 1e-12))
    }
  }
  group_prob <- 1 - exp(log_one_minus)
  group_prob[group_prob < 0] <- 0

  # Apply dependency rule in expectation: a child contributes
  # P(child) * prod_{p in parents}(1 - P(parent)) * w_child.
  cont <- numeric(length(gn)); names(cont) <- gn
  for (gk in gn) {
    if (group_prob[gk] <= 0) next
    deps <- dm[[gk]]
    factor <- 1
    if (length(deps) > 0L) {
      pp <- group_prob[deps]
      pp <- pp[!is.na(pp)]
      if (length(pp) > 0L) factor <- prod(1 - pp)
    }
    cont[gk] <- group_prob[gk] * factor * wm[gk]
  }
  list(e_cci = sum(cont), group_prob = group_prob)
}

#' Vectorised batch S2.
#'
#' Returns the per-encounter expected CCI as a numeric vector. If
#' `return_group_prob = TRUE`, returns a list with both the e_cci vector
#' and the per-(idx, gk) probability table - this is what the QA stage
#' needs for honest mass-conservation checks.
#'
#' @export
cci_probabilistic_batch <- function(dt, quan_map, cache,
                                    return_group_prob = FALSE) {
  gn <- get(".group_names", envir = cache)
  wm <- get(".wt_map",      envir = cache)
  dl <- get(".dep_lookup_dt", envir = cache)

  # Per-prefix per-group hit probability q_{pref, gk}.
  all_pref <- ls(cache)
  all_pref <- all_pref[!startsWith(all_pref, ".")]
  pg_rows <- list()
  for (pref in all_pref) {
    pc  <- get(pref, envir = cache, inherits = FALSE)
    cdt <- pc$child_freqs
    cg  <- pc$child_groups
    for (gk in gn) {
      q <- sum(cdt$prob[vapply(cdt$code_nodot,
                               function(cn) gk %in% (cg[[cn]] %||% character(0L)),
                               logical(1L))],
               na.rm = TRUE)
      if (q > 0)
        pg_rows[[length(pg_rows) + 1L]] <- list(pref = pref, gk = gk, q = q)
    }
  }
  pgp <- if (length(pg_rows) > 0L) rbindlist(pg_rows)
         else data.table(pref = character(0L), gk = character(0L), q = numeric(0L))

  # Long form (idx, pref).
  dx_list <- strsplit(as.character(dt$diagnosen), "\\|+")
  long <- data.table(
    idx  = rep.int(seq_len(nrow(dt)), lengths(dx_list)),
    code = toupper(gsub("[^A-Z0-9]", "", unlist(dx_list, use.names = FALSE)))
  )
  long <- long[nchar(code) >= 3L]
  long[, pref := substr(code, 1L, 3L)]
  long <- unique(long[, .(idx, pref)])

  hits <- merge(long, pgp, by = "pref", allow.cartesian = TRUE)
  if (nrow(hits) == 0L) {
    out <- numeric(nrow(dt))
    if (return_group_prob)
      return(list(e_cci = out,
                  group_prob = data.table(idx = integer(0L),
                                          gk = character(0L),
                                          p = numeric(0L))))
    return(out)
  }

  # Inclusion-exclusion: P(group g) = 1 - prod_i (1 - q_i).
  # Numerically stable via sum of log1p(-q).
  hits[, q_clip := pmin(q, 1 - 1e-12)]
  enc_gp <- hits[, .(p = 1 - exp(sum(log1p(-q_clip)))), by = .(idx, gk)]
  enc_gp[p < 0, p := 0]

  # Dependency in expectation: weight by prod(1 - P(parent)).
  if (nrow(dl) > 0L && nrow(enc_gp) > 0L) {
    parents <- merge(enc_gp[, .(idx, child_gk = gk)], dl,
                     by.x = "child_gk", by.y = "child", allow.cartesian = TRUE)
    parents <- merge(parents, enc_gp,
                     by.x = c("idx", "parent"), by.y = c("idx", "gk"),
                     all.x = TRUE)
    setnames(parents, "p", "p_parent")
    parents[is.na(p_parent), p_parent := 0]
    factor_dt <- parents[, .(factor = prod(1 - p_parent)),
                         by = .(idx, gk = child_gk)]
    enc_gp <- merge(enc_gp, factor_dt, by = c("idx", "gk"), all.x = TRUE)
    enc_gp[is.na(factor), factor := 1]
  } else {
    enc_gp[, factor := 1]
  }

  wt_dt  <- data.table(gk = names(wm), w = unname(wm))
  scored <- merge(enc_gp, wt_dt, by = "gk")
  scored[, contribution := p * factor * w]
  ecci <- scored[, .(e_cci = sum(contribution)), by = idx]

  out <- numeric(nrow(dt))
  out[ecci$idx] <- ecci$e_cci

  if (return_group_prob) {
    return(list(e_cci = out,
                group_prob = enc_gp[, .(idx, gk, p)]))
  }
  out
}
