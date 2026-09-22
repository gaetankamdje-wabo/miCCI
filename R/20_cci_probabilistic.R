# =============================================================================
# miCCI / 20_cci_probabilistic.R
# S2 - Probabilistic E-CCI.
#
# For an encounter with truncated prefixes p_1, ..., p_K, the probability
# that Charlson group g is active is:
#
#   P(g active) = 1 - prod_{i=1..K} (1 - q_{i,g})^{k_i}
#
# where q_{i,g} = sum over children c of prefix p_i of P(c | p_i) * 1{c maps to
# g}, and k_i is how often prefix p_i was coded for the encounter. This is the
# inclusion-exclusion union under the assumption that the realisations of
# distinct coded diagnosis positions are independent.
#
# HIERARCHICAL EXCLUSION (corrected in this revision)
# --------------------------------------------
# The Quan rule says a child group does not score when its severity-coded
# parent is also active. Previously that was applied as
#
#   contribution(c) = P(c) * prod_{p in parents(c)} (1 - P(p)) * w_c,
#
# which treats "child active" and "parent active" as independent events. They
# are not. When both branches live under the same truncated prefix, a single
# drawn subcode lands in at most one of them, so the events are mutually
# exclusive within that prefix and the product form understates the joint
# suppression.
#
# This revision computes the quantity exactly. Per prefix p define
#
#   q_par[p] = P(drawn subcode of p maps to any parent of c)
#   q_uni[p] = P(drawn subcode of p maps to c OR to any parent of c)
#
# Writing A_p = "subcode of p maps to c" and B_p = "subcode of p maps to a
# parent of c", and keeping only the across-prefix independence assumption:
#
#   P(c active AND no parent active)
#     = P(∩_p ¬B_p) - P(∩_p (¬A_p ∩ ¬B_p))
#     = prod_p (1 - q_par[p])^{k_p} - prod_p (1 - q_uni[p])^{k_p}
#
# because P(¬A_p ∩ ¬B_p) = 1 - P(A_p ∪ B_p) = 1 - q_uni[p]. No within-prefix
# independence is assumed anywhere. Mutual exclusivity is handled exactly,
# since q_uni differs from q_c + q_par precisely when a subcode maps to both
# branches. The result is an exact expected CCI under the stated model rather
# than an approximation.
#
# The per-group probabilities returned for the QA stage are the *effective*
# (post-suppression) activation probabilities, so that
# sum_g p_eff(g) * w_g reproduces e_cci exactly and the mass-conservation
# table decomposes the reported score additively.
# =============================================================================

#' Group membership of every code in a donor pool.
#'
#' For a resolved prefix this is the cached child-to-group map. For a
#' degenerate prefix the pool is the three-character code itself, whose groups
#' come from matching it directly against the Quan patterns.
#' @keywords internal
.pool_group_map <- function(pool, pref, cache, gn, ap) {
  if (identical(attr(pool, "status", exact = TRUE), "degenerate")) {
    d <- .degenerate_status(pref, gn, ap)
    return(setNames(list(names(d)[d == "certain"]), pref))
  }
  pc <- get_prefix_cache(pref, cache)
  if (is.null(pc)) return(setNames(vector("list", 0L), character(0L)))
  pc$child_groups
}

#' Per-(prefix, age band) hit probabilities for every Charlson group.
#'
#' Returns one row per (pref, age_key, gk) with `q` (group hit probability),
#' `q_par` (any parent of gk) and `q_uni` (gk or any parent). Rows are emitted
#' whenever `q` or `q_par` is positive, which is exactly the set of prefixes
#' that can move the suppression algebra above.
#' @keywords internal
.prefix_group_q <- function(combos, cache, gn, ap, dep_map) {
  rows <- vector("list", nrow(combos))
  for (i in seq_len(nrow(combos))) {
    pref <- combos$pref[i]
    ak   <- combos$age_key[i]
    pool <- prefix_pool(pref, cache, age_idx = if (ak == 0L) NULL else ak)
    if (nrow(pool) == 0L) next
    cg <- .pool_group_map(pool, pref, cache, gn, ap)

    # Logical membership matrix, one column per group.
    memb <- lapply(setNames(gn, gn), function(gk)
      vapply(pool$code_nodot,
             function(cn) gk %in% (cg[[cn]] %||% character(0L)),
             logical(1L)))

    qv <- vapply(gn, function(gk) sum(pool$prob[memb[[gk]]], na.rm = TRUE),
                 numeric(1L))

    q_par <- setNames(numeric(length(gn)), gn)
    q_uni <- setNames(numeric(length(gn)), gn)
    for (gk in gn) {
      pars <- dep_map[[gk]]
      if (length(pars) == 0L) next
      inP <- Reduce(`|`, lapply(pars, function(p) memb[[p]]))
      q_par[gk] <- sum(pool$prob[inP], na.rm = TRUE)
      q_uni[gk] <- sum(pool$prob[memb[[gk]] | inP], na.rm = TRUE)
    }

    keep <- qv > 0 | q_par > 0
    if (!any(keep)) next
    rows[[i]] <- data.table(
      pref    = pref,
      age_key = ak,
      gk      = gn[keep],
      q       = unname(qv[keep]),
      q_par   = unname(q_par[keep]),
      q_uni   = unname(q_uni[keep]),
      is_child = gn[keep] %in% names(dep_map)[lengths(dep_map) > 0L]
    )
  }
  rows <- rows[!vapply(rows, is.null, logical(1L))]
  if (length(rows) == 0L)
    return(data.table(pref = character(0L), age_key = integer(0L),
                      gk = character(0L), q = numeric(0L),
                      q_par = numeric(0L), q_uni = numeric(0L),
                      is_child = logical(0L)))
  rbindlist(rows)
}

#' Single-encounter S2.
#'
#' Defers to `cci_probabilistic_batch()` on length-1 input so the two paths
#' cannot disagree.
#'
#' @param icd_anon character vector of truncated ICD codes.
#' @param quan_map output of `load_quan_map()`.
#' @param cache    output of `precompute_lookups()`.
#' @param age      optional numeric age. When supplied, the age-conditioned
#'   Destatis prior is used for this encounter.
#' @param preserve_multiplicity keep repeated prefixes as separate diagnosis
#'   positions (default TRUE).
#' @return list with elements `e_cci` and `group_prob` (named numeric).
#' @export
cci_probabilistic <- function(icd_anon, quan_map, cache, age = NULL,
                              preserve_multiplicity = TRUE) {
  codes <- truncate_icd(icd_anon, preserve_multiplicity = preserve_multiplicity)
  one <- data.table(diagnosen = paste(codes, collapse = "|"))
  res <- cci_probabilistic_batch(
    one, quan_map, cache, return_group_prob = TRUE,
    age_idx = if (is.null(age)) NULL else age_to_bin_index(age),
    preserve_multiplicity = preserve_multiplicity)
  gp <- setNames(rep(0, length(get(".group_names", envir = cache))),
                 get(".group_names", envir = cache))
  if (nrow(res$group_prob) > 0L) gp[res$group_prob$gk] <- res$group_prob$p
  list(e_cci = unname(res$e_cci[1L]), group_prob = gp)
}

#' Vectorised batch S2.
#'
#' Returns the per-encounter expected CCI as a numeric vector. If
#' `return_group_prob = TRUE`, returns a list with both the e_cci vector
#' and the per-(idx, gk) effective activation probability table used by the
#' QA mass-conservation check.
#'
#' @param dt data.table with column `diagnosen`.
#' @param quan_map output of `load_quan_map()`.
#' @param cache output of `precompute_lookups()`.
#' @param return_group_prob also return the per-(idx, gk) probability table.
#' @param age_idx optional integer vector of length `nrow(dt)` with Destatis
#'   age-band indices (see `age_to_bin_index()`). NULL, or NA for an
#'   individual encounter, selects the marginal prior.
#' @param preserve_multiplicity treat a prefix coded k times as k independent
#'   diagnosis positions (default TRUE).
#' @export
cci_probabilistic_batch <- function(dt, quan_map, cache,
                                    return_group_prob = FALSE,
                                    age_idx = NULL,
                                    preserve_multiplicity = TRUE) {
  gn <- get(".group_names", envir = cache)
  wm <- get(".wt_map",      envir = cache)
  dm <- get(".dep_map",     envir = cache)
  ap <- get(".all_pats",    envir = cache)
  n  <- nrow(dt)

  empty_gp <- data.table(idx = integer(0L), gk = character(0L), p = numeric(0L))

  # Long form (idx, pref) with multiplicity.
  dx_list <- strsplit(as.character(dt$diagnosen), "\\|+")
  long <- data.table(
    idx  = rep.int(seq_len(n), lengths(dx_list)),
    code = toupper(gsub("[^A-Z0-9]", "", unlist(dx_list, use.names = FALSE)))
  )
  long <- long[nchar(code) >= 3L]
  long[, pref := substr(code, 1L, 3L)]
  long <- if (preserve_multiplicity) long[, .(k = .N), by = .(idx, pref)]
          else unique(long[, .(idx, pref)])[, k := 1L][]

  if (nrow(long) == 0L) {
    out <- numeric(n)
    if (return_group_prob) return(list(e_cci = out, group_prob = empty_gp))
    return(out)
  }

  # Age band per encounter; 0 means "use the marginal prior".
  ak_vec <- if (is.null(age_idx)) rep(0L, n) else {
    a <- as.integer(age_idx); a[is.na(a)] <- 0L; a
  }
  long[, age_key := ak_vec[idx]]

  qq <- .prefix_group_q(unique(long[, .(pref, age_key)]),
                        cache, gn, ap, dm)
  if (nrow(qq) == 0L) {
    out <- numeric(n)
    if (return_group_prob) return(list(e_cci = out, group_prob = empty_gp))
    return(out)
  }

  hits <- merge(long, qq, by = c("pref", "age_key"), allow.cartesian = TRUE)
  if (nrow(hits) == 0L) {
    out <- numeric(n)
    if (return_group_prob) return(list(e_cci = out, group_prob = empty_gp))
    return(out)
  }

  # Numerically stable products via log1p; clipped so q = 1 is representable.
  clip <- function(x) pmin(pmax(x, 0), 1 - 1e-12)
  hits[, `:=`(lq   = k * log1p(-clip(q)),
              lpar = k * log1p(-clip(q_par)),
              luni = k * log1p(-clip(q_uni)))]

  enc <- hits[, .(p_raw    = 1 - exp(sum(lq)),
                  keep_par = exp(sum(lpar)),
                  keep_uni = exp(sum(luni)),
                  is_child = any(is_child)),
              by = .(idx, gk)]

  # Effective (post-suppression) activation probability.
  enc[, p := fifelse(is_child, keep_par - keep_uni, p_raw)]
  enc[p < 0, p := 0]

  wt_dt <- data.table(gk = names(wm), w = unname(wm))
  scored <- merge(enc, wt_dt, by = "gk")
  scored[, contribution := p * w]
  ecci <- scored[, .(e_cci = sum(contribution)), by = idx]

  out <- numeric(n)
  out[ecci$idx] <- ecci$e_cci

  if (return_group_prob)
    return(list(e_cci = out, group_prob = enc[p > 0, .(idx, gk, p)]))
  out
}
