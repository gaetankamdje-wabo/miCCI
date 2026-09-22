# =============================================================================
# miCCI / 30_cci_mi.R
# S3 - Multiple-Imputation CCI.
#
# For each truncated prefix, draw a full-length subcode from its empirical
# subcode distribution m times, recompute the gold CCI on each imputed
# code set, and return the mean across the m imputations. Optionally,
# track per-group activation rates across draws for the QA stage.
# =============================================================================

#' Single-encounter S3.
#'
#' For length-1 input, defers to `cci_mi_batch()` so behaviour is identical
#' for one or many encounters.
#'
#' @param icd_anon character vector of truncated ICD codes for one encounter.
#' @param quan_map output of `load_quan_map()`.
#' @param cache    output of `precompute_lookups()`.
#' @param age      optional numeric age (currently unused; reserved for
#'   age-stratified imputation in later versions).
#' @param m        number of imputations.
#' @param seed     RNG seed (local; never alters the caller's RNG state).
#' @return list with element `mi_cci`.
#' @export
cci_mi <- function(icd_anon, quan_map, cache,
                   age = NULL, m = 20L, seed = 42L,
                   preserve_multiplicity = TRUE) {
  # Defer to batch on length-1 input so the two paths cannot disagree.
  codes <- truncate_icd(icd_anon, preserve_multiplicity = preserve_multiplicity)
  one <- data.table(diagnosen = paste(codes, collapse = "|"))
  v <- cci_mi_batch(one, quan_map, cache, m = m, seed = seed,
                    age_idx = if (is.null(age)) NULL else age_to_bin_index(age),
                    preserve_multiplicity = preserve_multiplicity)
  list(mi_cci = unname(v[1L]))
}

#' Vectorised S3 batch.
#'
#' @param dt   data.table with column `diagnosen`.
#' @param quan_map output of `load_quan_map()`.
#' @param cache output of `precompute_lookups()`.
#' @param m number of imputations.
#' @param seed RNG seed (local).
#' @param return_group_count If TRUE, also return a per-(idx, gk) trigger
#'   count across the m draws, used by the QA mass-conservation check.
#' @return numeric vector of length nrow(dt) (or a list if
#'   `return_group_count = TRUE`).
#' @export
cci_mi_batch <- function(dt, quan_map, cache,
                         m = 20L, seed = 42L,
                         return_group_count = FALSE,
                         age_idx = NULL,
                         preserve_multiplicity = TRUE,
                         return_rounds = FALSE) {
  pl <- build_pattern_lookup(quan_map)
  dl <- build_dep_lookup(quan_map)
  n  <- nrow(dt)

  # One prefix slot per coded diagnosis position when multiplicity is kept, so
  # an encounter carrying E11 twice draws two subcodes and can activate the
  # group from either. Collapsing to unique() caps the draw at one subcode per
  # prefix and undercounts comorbidities from the same ICD-10 block.
  dx_list <- strsplit(as.character(dt$diagnosen), "\\|+")
  code_sets <- lapply(dx_list, function(x) {
    p <- substr(toupper(gsub("[^A-Z0-9]", "", x)), 1L, 3L)
    p <- p[nchar(p) >= 3L]
    if (preserve_multiplicity) p else unique(p)
  })

  # Donor pools are keyed by prefix AND age band, so an encounter with a known
  # age draws from P(subcode | prefix, age band) instead of the marginal.
  ak_vec <- if (is.null(age_idx)) rep(0L, n) else {
    a <- as.integer(age_idx); a[is.na(a)] <- 0L; a
  }
  key_sets <- lapply(seq_len(n), function(i) {
    if (length(code_sets[[i]]) == 0L) character(0L)
    else paste0(code_sets[[i]], "@", ak_vec[i])
  })

  all_keys <- unique(unlist(key_sets, use.names = FALSE))
  donor_pools <- list()
  n_degenerate <- 0L; n_agefall <- 0L
  for (key in all_keys) {
    parts <- strsplit(key, "@", fixed = TRUE)[[1L]]
    pref  <- parts[1L]; ak <- as.integer(parts[2L])
    pool  <- prefix_pool(pref, cache, age_idx = if (ak == 0L) NULL else ak)
    if (nrow(pool) == 0L) pool <- data.table(code_nodot = pref, prob = 1)
    if (identical(attr(pool, "status", exact = TRUE), "degenerate"))
      n_degenerate <- n_degenerate + 1L
    if (isTRUE(attr(pool, "age_fallback", exact = TRUE)))
      n_agefall <- n_agefall + 1L
    donor_pools[[key]] <- pool
  }

  message(sprintf(paste0("S3: n=%d encounters, %d prefix-by-age pools, m=%d",
                         " (%d degenerate, %d age fallbacks)"),
                  n, length(all_keys), m, n_degenerate, n_agefall))

  cci_sum <- numeric(n)
  # Per-round scores, kept only on request. Round r consumes the same random
  # numbers whatever m is, so the first k columns of an m-round run ARE the
  # k-round run. The sensitivity analysis exploits that to evaluate a whole
  # grid of m from a single run of the largest one.
  rounds <- if (return_rounds) matrix(0, nrow = n, ncol = m) else NULL
  group_count <- if (return_group_count)
    integer(0L) else NULL  # accumulated below
  if (return_group_count) {
    group_count_dt <- data.table(idx = integer(0L), gk = character(0L), n_active = integer(0L))
  }

  .with_local_seed(seed, {
    for (imp in seq_len(m)) {
      imputed_dx <- vapply(seq_len(n), function(i) {
        keys <- key_sets[[i]]
        if (length(keys) == 0L) return("")
        drawn <- vapply(keys, function(key) {
          pool <- donor_pools[[key]]
          add_dot4(pool$code_nodot[sample.int(nrow(pool), 1L, prob = pool$prob)])
        }, character(1L))
        paste(drawn, collapse = "|")
      }, character(1L))

      imp_dt <- data.table(diagnosen = imputed_dx)

      if (return_group_count) {
        ga <- .gold_active_long_internal(imp_dt, quan_map, pl, dl)
        if (nrow(ga) > 0L) {
          ga[, n_active := 1L]
          group_count_dt <- rbindlist(list(group_count_dt, ga), use.names = TRUE, fill = TRUE)
        }
      }
      sc <- cci_gold_batch(imp_dt, quan_map, pl, dl)
      cci_sum <- cci_sum + sc
      if (return_rounds) rounds[, imp] <- sc
    }
  })

  if (return_rounds) return(rounds)
  out <- cci_sum / m
  if (return_group_count) {
    gp <- group_count_dt[, .(p = sum(n_active) / m), by = .(idx, gk)]
    return(list(mi_cci = out, group_prob = gp))
  }
  out
}

# Internal helper: same logic as cci_gold_batch but returns the (idx, gk)
# long form after dependency suppression. Kept here to avoid a circular
# dependency with the QA module.
.gold_active_long_internal <- function(dt, quan_map,
                                       pattern_lookup = NULL,
                                       dep_lookup     = NULL) {
  if (is.null(pattern_lookup)) pattern_lookup <- build_pattern_lookup(quan_map)
  if (is.null(dep_lookup))     dep_lookup     <- build_dep_lookup(quan_map)
  pl <- pattern_lookup; dl <- dep_lookup

  dx_list <- strsplit(as.character(dt$diagnosen), "\\|+")
  long <- data.table(
    idx      = rep.int(seq_len(nrow(dt)), lengths(dx_list)),
    code_raw = unlist(dx_list, use.names = FALSE)
  )
  long[, cn := toupper(gsub("[^A-Z0-9]", "", code_raw))]
  long <- long[nchar(cn) >= 3L]
  if (nrow(long) == 0L) return(data.table(idx = integer(0L), gk = character(0L)))

  unique_codes <- unique(long$cn)
  pat_vec <- pl$pat
  gk_vec  <- pl$gk
  rows <- list()
  for (uc in unique_codes) {
    hits <- gk_vec[startsWith(uc, pat_vec) | startsWith(pat_vec, uc)]
    if (length(hits) > 0L)
      rows[[length(rows) + 1L]] <- data.table(cn = uc, gk = unique(hits))
  }
  if (length(rows) == 0L)
    return(data.table(idx = integer(0L), gk = character(0L)))

  cgm <- rbindlist(rows); setkey(cgm, cn); setkey(long, cn)
  hits <- cgm[long, on = "cn", nomatch = NULL, allow.cartesian = TRUE]
  triggered <- unique(hits[, .(idx, gk)])
  if (nrow(dl) > 0L && nrow(triggered) > 0L) {
    suppress <- merge(triggered, dl, by.x = "gk", by.y = "child",
                      allow.cartesian = TRUE)
    suppress <- merge(suppress, triggered,
                      by.x = c("idx", "parent"), by.y = c("idx", "gk"))
    if (nrow(suppress) > 0L) {
      suppress_keys <- unique(suppress[, .(idx, gk)])
      triggered <- triggered[!suppress_keys, on = c("idx", "gk")]
    }
  }
  triggered
}
