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
                   age = NULL, m = 20L, seed = 42L) {
  # Defer to batch on length-1 input so the two paths cannot disagree.
  one <- data.table(diagnosen = paste(unique(icd3(icd_anon)), collapse = "|"))
  v <- cci_mi_batch(one, quan_map, cache, m = m, seed = seed)
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
                         return_group_count = FALSE) {
  pl <- build_pattern_lookup(quan_map)
  dl <- build_dep_lookup(quan_map)
  n  <- nrow(dt)

  dx_list <- strsplit(as.character(dt$diagnosen), "\\|+")
  code_sets <- lapply(dx_list, function(x) {
    unique(substr(toupper(gsub("[^A-Z0-9]", "", x)), 1L, 3L))
  })
  code_sets <- lapply(code_sets, function(x) x[nchar(x) >= 3L])

  all_prefs <- unique(unlist(code_sets, use.names = FALSE))
  donor_pools <- list()
  for (pref in all_prefs) {
    pc <- get_prefix_cache(pref, cache)
    if (is.null(pc)) {
      donor_pools[[pref]] <- data.table(code_nodot = pref, prob = 1)
    } else {
      p <- pc$child_freqs[prob > 0]
      donor_pools[[pref]] <- if (nrow(p) > 0L) p else data.table(code_nodot = pref, prob = 1)
    }
  }

  message(sprintf("S3: n=%d encounters, %d unique prefixes, m=%d",
                  n, length(all_prefs), m))

  cci_sum <- numeric(n)
  group_count <- if (return_group_count)
    integer(0L) else NULL  # accumulated below
  if (return_group_count) {
    group_count_dt <- data.table(idx = integer(0L), gk = character(0L), n_active = integer(0L))
  }

  .with_local_seed(seed, {
    for (imp in seq_len(m)) {
      imputed_dx <- vapply(seq_len(n), function(i) {
        prefs <- code_sets[[i]]
        if (length(prefs) == 0L) return("")
        drawn <- vapply(prefs, function(pref) {
          pool <- donor_pools[[pref]]
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
      cci_sum <- cci_sum + cci_gold_batch(imp_dt, quan_map, pl, dl)
    }
  })

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
