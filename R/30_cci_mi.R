# =============================================================================
# miCCI v1.0.0 — 30_cci_mi.R
# S3: MI-CCI — bulk imputation, vectorised gold per draw, mean across draws
# v1.0 change: returns raw float vector (no round())
# =============================================================================

#' @export
cci_mi <- function(icd_anon, quan_map, cache, age = NULL, m = 20L, seed = 42L) {
  set.seed(seed)
  codes <- unique(icd3(icd_anon))
  if (length(codes) == 0) return(list(mi_cci = 0))
  pools <- lapply(codes, function(code) {
    pc <- get_prefix_cache(code, cache)
    if (is.null(pc)) return(data.table(code_nodot = drop_dot(code), prob = 1))
    p <- pc$child_freqs[prob > 0]
    if (nrow(p) == 0) return(data.table(code_nodot = drop_dot(code), prob = 1))
    p
  })
  cci_values <- integer(m)
  for (i in seq_len(m)) {
    imputed <- vapply(pools, function(pool)
      add_dot4(pool$code_nodot[sample.int(nrow(pool), 1L, prob = pool$prob)]),
      character(1))
    cci_values[i] <- cci_gold(imputed, quan_map)$cci
  }
  list(mi_cci = mean(cci_values))
}

#' VECTORISED S3 — bulk imputation; returns float vector
cci_mi_batch <- function(dt, quan_map, cache, m = 20L, seed = 42L) {
  set.seed(seed)
  pl <- build_pattern_lookup(quan_map)
  dl <- build_dep_lookup(quan_map)
  n <- nrow(dt)

  dx_list <- strsplit(as.character(dt$diagnosen), "\\|+")
  code_sets <- lapply(dx_list, function(x) {
    unique(substr(toupper(gsub("[^A-Z0-9]", "", x)), 1, 3))
  })

  all_prefs <- unique(unlist(code_sets))
  donor_pools <- list()
  for (pref in all_prefs) {
    pc <- get_prefix_cache(pref, cache)
    if (is.null(pc)) {
      donor_pools[[pref]] <- data.table(code_nodot = pref, prob = 1)
    } else {
      p <- pc$child_freqs[prob > 0]
      donor_pools[[pref]] <- if (nrow(p) > 0) p else data.table(code_nodot = pref, prob = 1)
    }
  }

  message(sprintf("  S3: %d encounters, %d unique prefixes, m=%d",
                  n, length(all_prefs), m))

  cci_sum <- numeric(n)
  for (imp in seq_len(m)) {
    imputed_dx <- vapply(seq_len(n), function(i) {
      prefs <- code_sets[[i]]
      if (length(prefs) == 0) return("")
      drawn <- vapply(prefs, function(pref) {
        pool <- donor_pools[[pref]]
        add_dot4(pool$code_nodot[sample.int(nrow(pool), 1L, prob = pool$prob)])
      }, character(1))
      paste(drawn, collapse = "|")
    }, character(1))

    imp_dt <- data.table(diagnosen = imputed_dx)
    cci_sum <- cci_sum + cci_gold_batch(imp_dt, quan_map, pl, dl)
  }

  cci_sum / m   # raw float, no round
}
