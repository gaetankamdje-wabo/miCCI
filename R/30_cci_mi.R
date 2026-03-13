# =============================================================================
# miCCI — 30_cci_mi.R
# S3: MI-CCI — VECTORISED: draw m imputed datasets in bulk, then batch-CCI
# Zero per-encounter loops. ~2 min for 500k @ m=20.
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
  list(mi_cci = round(mean(cci_values), 2))
}

#' VECTORISED S3: bulk imputation + vectorised CCI
cci_mi_batch <- function(dt, quan_map, cache, m = 20L, seed = 42L) {
  set.seed(seed)
  pl <- build_pattern_lookup(quan_map)
  dl <- build_dep_lookup(quan_map)
  n <- nrow(dt)

  # Step 1: For each encounter, extract 3-char prefixes and build donor pools
  dx_list <- strsplit(as.character(dt$diagnosen), "\\|+")
  code_sets <- lapply(dx_list, function(x) {
    unique(substr(toupper(gsub("[^A-Z0-9]", "", x)), 1, 3))
  })

  # Step 2: Build per-prefix donor pools (from cache, ONCE)
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

  message(sprintf("  S3: %d encounters, %d unique prefixes, m=%d", n, length(all_prefs), m))

  # Step 3: For each imputation, draw subcodes and compute CCI in bulk
  cci_sum <- numeric(n)
  for (imp in seq_len(m)) {
    # Draw one subcode per prefix per encounter
    imputed_dx <- vapply(seq_len(n), function(i) {
      prefs <- code_sets[[i]]
      drawn <- vapply(prefs, function(pref) {
        pool <- donor_pools[[pref]]
        add_dot4(pool$code_nodot[sample.int(nrow(pool), 1L, prob = pool$prob)])
      }, character(1))
      paste(drawn, collapse = "|")
    }, character(1))

    # Compute CCI for all encounters at once using vectorised gold batch
    imp_dt <- data.table(diagnosen = imputed_dx)
    cci_sum <- cci_sum + cci_gold_batch(imp_dt, quan_map, pl, dl)
  }

  round(cci_sum / m, 2)
}
