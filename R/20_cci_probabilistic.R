# =============================================================================
# miCCI — 20_cci_probabilistic.R
# S2: E-CCI — FULLY VECTORISED via data.table joins
# =============================================================================

#' @export
cci_probabilistic <- function(icd_anon, quan_map, cache, age = NULL) {
  codes <- unique(icd3(icd_anon))
  gn <- get(".group_names", envir = cache)
  dm <- get(".dep_map", envir = cache)
  wm <- get(".wt_map", envir = cache)
  gp <- setNames(numeric(length(gn)), gn)
  for (code in codes) {
    pc <- get_prefix_cache(code, cache)
    if (is.null(pc)) next
    probs <- pc$child_freqs
    for (gk in gn) {
      matching <- vapply(probs$code_nodot, function(cn)
        gk %in% (pc$child_groups[[cn]] %||% character(0)), logical(1))
      gp[gk] <- max(gp[gk], sum(probs$prob[matching], na.rm = TRUE))
    }
  }
  cont <- numeric(length(gn)); names(cont) <- gn
  for (gk in gn) {
    if (gp[gk] <= 0) next
    deps <- dm[[gk]]
    if (length(deps) > 0 && any(gp[deps] > 0, na.rm = TRUE)) next
    cont[gk] <- gp[gk] * wm[gk]
  }
  list(e_cci = round(sum(cont), 2))
}

#' FULLY VECTORISED S2 batch
cci_probabilistic_batch <- function(dt, quan_map, cache) {
  gn <- get(".group_names", envir = cache)
  wm <- get(".wt_map", envir = cache)
  dl <- get(".dep_lookup_dt", envir = cache)

  # Build master (prefix, group, prob) from cache — ONCE
  all_pref <- ls(cache)
  all_pref <- all_pref[!startsWith(all_pref, ".")]
  pg_rows <- list()
  for (pref in all_pref) {
    pc <- get(pref, envir = cache)
    cdt <- pc$child_freqs
    cg  <- pc$child_groups
    for (gk in gn) {
      p <- sum(cdt$prob[vapply(cdt$code_nodot, function(cn)
        gk %in% (cg[[cn]] %||% character(0)), logical(1))], na.rm = TRUE)
      if (p > 0) pg_rows[[length(pg_rows) + 1L]] <- list(pref = pref, gk = gk, gp = p)
    }
  }
  pgp <- rbindlist(pg_rows)
  message(sprintf("  S2: prefix-group prob table: %d entries", nrow(pgp)))

  # Explode encounters -> prefixes
  dx_list <- strsplit(as.character(dt$diagnosen), "\\|+")
  long <- data.table(
    idx  = rep(seq_len(nrow(dt)), lengths(dx_list)),
    code = toupper(gsub("[^A-Z0-9]", "", unlist(dx_list, use.names = FALSE)))
  )
  long <- long[nchar(code) >= 3]
  long[, pref := substr(code, 1, 3)]
  long <- unique(long[, .(idx, pref)])

  # Join -> (idx, gk, gp)
  hits <- merge(long, pgp, by = "pref", allow.cartesian = TRUE)
  enc_gp <- hits[, .(gp = max(gp)), by = .(idx, gk)]

  # depends_on: remove child where parent also present
  if (nrow(dl) > 0 && nrow(enc_gp) > 0) {
    # Find (idx, child_gk) where parent_gk also in enc_gp
    check <- merge(enc_gp[, .(idx, child_gk = gk)], dl, by.x = "child_gk", by.y = "child",
                   allow.cartesian = TRUE)
    # Keep only where parent actually exists for same idx
    check <- check[enc_gp, on = .(idx = idx, parent = gk), nomatch = 0L]
    if (nrow(check) > 0) {
      suppress_keys <- unique(check[, .(idx, gk = child_gk)])
      enc_gp <- enc_gp[!suppress_keys, on = c("idx", "gk")]
    }
  }

  # Weight * prob, sum per encounter
  wt_dt <- data.table(gk = names(wm), w = unname(wm))
  scored <- merge(enc_gp, wt_dt, by = "gk")
  scored[, contribution := gp * w]
  ecci <- scored[, .(e_cci = round(sum(contribution), 2)), by = idx]

  out <- rep(0, nrow(dt))
  out[ecci$idx] <- ecci$e_cci
  out
}
