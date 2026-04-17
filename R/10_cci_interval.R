# =============================================================================
# miCCI v1.0.0 — 10_cci_interval.R
# S1: Interval CCI — VECTORISED via data.table joins
# v1.0 change: cci_mid is returned as a true float (no integer coercion)
# =============================================================================

#' @export
cci_interval <- function(icd_anon, quan_map, cache) {
  codes <- unique(icd3(icd_anon))
  gn <- get(".group_names", envir = cache)
  dm <- get(".dep_map",     envir = cache)
  wm <- get(".wt_map",      envir = cache)
  status <- setNames(rep("none", length(gn)), gn)
  for (code in codes) {
    pc <- get_prefix_cache(code, cache)
    if (is.null(pc)) next
    for (gk in gn) {
      s <- pc$group_status[gk]
      if (s == "certain") status[gk] <- "certain"
      else if (s == "possible" && status[gk] == "none") status[gk] <- "possible"
    }
  }
  certain <- names(status)[status == "certain"]
  for (gk in certain) {
    deps <- dm[[gk]]
    if (length(deps) > 0 && any(deps %in% certain)) certain <- setdiff(certain, gk)
  }
  any_hit <- names(status)[status %in% c("possible", "certain")]
  cci_min <- sum(wm[certain])
  cci_max <- sum(wm[any_hit])
  list(cci_min        = as.numeric(cci_min),
       cci_max        = as.numeric(cci_max),
       cci_mid        = (cci_min + cci_max) / 2,
       interval_width = as.numeric(cci_max - cci_min))
}

#' VECTORISED S1 batch via memoization by unique code-set signature
cci_interval_batch <- function(dt, quan_map, cache) {
  gn <- get(".group_names", envir = cache)
  dm <- get(".dep_map",     envir = cache)
  wm <- get(".wt_map",      envir = cache)

  all_pref <- ls(cache)
  all_pref <- all_pref[!startsWith(all_pref, ".")]
  pg_rows <- list()
  for (pref in all_pref) {
    pc <- get(pref, envir = cache)
    for (gk in gn) {
      if (pc$group_status[gk] != "none") {
        pg_rows[[length(pg_rows) + 1L]] <- list(pref = pref, gk = gk, st = pc$group_status[gk])
      }
    }
  }
  pg <- rbindlist(pg_rows)

  dx_list <- strsplit(as.character(dt$diagnosen), "\\|+")
  lens <- lengths(dx_list)
  long <- data.table(
    idx  = rep(seq_len(nrow(dt)), lens),
    code = toupper(gsub("[^A-Z0-9]", "", unlist(dx_list, use.names = FALSE)))
  )
  long <- long[nchar(code) >= 3]
  long[, pref := substr(code, 1, 3)]
  long <- unique(long[, .(idx, pref)])

  hits <- merge(long, pg, by = "pref", allow.cartesian = TRUE)
  enc_grp <- hits[, .(st = if (any(st == "certain")) "certain" else "possible"),
                  by = .(idx, gk)]

  # CCI_max — all triggered groups
  max_dt <- merge(enc_grp, data.table(gk = names(wm), w = unname(wm)), by = "gk")
  cci_max_dt <- max_dt[, .(cci_max = sum(w)), by = idx]

  # CCI_min — only certain groups, with depends_on
  certain_dt <- enc_grp[st == "certain"]
  if (nrow(certain_dt) > 0 && nrow(get(".dep_lookup_dt", envir = cache)) > 0) {
    dl <- get(".dep_lookup_dt", envir = cache)
    suppress <- merge(certain_dt[, .(idx, gk)], dl, by.x = "gk", by.y = "child",
                      allow.cartesian = TRUE)
    suppress <- merge(suppress, certain_dt[, .(idx, gk)],
                      by.x = c("idx", "parent"), by.y = c("idx", "gk"))
    if (nrow(suppress) > 0) {
      certain_dt <- fsetdiff(certain_dt[, .(idx, gk)],
                             unique(suppress[, .(idx, gk)]))
    } else {
      certain_dt <- certain_dt[, .(idx, gk)]
    }
  } else {
    certain_dt <- certain_dt[, .(idx, gk)]
  }
  min_wt <- merge(certain_dt, data.table(gk = names(wm), w = unname(wm)), by = "gk")
  cci_min_dt <- min_wt[, .(cci_min = sum(w)), by = idx]

  out <- data.table(idx = seq_len(nrow(dt)))
  out <- merge(out, cci_min_dt, by = "idx", all.x = TRUE)
  out <- merge(out, cci_max_dt, by = "idx", all.x = TRUE)
  out[is.na(cci_min), cci_min := 0]
  out[is.na(cci_max), cci_max := 0]
  out[, cci_min := as.numeric(cci_min)]
  out[, cci_max := as.numeric(cci_max)]
  out[, cci_mid := (cci_min + cci_max) / 2]
  out[, interval_width := cci_max - cci_min]
  setorder(out, idx)
  out[, .(cci_min, cci_max, cci_mid, interval_width)]
}
