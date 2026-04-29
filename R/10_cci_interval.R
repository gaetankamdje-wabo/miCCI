# =============================================================================
# miCCI / 10_cci_interval.R
# S1 - Interval CCI.
#
# A truncated three-character ICD prefix may map to a Charlson group with
# certainty (every full-length child of the prefix is in the group) or only
# possibly (some children are, some aren't). S1 returns:
#
#   cci_min : sum of weights for groups triggered with certainty.
#   cci_max : sum of weights for groups triggered with certainty OR possibly.
#   cci_mid : (cci_min + cci_max) / 2, kept as a true float.
#   width   : cci_max - cci_min.
#
# By construction, cci_min <= cci_gold <= cci_max for any encounter whose
# original full-length codes appear in the frequency table that built the
# cache. The test suite asserts this envelope explicitly.
# =============================================================================

#' Single-encounter S1 interval CCI.
#'
#' For a cohort, prefer `cci_interval_batch()`. The two share the same
#' lookup engine and produce identical results.
#'
#' @param icd_anon Character vector of (already truncated) ICD codes.
#' @param quan_map Output of `load_quan_map()`.
#' @param cache    Output of `precompute_lookups()`.
#' @return list with `cci_min`, `cci_max`, `cci_mid`, `interval_width`.
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
    if (length(deps) > 0L && any(deps %in% certain))
      certain <- setdiff(certain, gk)
  }
  any_hit <- names(status)[status %in% c("possible", "certain")]
  cci_min <- sum(wm[certain])
  cci_max <- sum(wm[any_hit])
  list(
    cci_min        = as.numeric(cci_min),
    cci_max        = as.numeric(cci_max),
    cci_mid        = (cci_min + cci_max) / 2,
    interval_width = as.numeric(cci_max - cci_min)
  )
}

#' Vectorised S1 batch.
#'
#' Returns a data.table with one row per encounter and columns
#' `cci_min`, `cci_max`, `cci_mid`, `interval_width`.
#' @export
cci_interval_batch <- function(dt, quan_map, cache) {
  gn <- get(".group_names", envir = cache)
  wm <- get(".wt_map",      envir = cache)
  dl <- get(".dep_lookup_dt", envir = cache)

  # Build a (prefix, group, status) table from the cache.
  all_pref <- ls(cache)
  all_pref <- all_pref[!startsWith(all_pref, ".")]
  pg_rows <- list()
  for (pref in all_pref) {
    pc <- get(pref, envir = cache, inherits = FALSE)
    for (gk in gn) {
      st <- pc$group_status[gk]
      if (st != "none") {
        pg_rows[[length(pg_rows) + 1L]] <-
          data.table(pref = pref, gk = gk, st = unname(st))
      }
    }
  }
  pg <- if (length(pg_rows) > 0L) rbindlist(pg_rows)
        else data.table(pref = character(0L), gk = character(0L), st = character(0L))

  # Long form (idx, prefix) for the cohort.
  dx_list <- strsplit(as.character(dt$diagnosen), "\\|+")
  long <- data.table(
    idx  = rep.int(seq_len(nrow(dt)), lengths(dx_list)),
    code = toupper(gsub("[^A-Z0-9]", "", unlist(dx_list, use.names = FALSE)))
  )
  long <- long[nchar(code) >= 3L]
  long[, pref := substr(code, 1L, 3L)]
  long <- unique(long[, .(idx, pref)])

  # Resolve (idx, group) status: certain wins over possible.
  hits <- merge(long, pg, by = "pref", allow.cartesian = TRUE)
  enc_grp <- hits[, .(st = if (any(st == "certain")) "certain" else "possible"),
                  by = .(idx, gk)]

  # cci_max - all hit groups (no dependency suppression for the upper bound).
  wt_dt  <- data.table(gk = names(wm), w = unname(wm)); setkey(wt_dt, gk)
  max_dt <- merge(enc_grp, wt_dt, by = "gk")
  cci_max_dt <- max_dt[, .(cci_max = sum(w)), by = idx]

  # cci_min - certain only, with dependency suppression via anti-join.
  certain_dt <- enc_grp[st == "certain", .(idx, gk)]
  if (nrow(certain_dt) > 0L && nrow(dl) > 0L) {
    suppress <- merge(certain_dt, dl,
                      by.x = "gk", by.y = "child", allow.cartesian = TRUE)
    suppress <- merge(suppress, certain_dt,
                      by.x = c("idx", "parent"), by.y = c("idx", "gk"))
    if (nrow(suppress) > 0L) {
      suppress_keys <- unique(suppress[, .(idx, gk)])
      certain_dt <- certain_dt[!suppress_keys, on = c("idx", "gk")]
    }
  }
  min_dt <- merge(certain_dt, wt_dt, by = "gk")
  cci_min_dt <- min_dt[, .(cci_min = sum(w)), by = idx]

  out <- data.table(idx = seq_len(nrow(dt)))
  out <- merge(out, cci_min_dt, by = "idx", all.x = TRUE)
  out <- merge(out, cci_max_dt, by = "idx", all.x = TRUE)
  out[is.na(cci_min), cci_min := 0]
  out[is.na(cci_max), cci_max := 0]
  out[, cci_min        := as.numeric(cci_min)]
  out[, cci_max        := as.numeric(cci_max)]
  out[, cci_mid        := (cci_min + cci_max) / 2]
  out[, interval_width := cci_max - cci_min]
  setorder(out, idx)
  out[, .(cci_min, cci_max, cci_mid, interval_width)]
}
