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
# cache.
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
  ap <- get(".all_pats",    envir = cache)

  status <- setNames(rep("none", length(gn)), gn)
  for (code in codes) {
    pc <- get_prefix_cache(code, cache)
    # A prefix with no children in the reference table stands for itself
    # (see prefix_pool() in 02_destatis.R), so its groups are certain, not
    # absent. Skipping it here is what used to zero out AIDS under S1.
    gs <- if (is.null(pc)) .degenerate_status(code, gn, ap) else pc$group_status
    for (gk in gn) {
      s <- gs[gk]
      if (s == "certain") status[gk] <- "certain"
      else if (s == "possible" && status[gk] == "none") status[gk] <- "possible"
    }
  }

  # Hierarchical exclusion is a property of the Charlson score, not of the
  # bound being computed, so it applies to both ends of the interval. In every
  # Quan dependency pair the parent outweighs the child (diabetes 2 over 1,
  # liver 3 over 1, malignancy 6 over 2), so the highest-scoring feasible
  # assignment is always "parent active, child suppressed". Suppressing the
  # child therefore yields the true maximum rather than an inflated sum of
  # mutually exclusive weights.
  .suppress <- function(active) {
    for (gk in active) {
      deps <- dm[[gk]]
      if (length(deps) > 0L && any(deps %in% active))
        active <- setdiff(active, gk)
    }
    active
  }
  certain <- .suppress(names(status)[status == "certain"])
  any_hit <- .suppress(names(status)[status %in% c("possible", "certain")])
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
#' If `return_group_prob = TRUE`, also returns a `group_prob` element:
#' a data.table with columns `idx`, `gk`, `p`, where `p = 1` for certain
#' groups and `p = 0.5` for possible groups (the mid-point assumption).
#' @export
cci_interval_batch <- function(dt, quan_map, cache,
                               return_group_prob = FALSE) {
  gn <- get(".group_names", envir = cache)
  wm <- get(".wt_map",      envir = cache)
  dl <- get(".dep_lookup_dt", envir = cache)
  ap <- get(".all_pats",    envir = cache)

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

  # Long form (idx, prefix) for the cohort. S1 reads only the group status of a
  # prefix, which does not change when the same prefix is coded twice, so the
  # deduplication here is a performance step and not a multiplicity decision.
  dx_list <- strsplit(as.character(dt$diagnosen), "\\|+")
  long <- data.table(
    idx  = rep.int(seq_len(nrow(dt)), lengths(dx_list)),
    code = toupper(gsub("[^A-Z0-9]", "", unlist(dx_list, use.names = FALSE)))
  )
  long <- long[nchar(code) >= 3L]
  long[, pref := substr(code, 1L, 3L)]
  long <- unique(long[, .(idx, pref)])

  # Degenerate prefixes: present in the cohort, absent from the reference
  # table. They stand for themselves, so their groups are certain.
  deg <- setdiff(unique(long$pref), all_pref)
  if (length(deg) > 0L) {
    deg_rows <- list()
    for (pref in deg) {
      st <- .degenerate_status(pref, gn, ap)
      hit <- names(st)[st == "certain"]
      if (length(hit) > 0L)
        deg_rows[[length(deg_rows) + 1L]] <-
          data.table(pref = pref, gk = hit, st = "certain")
    }
    if (length(deg_rows) > 0L)
      pg <- rbindlist(c(list(pg), deg_rows), use.names = TRUE)
  }

  # Resolve (idx, group) status: certain wins over possible.
  hits <- merge(long, pg, by = "pref", allow.cartesian = TRUE)
  enc_grp <- hits[, .(st = if (any(st == "certain")) "certain" else "possible"),
                  by = .(idx, gk)]

  wt_dt  <- data.table(gk = names(wm), w = unname(wm)); setkey(wt_dt, gk)

  # Hierarchical exclusion applies to BOTH bounds. Previously the upper
  # bound summed every certain-or-possible group regardless of the Quan
  # exclusion rules, so an encounter whose prefix could reach both mild and
  # severe liver disease was credited 1 + 3 even though no single subcode
  # assignment scores that way. Because the parent outweighs the child in all
  # three Quan pairs, the suppressed sum is the true reachable maximum.
  .suppress_dt <- function(x) {
    if (nrow(x) == 0L || nrow(dl) == 0L) return(x)
    sup <- merge(x, dl, by.x = "gk", by.y = "child", allow.cartesian = TRUE)
    sup <- merge(sup, x, by.x = c("idx", "parent"), by.y = c("idx", "gk"))
    if (nrow(sup) == 0L) return(x)
    x[!unique(sup[, .(idx, gk)]), on = c("idx", "gk")]
  }

  # cci_max - certain or possible, after exclusion suppression.
  any_dt <- .suppress_dt(enc_grp[, .(idx, gk)])
  max_dt <- merge(any_dt, wt_dt, by = "gk")
  cci_max_dt <- max_dt[, .(cci_max = sum(w)), by = idx]

  # cci_min - certain only, after exclusion suppression.
  certain_dt <- .suppress_dt(enc_grp[st == "certain", .(idx, gk)])
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
  if (!return_group_prob)
    return(out[, .(cci_min, cci_max, cci_mid, interval_width)])

  # Build (idx, gk, p) from the two suppressed bound sets, so that
  # sum_g p(g) * w_g reproduces cci_mid exactly. A group that survives
  # suppression in both bounds contributes 1, one that survives only in the
  # upper bound contributes the mid-point 0.5, and a group suppressed out of
  # both contributes nothing. Reading the raw certain/possible status instead
  # would credit mass to groups that the exclusion rules have already removed
  # from the score.
  gp <- merge(any_dt[, .(idx, gk, in_any = 1L)],
              certain_dt[, .(idx, gk, in_cert = 1L)],
              by = c("idx", "gk"), all.x = TRUE)
  gp[is.na(in_cert), in_cert := 0L]
  gp[, p := (in_cert + in_any) / 2]
  gp <- gp[, .(idx, gk, p)]
  list(
    cci_min    = out$cci_min,
    cci_max    = out$cci_max,
    cci_mid    = out$cci_mid,
    interval_width = out$interval_width,
    group_prob = gp
  )
}