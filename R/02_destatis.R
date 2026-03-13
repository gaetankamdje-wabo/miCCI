# =============================================================================
# miCCI — 02_destatis.R
# Destatis 23131-01 + PRECOMPUTATION ENGINE
# All per-prefix distributions computed ONCE, then stored in a hash for O(1) lookup
# =============================================================================

.destatis_env <- new.env(parent = emptyenv())

#' Download and parse Destatis 23131-01 (sex-aggregated)
#' @export
load_destatis <- function(cache = TRUE) {
  if (cache && exists("dt", envir = .destatis_env))
    return(get("dt", envir = .destatis_env))

  url <- "https://www.destatis.de/static/DE/dokumente/5231301237015_SB.xlsx"
  tmp <- tempfile(fileext = ".xlsx")
  download.file(url, tmp, mode = "wb", quiet = TRUE)
  on.exit(unlink(tmp))

  age_labels <- c(
    "unter 1","1 - 5","5 - 10","10 - 15","15-18","18-20","20 - 25",
    "25 - 30","30 - 35","35 - 40","40 - 45","45 - 50","50 - 55",
    "55 - 60","60 - 65","65 - 70","70 - 75","75 - 80","80 - 85",
    "85 - 90","90 - 95","95 u. \u00e4lter")

  raw <- read_excel(tmp, sheet = "23131-01", col_names = FALSE)
  dat <- as.data.frame(raw[-(1:4), ])
  colnames(dat) <- c("code", "sex_raw", "total", age_labels)
  dat <- dat[grepl("^[A-Z]\\d{2}\\.", dat$code), ]

  parse_num <- function(x) {
    x <- as.character(x); x[x %in% c("-",".","/","NA")] <- NA
    x <- gsub("\\.", "", x); x <- gsub(",", ".", x)
    suppressWarnings(as.numeric(x))
  }

  dt <- data.table(code = toupper(trimws(dat$code)), freq_total = parse_num(dat$total))
  for (ag in age_labels) dt[, (ag) := parse_num(dat[[ag]])]

  # Aggregate across sexes
  dt <- dt[, lapply(.SD, sum, na.rm = TRUE), by = code,
           .SDcols = c("freq_total", age_labels)]
  dt[, code_nodot := gsub("\\.", "", code)]
  dt[, code3 := substr(code_nodot, 1, 3)]
  setkey(dt, code_nodot)

  message(sprintf("\u2705 Destatis: %d codes loaded", uniqueN(dt$code_nodot)))
  if (cache) assign("dt", dt, envir = .destatis_env)
  dt
}

# =============================================================================
# PRECOMPUTATION: Build all lookups ONCE
# =============================================================================

#' Precompute ALL per-prefix Charlson group mappings + probability distributions
#'
#' This is the single most important performance function.
#' Called ONCE at pipeline start. Creates:
#'   1) prefix_group_map: for each 3-char prefix, which groups are certain/possible
#'   2) prefix_subcode_map: for each 3-char prefix, subcode probability distribution
#'                          + which Charlson group each subcode maps to
#'
#' @param dt_destatis Destatis data.table
#' @param quan_map Quan mapping
#' @return Environment with lookup tables (hash-based, O(1) access)
precompute_lookups <- function(dt_destatis, quan_map) {
  cache <- new.env(parent = emptyenv(), hash = TRUE)

  # All unique 3-char prefixes in Destatis
  all_prefixes <- unique(dt_destatis$code3)
  n_pref <- length(all_prefixes)
  message(sprintf("  Precomputing %d prefix lookups...", n_pref))

  # Precompute pattern -> group mappings (flat vectors for speed)
  group_names <- names(quan_map)
  all_pats <- list()
  for (gk in group_names) {
    all_pats[[gk]] <- drop_dot(.extract_patterns(quan_map[[gk]]))
  }

  # depends_on + weight lookup
  dep_map <- lapply(setNames(group_names, group_names), function(gk) unlist(quan_map[[gk]]$depends_on))
  wt_map  <- vapply(group_names, function(gk) quan_map[[gk]]$weight, numeric(1))

  # For each prefix: compute group classification + subcode distributions
  for (pref in all_prefixes) {
    children <- dt_destatis[code3 == pref, code_nodot]
    if (length(children) == 0) next

    # ── Group classification (certain/possible/none) ──
    group_status <- setNames(rep("none", length(group_names)), group_names)
    # ── Per-subcode group mapping ──
    child_groups <- setNames(vector("list", length(children)), children)

    for (gk in group_names) {
      pats_gk <- all_pats[[gk]]
      child_matches <- vapply(children, function(ch) {
        any(startsWith(ch, pats_gk) | startsWith(pats_gk, ch))
      }, logical(1))

      n_match <- sum(child_matches)
      if (n_match == length(children)) {
        group_status[gk] <- "certain"
      } else if (n_match > 0) {
        group_status[gk] <- "possible"
      }

      # Tag each matching child
      for (ch in children[child_matches]) {
        child_groups[[ch]] <- c(child_groups[[ch]], gk)
      }
    }

    # ── Marginal subcode distribution (freq_total, no age filter) ──
    child_freqs <- dt_destatis[code_nodot %in% children, .(code_nodot, freq_total)]
    total <- sum(child_freqs$freq_total, na.rm = TRUE)
    child_freqs[, prob := if (total > 0) freq_total / total else 0]

    # Store in cache
    assign(pref, list(
      group_status = group_status,
      child_freqs  = child_freqs,      # code_nodot, freq_total, prob
      child_groups = child_groups,       # code_nodot -> vector of group_keys
      n_children   = length(children)
    ), envir = cache)
  }

  # Store meta
  assign(".group_names", group_names, envir = cache)
  assign(".dep_map", dep_map, envir = cache)
  assign(".wt_map", wt_map, envir = cache)
  assign(".all_pats", all_pats, envir = cache)
  assign(".dt_destatis", dt_destatis, envir = cache)

  # depends_on as data.table for vectorised S1
  dep_rows <- list()
  for (gk in group_names) {
    deps <- dep_map[[gk]]
    if (length(deps) > 0) {
      for (d in deps) dep_rows[[length(dep_rows)+1L]] <- list(child = gk, parent = d)
    }
  }
  dep_dt <- if (length(dep_rows) > 0) rbindlist(dep_rows) else
            data.table(child = character(0), parent = character(0))
  assign(".dep_lookup_dt", dep_dt, envir = cache)

  message(sprintf("  \u2705 %d prefix lookups precomputed", n_pref))
  cache
}

#' Get cached lookup for a 3-char prefix
get_prefix_cache <- function(prefix_3char, cache) {
  pref <- substr(drop_dot(prefix_3char), 1, 3)
  if (exists(pref, envir = cache, inherits = FALSE)) {
    get(pref, envir = cache)
  } else {
    NULL
  }
}

#' Get age-specific subcode probs (on-the-fly, using cached Destatis)
get_subcode_probs_fast <- function(prefix_3char, cache, age = NULL) {
  pc <- get_prefix_cache(prefix_3char, cache)
  if (is.null(pc)) return(data.table(code_nodot = character(0), prob = numeric(0)))

  if (is.null(age) || is.na(age)) return(pc$child_freqs[, .(code_nodot, prob)])

  # Age-specific: recompute from cached Destatis
  age_bin <- age_to_bin(age)
  if (is.na(age_bin)) return(pc$child_freqs[, .(code_nodot, prob)])

  dt_d <- get(".dt_destatis", envir = cache)
  children <- pc$child_freqs$code_nodot
  d <- dt_d[code_nodot %in% children]

  if (!age_bin %in% names(d)) return(pc$child_freqs[, .(code_nodot, prob)])

  agg <- d[, .(freq = sum(get(age_bin), na.rm = TRUE)), by = code_nodot]
  total <- sum(agg$freq, na.rm = TRUE)
  agg[, prob := if (total > 0) freq / total else 0]
  agg[, .(code_nodot, prob)]
}
