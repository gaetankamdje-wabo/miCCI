# =============================================================================
# miCCI / 02_destatis.R
# Destatis 23131-01 loader + per-prefix lookup precomputation.
#
# The user can either:
#   (a) call load_destatis() to download the German Federal Statistical
#       Office's age-stratified diagnosis-frequency table, or
#   (b) supply a custom data.table of the same schema via the `freq_table`
#       argument of precompute_lookups(), which is how a site uses its own
#       population frequencies instead of Destatis.
# =============================================================================

# Cache for the parsed Destatis table. Keyed by URL so a site that swaps
# the URL gets a clean reload without restarting the session.
.destatis_env <- new.env(parent = emptyenv())

#' Download and parse Destatis 23131-01 (sex-aggregated, age-stratified).
#'
#' The table contains, per ICD-10-GM 4-character code, the national
#' inpatient frequency by age band. miCCI uses it as the prior for S2 to S5.
#'
#' Network access is required on first call. Subsequent calls return the
#' cached table unless `force = TRUE`.
#'
#' @param url URL of the Destatis xlsx workbook. Defaults to the official
#'   23131-01 publication URL.
#' @param cache If TRUE, store the parsed table in a package-level
#'   environment so further calls in the same session are free.
#' @param force If TRUE, ignore the cache and reload (and reparse) from
#'   the URL. Useful when Destatis releases an updated annual table.
#' @param quiet Suppress download progress output.
#'
#' @return A data.table with columns `code`, `code_nodot`, `code3`,
#'   `freq_total`, and one column per Destatis age band.
#'
#' @export
load_destatis <- function(url   = "https://www.destatis.de/static/DE/dokumente/5231301237015_SB.xlsx",
                          cache = TRUE,
                          force = FALSE,
                          quiet = TRUE) {
  cache_key <- paste0("dt::", url)
  if (cache && !force && exists(cache_key, envir = .destatis_env, inherits = FALSE))
    return(get(cache_key, envir = .destatis_env, inherits = FALSE))

  tmp <- tempfile(fileext = ".xlsx")
  on.exit(unlink(tmp), add = TRUE)
  utils::download.file(url, tmp, mode = "wb", quiet = quiet)

  # Single source of truth, shared with age_to_bin_index() (see 00_utils.R).
  # Sheet 23131-01 is (ICD-10 four-character code) x (sex) x (total + 22 age
  # bands); the sex rows are aggregated below, the age bands are retained and
  # become the age-conditioned prior.
  age_labels <- .MICCI_AGE_BANDS

  raw <- readxl::read_excel(tmp, sheet = "23131-01", col_names = FALSE)
  dat <- as.data.frame(raw[-(seq_len(4L)), ])
  colnames(dat) <- c("code", "sex_raw", "total", age_labels)
  # Keep only rows that look like 4-character ICD-10 codes (e.g. "E11.4").
  dat <- dat[grepl("^[A-Z]\\d{2}\\.", dat$code), , drop = FALSE]

  parse_num <- function(x) {
    x <- as.character(x)
    x[x %in% c("-", ".", "/", "NA")] <- NA
    x <- gsub(".", "", x, fixed = TRUE)
    x <- gsub(",", ".", x, fixed = TRUE)
    suppressWarnings(as.numeric(x))
  }

  dt <- data.table(code = toupper(trimws(dat$code)),
                   freq_total = parse_num(dat$total))
  for (ag in age_labels) dt[, (ag) := parse_num(dat[[ag]])]

  # Aggregate across the per-sex rows.
  dt <- dt[, lapply(.SD, sum, na.rm = TRUE),
           by = code,
           .SDcols = c("freq_total", age_labels)]
  dt[, code_nodot := gsub(".", "", code, fixed = TRUE)]
  dt[, code3      := substr(code_nodot, 1L, 3L)]
  setkey(dt, code_nodot)

  message(sprintf("Destatis loaded: %d codes (%d distinct three-char prefixes)",
                  uniqueN(dt$code_nodot), uniqueN(dt$code3)))

  if (cache) assign(cache_key, dt, envir = .destatis_env)
  dt
}

# =============================================================================
# Per-prefix lookup precomputation
# =============================================================================

#' Precompute every per-prefix Charlson mapping and probability distribution.
#'
#' This is the work that makes the batch strategies fast: for every
#' three-character ICD prefix observed in the frequency table, miCCI
#' enumerates the matching four/five-character children, classifies the
#' prefix relative to each Charlson group as `certain` (every child code
#' triggers the group), `possible` (some children do, some don't), or
#' `none`, and stores the marginal subcode distribution. Once this
#' environment exists, S1 to S5 read O(1) per code.
#'
#' @param freq_table A data.table with the same schema as the output of
#'   `load_destatis()`: columns `code`, `code_nodot`, `code3`, `freq_total`,
#'   and optional age-band columns. To plug in your own population
#'   frequencies, pass them here.
#' @param quan_map   Output of `load_quan_map()`.
#' @return An environment containing one entry per three-character prefix
#'   plus the meta entries `.group_names`, `.dep_map`, `.wt_map`,
#'   `.all_pats`, `.dt_destatis`, `.dep_lookup_dt`.
#' @export
precompute_lookups <- function(freq_table, quan_map) {
  stopifnot(is.data.table(freq_table) || is.data.frame(freq_table))
  if (!is.data.table(freq_table)) freq_table <- as.data.table(freq_table)
  needed <- c("code", "code_nodot", "code3", "freq_total")
  miss <- setdiff(needed, names(freq_table))
  if (length(miss))
    stop("freq_table is missing required columns: ", paste(miss, collapse = ", "))

  cache <- new.env(parent = emptyenv(), hash = TRUE)

  all_prefixes <- unique(freq_table$code3)
  n_pref <- length(all_prefixes)

  group_names <- names(quan_map)
  all_pats <- lapply(setNames(group_names, group_names),
                     function(gk) drop_dot(.extract_patterns(quan_map[[gk]])))

  dep_map <- lapply(setNames(group_names, group_names),
                    function(gk) unlist(quan_map[[gk]]$depends_on))
  wt_map  <- vapply(group_names,
                    function(gk) as.numeric(quan_map[[gk]]$weight), numeric(1L))

  for (pref in all_prefixes) {
    children <- freq_table[code3 == pref, code_nodot]
    if (length(children) == 0L) next

    group_status <- setNames(rep("none", length(group_names)), group_names)
    child_groups <- setNames(vector("list", length(children)), children)

    for (gk in group_names) {
      pats_gk <- all_pats[[gk]]
      if (length(pats_gk) == 0L) next
      child_matches <- vapply(children, function(ch) {
        any(startsWith(ch, pats_gk) | startsWith(pats_gk, ch))
      }, logical(1L))

      n_match <- sum(child_matches)
      if (n_match == length(children)) {
        group_status[gk] <- "certain"
      } else if (n_match > 0L) {
        group_status[gk] <- "possible"
      }
      for (ch in children[child_matches]) {
        child_groups[[ch]] <- c(child_groups[[ch]], gk)
      }
    }

    child_freqs <- freq_table[code_nodot %in% children, .(code_nodot, freq_total)]
    total <- sum(child_freqs$freq_total, na.rm = TRUE)
    child_freqs[, prob := if (total > 0) freq_total / total else 0]

    # Age-conditioned prior. One column per Destatis age band, each column
    # normalised over the children of this prefix, so column a holds
    # P(subcode | prefix, age band a). A band with no recorded case for this
    # prefix yields an all-zero column; `age_ok` marks those so callers fall
    # back to the marginal instead of dividing by zero.
    child_prob_age <- NULL
    age_ok <- rep(FALSE, length(.MICCI_AGE_BANDS))
    if (all(.MICCI_AGE_BANDS %in% names(freq_table))) {
      am <- as.matrix(
        freq_table[match(child_freqs$code_nodot, code_nodot),
                   .SD, .SDcols = .MICCI_AGE_BANDS]
      )
      am[is.na(am)] <- 0
      col_tot <- colSums(am)
      age_ok  <- col_tot > 0
      if (any(age_ok)) {
        am[, age_ok] <- sweep(am[, age_ok, drop = FALSE], 2L,
                              col_tot[age_ok], "/")
      }
      am[, !age_ok] <- 0
      dimnames(am) <- list(child_freqs$code_nodot, .MICCI_AGE_BANDS)
      child_prob_age <- am
    }

    assign(pref,
           list(group_status   = group_status,
                child_freqs    = child_freqs,
                child_groups   = child_groups,
                child_prob_age = child_prob_age,
                age_ok         = age_ok,
                n_children     = length(children)),
           envir = cache)
  }

  assign(".group_names", group_names, envir = cache)
  assign(".dep_map",     dep_map,     envir = cache)
  assign(".wt_map",      wt_map,      envir = cache)
  assign(".all_pats",    all_pats,    envir = cache)
  assign(".dt_destatis", freq_table,  envir = cache)
  assign(".dep_lookup_dt", build_dep_lookup(quan_map), envir = cache)

  message(sprintf("Precomputed lookups for %d ICD-3 prefixes", n_pref))
  cache
}

#' Retrieve the cached lookup entry for a single three-character prefix.
#' Returns NULL if the prefix is not in the cache.
#' @keywords internal
get_prefix_cache <- function(prefix_3char, cache) {
  pref <- substr(drop_dot(prefix_3char), 1L, 3L)
  if (exists(pref, envir = cache, inherits = FALSE))
    get(pref, envir = cache, inherits = FALSE)
  else
    NULL
}

# =============================================================================
# Donor pools: the single place where "what can this prefix stand for?" is
# answered. Previously each strategy answered it separately and they
# disagreed: S1 and S2 dropped a prefix that had no children in the reference
# table, while S3 and S4 fell back to the prefix itself. That divergence is
# exactly why AIDS (B20-B24, which carry no four-character children in Destatis)
# scored zero under S1 and S2 but scored correctly under S3 and S4. It is also
# the undefined behaviour when every subcode of a prefix falls in the silent
# gap, since a prefix whose children are all absent from Destatis is precisely
# a prefix with no rows in Destatis.
#
# One helper now owns the decision, with three explicit outcomes:
#
#   resolved   - the reference table lists children; draw from them.
#   degenerate - the reference table lists no children, so the three-character
#                prefix stands for itself with probability 1. This is not an
#                error and not a zero: a three-character ICD-10-GM code is
#                itself codeable and scoreable, and dropping it would silently
#                discard a real diagnosis.
#   empty      - no codeable content at all (blank input). Contributes nothing.
# =============================================================================

#' Donor pool for one three-character prefix.
#'
#' @param prefix_3char three-character ICD prefix.
#' @param cache   output of `precompute_lookups()`.
#' @param age_idx optional 1-based Destatis age-band index (see
#'   `age_to_bin_index()`). When supplied and the band carries mass for this
#'   prefix, the age-conditioned probabilities are returned; otherwise the
#'   marginal is used and the fallback is flagged in the `age_fallback`
#'   attribute so the pipeline can count how often it happened.
#' @return data.table with `code_nodot` and `prob`, carrying attributes
#'   `status` and `age_fallback`.
#' @export
prefix_pool <- function(prefix_3char, cache, age_idx = NULL) {
  pref <- substr(drop_dot(prefix_3char), 1L, 3L)
  out_degenerate <- function() {
    d <- data.table(code_nodot = pref, prob = 1)
    setattr(d, "status", "degenerate")
    setattr(d, "age_fallback", FALSE)
    d
  }
  if (nchar(pref) < 3L) {
    d <- data.table(code_nodot = character(0L), prob = numeric(0L))
    setattr(d, "status", "empty"); setattr(d, "age_fallback", FALSE)
    return(d)
  }

  pc <- get_prefix_cache(pref, cache)
  if (is.null(pc)) return(out_degenerate())

  marg <- pc$child_freqs[prob > 0, .(code_nodot, prob)]
  if (nrow(marg) == 0L) return(out_degenerate())

  fallback <- FALSE
  if (!is.null(age_idx) && !is.na(age_idx) && !is.null(pc$child_prob_age) &&
      isTRUE(pc$age_ok[age_idx])) {
    # drop = FALSE and an explicit rownames lookup. Subsetting a column of a
    # ONE-ROW matrix silently discards the dimnames, so relying on names() here
    # would hand back a zero-length code vector beside a length-one
    # probability for every prefix that has a single child in the reference
    # table. That is a large share of prefixes, and the result is a malformed
    # donor pool rather than an error.
    p  <- as.numeric(pc$child_prob_age[, age_idx, drop = FALSE])
    cn <- rownames(pc$child_prob_age)
    keep <- p > 0 & !is.na(p)
    if (any(keep) && length(cn) == length(p)) {
      d <- data.table(code_nodot = cn[keep], prob = p[keep])
      setattr(d, "status", "resolved"); setattr(d, "age_fallback", FALSE)
      return(d)
    }
    fallback <- TRUE
  } else if (!is.null(age_idx) && !is.na(age_idx)) {
    fallback <- TRUE
  }

  d <- copy(marg)
  setattr(d, "status", "resolved")
  setattr(d, "age_fallback", fallback)
  d
}

#' Charlson group status of a degenerate prefix.
#'
#' A degenerate prefix has exactly one candidate code, itself. So a group is
#' either certain (the prefix matches it) or absent. There is no `possible`
#' state, because there is nothing to be uncertain between.
#' @keywords internal
.degenerate_status <- function(pref, group_names, all_pats) {
  st <- setNames(rep("none", length(group_names)), group_names)
  p0 <- drop_dot(pref)
  for (gk in group_names) {
    pats <- all_pats[[gk]]
    if (length(pats) == 0L) next
    if (any(startsWith(p0, pats) | startsWith(pats, p0))) st[gk] <- "certain"
  }
  st
}

#' Per-subcode probabilities for a three-character prefix, optionally
#' age-stratified using the Destatis age-band columns.
#' @keywords internal
get_subcode_probs_fast <- function(prefix_3char, cache, age = NULL) {
  pc <- get_prefix_cache(prefix_3char, cache)
  if (is.null(pc)) return(data.table(code_nodot = character(0L), prob = numeric(0L)))

  if (is.null(age) || is.na(age)) return(pc$child_freqs[, .(code_nodot, prob)])

  age_bin <- age_to_bin(age)
  if (is.na(age_bin)) return(pc$child_freqs[, .(code_nodot, prob)])

  dt_d <- get(".dt_destatis", envir = cache)
  if (!age_bin %in% names(dt_d)) return(pc$child_freqs[, .(code_nodot, prob)])

  children <- pc$child_freqs$code_nodot
  d <- dt_d[code_nodot %in% children]
  agg <- d[, .(freq = sum(get(age_bin), na.rm = TRUE)), by = code_nodot]
  total <- sum(agg$freq, na.rm = TRUE)
  agg[, prob := if (total > 0) freq / total else 0]
  agg[, .(code_nodot, prob)]
}
