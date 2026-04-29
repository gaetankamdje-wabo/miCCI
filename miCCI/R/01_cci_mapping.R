# =============================================================================
# miCCI / 01_cci_mapping.R
# Charlson group mapping: pattern lookup, dependency suppression, gold CCI.
# Single-encounter and vectorised batch versions stay numerically identical.
# =============================================================================

#' Load the bundled Quan et al. enhanced ICD-10 mapping.
#'
#' Reads the JSON shipped under `inst/extdata/codes_quan_orig.json`. If the
#' package is not installed and the working directory is the package source,
#' falls back to `inst/extdata/...` so that the function works during
#' `devtools::load_all()`.
#'
#' @param path Optional path to an alternative JSON file with the same
#'   schema as the bundled mapping. Useful if a site keeps its own
#'   curated Charlson mapping.
#' @return A named list, one entry per Charlson group, each with
#'   `name`, `weight`, `codes` (block list with `condition` + `codes`)
#'   and optionally `depends_on`.
#' @export
load_quan_map <- function(path = NULL) {
  if (is.null(path)) {
    path <- system.file("extdata", "codes_quan_orig.json", package = "miCCI")
    if (!nzchar(path))
      path <- file.path(getwd(), "inst", "extdata", "codes_quan_orig.json")
  }
  if (!file.exists(path))
    stop("Quan map JSON not found: ", path)
  jsonlite::fromJSON(path, simplifyVector = FALSE)$cci_icd_quan_orig
}

# Internal: pattern strings from a single group definition.
.extract_patterns <- function(group_def) {
  pats <- character(0L)
  for (blk in group_def$codes) {
    if (identical(blk$condition, "any")) pats <- c(pats, unlist(blk$codes))
  }
  unique(pats)
}

#' Build a flat (pattern -> group, weight) lookup from a Quan map.
#' @keywords internal
build_pattern_lookup <- function(quan_map) {
  rows <- list()
  for (gk in names(quan_map)) {
    pats <- .extract_patterns(quan_map[[gk]])
    if (length(pats) == 0L) next
    rows[[length(rows) + 1L]] <- data.table(
      pat = drop_dot(pats),
      gk  = gk,
      w   = quan_map[[gk]]$weight
    )
  }
  rbindlist(rows)
}

#' Build a (child -> parent) dependency lookup from a Quan map.
#'
#' A row says: if `parent` is also active for the same encounter, the
#' weight of `child` must not be counted.
#' @keywords internal
build_dep_lookup <- function(quan_map) {
  rows <- list()
  for (gk in names(quan_map)) {
    deps <- unlist(quan_map[[gk]]$depends_on)
    if (length(deps) == 0L) next
    rows[[length(rows) + 1L]] <- data.table(child = gk, parent = deps)
  }
  if (length(rows) == 0L)
    return(data.table(child = character(0L), parent = character(0L)))
  rbindlist(rows)
}

#' Single-encounter gold-standard Charlson Comorbidity Index.
#'
#' Use this for unit tests, REPL exploration and one-off lookups. For a
#' cohort, call `compute_predictions()` (which uses the vectorised batch
#' path) or `cci_gold_batch()` directly. Length-1 inputs are routed to
#' the batch path internally so the two functions never disagree.
#'
#' @param icd_codes character vector of full ICD-10-GM codes for one
#'   encounter, dotted or not.
#' @param quan_map  output of `load_quan_map()`.
#' @return list with elements
#'   \describe{
#'     \item{cci}{integer Charlson sum after dependency suppression.}
#'     \item{active}{named numeric vector of contributing weights.}
#'   }
#' @export
cci_gold <- function(icd_codes, quan_map) {
  codes_nd <- unique(drop_dot(icd_codes))
  codes_nd <- codes_nd[nchar(codes_nd) >= 3L]
  if (length(codes_nd) == 0L)
    return(list(cci = 0L, active = numeric(0L)))

  # Identify triggered groups via bidirectional prefix match.
  triggered <- character(0L)
  for (gk in names(quan_map)) {
    pats <- drop_dot(.extract_patterns(quan_map[[gk]]))
    if (length(pats) == 0L) next
    hit <- FALSE
    for (cn in codes_nd) {
      if (any(startsWith(cn, pats) | startsWith(pats, cn))) { hit <- TRUE; break }
    }
    if (hit) triggered <- c(triggered, gk)
  }
  if (length(triggered) == 0L) return(list(cci = 0L, active = numeric(0L)))

  # Dependency suppression: child is dropped if any of its parents is also triggered.
  active <- triggered
  for (gk in triggered) {
    deps <- unlist(quan_map[[gk]]$depends_on)
    if (length(deps) > 0L && any(deps %in% triggered))
      active <- setdiff(active, gk)
  }
  weights <- vapply(active, function(gk) as.numeric(quan_map[[gk]]$weight), numeric(1L))
  names(weights) <- active
  list(cci = as.integer(sum(weights)), active = weights)
}

#' Vectorised batch gold CCI for a full cohort.
#'
#' Equivalent to `vapply(split_codes, cci_gold, ..., quan_map)` but two to
#' three orders of magnitude faster on cohorts of 100k+ encounters.
#'
#' @param dt data.table with a character column `diagnosen` containing
#'   pipe-separated ICD codes.
#' @param quan_map        output of `load_quan_map()`.
#' @param pattern_lookup  optional precomputed pattern lookup.
#' @param dep_lookup      optional precomputed dependency lookup.
#' @return integer vector, length `nrow(dt)`, with the gold CCI per encounter.
#' @export
cci_gold_batch <- function(dt, quan_map,
                           pattern_lookup = NULL,
                           dep_lookup     = NULL) {
  if (is.null(pattern_lookup)) pattern_lookup <- build_pattern_lookup(quan_map)
  if (is.null(dep_lookup))     dep_lookup     <- build_dep_lookup(quan_map)
  pl <- pattern_lookup
  dl <- dep_lookup

  # Weight lookup, one row per group.
  wt <- unique(pl[, .(gk, w)])
  setkey(wt, gk)

  # Step 1: explode diagnoses into long form (idx, cn).
  dx_list <- strsplit(as.character(dt$diagnosen), "\\|+")
  lens    <- lengths(dx_list)
  long <- data.table(
    idx      = rep.int(seq_len(nrow(dt)), lens),
    code_raw = unlist(dx_list, use.names = FALSE)
  )
  long[, cn := toupper(gsub("[^A-Z0-9]", "", code_raw))]
  long <- long[nchar(cn) >= 3L]

  if (nrow(long) == 0L) return(rep(0L, nrow(dt)))

  # Step 2: match unique codes to Charlson groups (bidirectional prefix).
  unique_codes <- unique(long$cn)
  pat_vec <- pl$pat
  gk_vec  <- pl$gk
  match_rows <- list()
  for (uc in unique_codes) {
    hits <- gk_vec[startsWith(uc, pat_vec) | startsWith(pat_vec, uc)]
    if (length(hits) > 0L) {
      match_rows[[length(match_rows) + 1L]] <-
        data.table(cn = uc, gk = unique(hits))
    }
  }
  if (length(match_rows) == 0L) return(rep(0L, nrow(dt)))
  code_group_map <- rbindlist(match_rows); setkey(code_group_map, cn)

  # Step 3: join code-to-group onto the long form, get unique (idx, gk) pairs.
  setkey(long, cn)
  hits <- code_group_map[long, on = "cn", nomatch = NULL, allow.cartesian = TRUE]
  triggered <- unique(hits[, .(idx, gk)])

  # Step 4: dependency suppression via anti-join (no fsetdiff column-order trap).
  if (nrow(dl) > 0L && nrow(triggered) > 0L) {
    suppress <- merge(triggered, dl,
                      by.x = "gk", by.y = "child",
                      allow.cartesian = TRUE)
    suppress <- merge(suppress, triggered,
                      by.x = c("idx", "parent"), by.y = c("idx", "gk"))
    if (nrow(suppress) > 0L) {
      suppress_keys <- unique(suppress[, .(idx, gk)])
      triggered <- triggered[!suppress_keys, on = c("idx", "gk")]
    }
  }

  # Step 5: sum weights per encounter.
  scored <- merge(triggered, wt, by = "gk")[, .(cci = sum(w)), by = idx]
  out <- rep(0L, nrow(dt))
  out[scored$idx] <- as.integer(scored$cci)
  out
}
