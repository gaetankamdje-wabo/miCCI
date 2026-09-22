# =============================================================================
# miCCI / 00_utils.R
# Shared helpers: ICD normalisation, prefix matching, age binning.
# =============================================================================

# Destatis age-bin breakpoints. Lower bound inclusive, upper bound exclusive.
# The first bin "unter 1" is age < 1; the last bin "95 u. aelter" is age >= 95.
.MICCI_AGE_CUTS <- c(1, 5, 10, 15, 18, 20, 25, 30, 35, 40, 45,
                     50, 55, 60, 65, 70, 75, 80, 85, 90, 95)
.MICCI_AGE_LABELS <- c(
  "1 - 5", "5 - 10", "10 - 15", "15-18", "18-20", "20 - 25",
  "25 - 30", "30 - 35", "35 - 40", "40 - 45", "45 - 50",
  "50 - 55", "55 - 60", "60 - 65", "65 - 70", "70 - 75",
  "75 - 80", "80 - 85", "85 - 90", "90 - 95"
)

# The 22 Destatis 23131-01 age-band column headers, in sheet order. This is the
# single source of truth: the loader names its columns from it and the
# age-conditioned prior indexes its probability matrix by it, so the two can
# never drift apart.
.MICCI_AGE_BANDS <- c("unter 1", .MICCI_AGE_LABELS, "95 u. \u00e4lter")

#' Normalise an ICD code (uppercase, comma -> dot, strip non-coding chars).
#'
#' @param x character vector of raw ICD-10-GM codes.
#' @return Cleaned uppercase codes with a single optional dot.
#' @keywords internal
normalize_icd <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub(",", ".", x, fixed = TRUE)
  x <- gsub("[^A-Z0-9.]", "", x)
  x <- gsub("\\.+", ".", x)
  sub("\\.$", "", x)
}

#' Strip the dot from an ICD code, e.g. E11.40 -> E1140.
#' @keywords internal
drop_dot <- function(x) gsub(".", "", normalize_icd(x), fixed = TRUE)

#' Bidirectional prefix match between two ICD codes (dot-insensitive).
#' Returns TRUE if either code is a prefix of the other.
#' @keywords internal
like_match <- function(code, pattern) {
  c0 <- drop_dot(code); p0 <- drop_dot(pattern)
  startsWith(c0, p0) || startsWith(p0, c0)
}

#' Does `code` match any pattern in `patterns` (bidirectional prefix)?
#' @keywords internal
matches_any_pattern <- function(code, patterns) {
  c0 <- drop_dot(code)
  p0 <- drop_dot(patterns)
  any(startsWith(c0, p0) | startsWith(p0, c0))
}

#' First three characters of an ICD code (after normalisation).
#' @keywords internal
icd3 <- function(x) substr(drop_dot(normalize_icd(x)), 1, 3)

#' Insert the dot into a 4-character no-dot ICD code, e.g. E114 -> E11.4.
#' Codes that are already dotted, or shorter/longer than 4 characters, are
#' returned unchanged.
#' @keywords internal
add_dot4 <- function(x) {
  x <- toupper(trimws(x))
  is_plain <- grepl("^[A-Z][0-9]{3}$", x)
  x[is_plain] <- paste0(substr(x[is_plain], 1, 3), ".", substr(x[is_plain], 4, 4))
  x
}

#' Map a numeric age to the Destatis 23131-01 age-bin label.
#'
#' Vectorised. NA for missing age. Boundary convention: lower inclusive,
#' upper exclusive, except the open-ended last bin "95 u. aelter".
#' @keywords internal
age_to_bin <- function(age) {
  a <- suppressWarnings(as.numeric(age))
  out <- rep(NA_character_, length(a))
  if (length(a) == 0L) return(out)

  ok <- !is.na(a)
  # under 1
  out[ok & a < 1]  <- "unter 1"
  # 95 and older
  out[ok & a >= 95] <- "95 u. \u00e4lter"
  # interior bins via findInterval (no R-level loop)
  mid <- ok & a >= 1 & a < 95
  if (any(mid)) {
    idx <- findInterval(a[mid], .MICCI_AGE_CUTS, rightmost.closed = FALSE)
    # idx is in 1..(length(cuts)-1); maps to .MICCI_AGE_LABELS
    out[mid] <- .MICCI_AGE_LABELS[idx]
  }
  out
}

#' Map a numeric age to the 1-based index of its Destatis age band.
#'
#' Same boundary convention as `age_to_bin()`, but returns the column position
#' in `.MICCI_AGE_BANDS` so the age-conditioned prior can index its probability
#' matrix without a string join. `NA_integer_` for a missing or unparseable age,
#' which every caller treats as "use the marginal prior".
#' @keywords internal
age_to_bin_index <- function(age) {
  lbl <- age_to_bin(age)
  match(lbl, .MICCI_AGE_BANDS)
}

#' Split a pipe-separated ICD-code string into a normalised character vector.
#'
#' Handles single and double pipe separators ("|" and "||"). Codes shorter
#' than three characters after normalisation are dropped.
#' @keywords internal
split_icd <- function(x) {
  if (is.na(x) || x == "" || x == "NA") return(character(0L))
  codes <- unlist(strsplit(as.character(x), "\\|+"), use.names = FALSE)
  codes <- normalize_icd(codes)
  codes[nchar(codes) >= 3L]
}

#' Truncate ICD codes to their three-character prefix, mimicking the
#' privacy-preserving truncation used in shared statutory data.
#'
#' @param codes character vector of full-length ICD-10-GM codes.
#' @param preserve_multiplicity if TRUE (new default), one prefix slot
#'   is kept per original diagnosis position, so an encounter coded E11.20 and
#'   E11.74 yields `c("E11", "E11")`. If FALSE, the previous behaviour is
#'   restored and repeated prefixes collapse to one.
#'
#'   Multiplicity matters because a truncated release still carries one row per
#'   coded diagnosis. Collapsing them discards the information that the prefix
#'   was coded twice, which caps how many subcodes S3 and S4 can draw from the
#'   same prefix and systematically undercounts comorbidities arising from the
#'   same ICD-10 block.
#' @keywords internal
truncate_icd <- function(codes, preserve_multiplicity = TRUE) {
  p <- icd3(codes)
  p <- p[nchar(p) >= 3L]
  if (preserve_multiplicity) p else unique(p)
}

#' Round half away from zero, the convention used in clinical reporting.
#'
#' `base::round()` rounds halves to even, so 0.5 and 1.5 both land on 0 and 2
#' respectively. Applied to the half-integer output of S1 that asymmetry
#' concentrates mass on even bins and shows up in a histogram as a
#' distributional pattern that is really an artefact of the rounding rule.
#' This helper rounds halves consistently upward in magnitude instead, and is
#' the only rounding used for binning anywhere in the package.
#' @keywords internal
.round_half_up <- function(x) {
  sign(x) * floor(abs(x) + 0.5)
}

#' Is x an exact integer (within floating-point tolerance)?
#' @keywords internal
.is_int_valued <- function(x, tol = 1e-8) {
  !is.na(x) & abs(x - round(x)) < tol
}

#' Null-coalescing operator: returns `a` unless it is NULL.
#' @keywords internal
`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Run an expression with a temporarily set RNG state, restoring the
#' caller's `.Random.seed` afterwards. The expression is captured by name
#' (lazy promise), so this never accidentally evaluates twice.
#'
#' This is dependency-free; no global RNG pollution.
#' @keywords internal
.with_local_seed <- function(seed, expr) {
  ex <- substitute(expr)
  pf <- parent.frame()
  has_old <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (has_old) old <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (has_old) {
      assign(".Random.seed", old, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  eval(ex, envir = pf)
}
