# =============================================================================
# miCCI — 00_utils.R
# Shared utilities: ICD normalization, prefix matching, age binning
# =============================================================================

#' Normalize ICD code: uppercase, comma->dot, remove invalid chars
normalize_icd <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub(",", ".", x)
  x <- gsub("[^A-Z0-9.]", "", x)
  x <- gsub("\\.+", ".", x)
  gsub("\\.$", "", x)
}

#' Remove dot from ICD code
drop_dot <- function(x) gsub("\\.", "", normalize_icd(x))

#' Bidirectional prefix match (dot-insensitive)
like_match <- function(code, pattern) {
  c0 <- drop_dot(code); p0 <- drop_dot(pattern)
  startsWith(c0, p0) || startsWith(p0, c0)
}

#' Check if code matches ANY pattern in a set
matches_any_pattern <- function(code, patterns) {
  c0 <- drop_dot(code)
  any(vapply(patterns, function(p) {
    p0 <- drop_dot(p)
    startsWith(c0, p0) || startsWith(p0, c0)
  }, logical(1)))
}

#' Extract 3-character prefix
icd3 <- function(x) substr(drop_dot(normalize_icd(x)), 1, 3)

#' Add dot to 4-digit nodot code: "E114" -> "E11.4"
add_dot4 <- function(x) {
  x <- toupper(trimws(x))
  is_plain <- grepl("^[A-Z][0-9]{3}$", x)
  x[is_plain] <- paste0(substr(x[is_plain], 1, 3), ".", substr(x[is_plain], 4, 4))
  x
}

#' Map numeric age to Destatis age bin
age_to_bin <- function(age) {
  if (is.null(age) || length(age) == 0 || is.na(age)) return(NA_character_)
  a <- as.numeric(age)
  if (is.na(a)) return(NA_character_)
  if (a < 1)   return("unter 1")
  if (a >= 95) return("95 u. \u00e4lter")
  cuts   <- c(1,5,10,15,18,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95)
  labels <- c("1 - 5","5 - 10","10 - 15","15-18","18-20","20 - 25",
              "25 - 30","30 - 35","35 - 40","40 - 45","45 - 50","50 - 55",
              "55 - 60","60 - 65","65 - 70","70 - 75","75 - 80","80 - 85",
              "85 - 90","90 - 95")
  for (i in seq_along(labels)) {
    if (a >= cuts[i] && a < cuts[i + 1]) return(labels[i])
  }
  NA_character_
}

#' Split pipe-separated ICD codes (handles "|" and "||")
split_icd <- function(x) {
  if (is.na(x) || x == "" || x == "NA") return(character(0))
  codes <- unlist(strsplit(as.character(x), "\\|+"))
  codes <- normalize_icd(codes)
  codes[nchar(codes) >= 3]
}

#' Truncate ICD codes to 3-character level (simulate anonymisation)
truncate_icd <- function(codes) unique(icd3(codes))

#' Null-coalescing operator
`%||%` <- function(a, b) if (!is.null(a)) a else b
