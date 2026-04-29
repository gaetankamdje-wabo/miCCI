# =============================================================================
# miCCI v0.5.0 — 60_cci_meta.R
# S6: Deterministic meta-estimator built on S1-S4 features.
#
# Rationale (manuscript-ready):
#   The meta-estimator combines the deterministic interval midpoint (S1) with
#   the consensus of the three distributional estimators (S2 expected, S3
#   multiple imputation, S4 Bayesian posterior median). The combination weight
#   is driven by the per-encounter S1 interval width: a narrow interval means
#   the truncated codes are diagnostic — S1 dominates; a wide interval means
#   ambiguity — S2/S3/S4 dominate.
#
# Definitions (per encounter):
#   m_dist  = mean(s2_ecci, s3_mi, s4_bayes)
#   w       = s1_width / (s1_width + kappa)        # 0 ≤ w < 1
#   raw     = (1 - w) * s1_mid + w * m_dist
#   s6_meta = clamp(raw, s1_min, s1_max)
#
# kappa = 2.0 is fixed by construction (twice the smallest non-trivial
#   Charlson weight). It is NOT learned from the data, so the estimator can
#   be applied to the entire 2010-2024 cohort without train/calibration split.
#
# Properties:
#   * Bounded by S1's certainty interval -> inherits S1 coverage.
#   * Reduces to s1_mid when s1_width = 0 (interval collapses).
#   * Reduces to mean(S2,S3,S4) (clipped) as s1_width -> infinity.
#   * Continuous, monotone in each input, fully training-free.
#   * Works on a single encounter or on a long vector (vectorised below).
# =============================================================================

#' Default kappa for the S6 meta-estimator
#' @export
MICCI_META_KAPPA <- 2.0

#' Compute S6 meta-estimator (deterministic, no training)
#'
#' @param s1_min   numeric vector — lower bound of S1 interval
#' @param s1_max   numeric vector — upper bound of S1 interval
#' @param s1_mid   numeric vector — midpoint of S1 interval
#' @param s2_ecci  numeric vector — S2 expected CCI
#' @param s3_mi    numeric vector — S3 multiple-imputation CCI
#' @param s4_bayes numeric vector — S4 Bayesian posterior median CCI
#' @param kappa    scalar smoothing constant; fixed at MICCI_META_KAPPA
#'
#' @return numeric vector of meta-estimates (float)
#' @export
cci_meta <- function(s1_min, s1_max, s1_mid,
                     s2_ecci, s3_mi, s4_bayes,
                     kappa = MICCI_META_KAPPA) {
  n <- length(s1_min)
  stopifnot(length(s1_max)  == n, length(s1_mid)  == n,
            length(s2_ecci) == n, length(s3_mi)   == n,
            length(s4_bayes) == n,
            is.numeric(kappa), length(kappa) == 1, kappa > 0)

  s1_min   <- as.numeric(s1_min)
  s1_max   <- as.numeric(s1_max)
  s1_mid   <- as.numeric(s1_mid)
  s2_ecci  <- as.numeric(s2_ecci)
  s3_mi    <- as.numeric(s3_mi)
  s4_bayes <- as.numeric(s4_bayes)

  s1_width <- pmax(s1_max - s1_min, 0)
  m_dist   <- (s2_ecci + s3_mi + s4_bayes) / 3

  w   <- s1_width / (s1_width + kappa)
  raw <- (1 - w) * s1_mid + w * m_dist

  pmin(pmax(raw, s1_min), s1_max)
}

#' Convenience wrapper that takes a data.table with the standard column names.
#' @param dt data.table with columns
#'   s1_min, s1_max, s1_mid, s2_ecci, s3_mi, s4_bayes
#' @export
cci_meta_dt <- function(dt, kappa = MICCI_META_KAPPA) {
  needed <- c("s1_min", "s1_max", "s1_mid", "s2_ecci", "s3_mi", "s4_bayes")
  miss <- setdiff(needed, names(dt))
  if (length(miss)) stop("cci_meta_dt: missing columns: ",
                         paste(miss, collapse = ", "))
  cci_meta(dt$s1_min, dt$s1_max, dt$s1_mid,
           dt$s2_ecci, dt$s3_mi, dt$s4_bayes,
           kappa = kappa)
}
