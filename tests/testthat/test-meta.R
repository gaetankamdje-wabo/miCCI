# =============================================================================
# miCCI v0.5.0 — Unit tests for S6 deterministic meta-estimator
# =============================================================================
library(testthat)

test_that("S6: collapsed interval returns s1_mid exactly", {
  # s1_width = 0  =>  w = 0  =>  raw = s1_mid  =>  clip = s1_mid
  out <- cci_meta(s1_min   = c(2, 5, 0),
                  s1_max   = c(2, 5, 0),
                  s1_mid   = c(2, 5, 0),
                  s2_ecci  = c(7, 1, 9),
                  s3_mi    = c(8, 2, 9),
                  s4_bayes = c(9, 0, 9))
  expect_equal(out, c(2, 5, 0))
})

test_that("S6: result is always inside the S1 interval", {
  set.seed(1)
  n <- 500
  lo <- runif(n, 0, 10)
  wd <- runif(n, 0, 8)
  hi <- lo + wd
  mid <- (lo + hi) / 2
  s2 <- runif(n, 0, 20)   # deliberately allowed to escape interval
  s3 <- runif(n, 0, 20)
  s4 <- runif(n, 0, 20)
  out <- cci_meta(lo, hi, mid, s2, s3, s4)
  expect_true(all(out >= lo - 1e-12))
  expect_true(all(out <= hi + 1e-12))
})

test_that("S6: monotone in m_dist when interval is wide", {
  # very wide interval => weight on m_dist is large
  base <- list(s1_min = 0, s1_max = 30, s1_mid = 15)
  a <- cci_meta(base$s1_min, base$s1_max, base$s1_mid, 5, 5, 5)
  b <- cci_meta(base$s1_min, base$s1_max, base$s1_mid, 10, 10, 10)
  c_ <- cci_meta(base$s1_min, base$s1_max, base$s1_mid, 20, 20, 20)
  expect_lt(a, b); expect_lt(b, c_)
})

test_that("S6: kappa controls weight allocation", {
  # smaller kappa -> trusts m_dist more for the same width
  out_k1 <- cci_meta(0, 4, 2, 4, 4, 4, kappa = 1)
  out_k2 <- cci_meta(0, 4, 2, 4, 4, 4, kappa = 2)
  out_k8 <- cci_meta(0, 4, 2, 4, 4, 4, kappa = 8)
  # m_dist (=4) > s1_mid (=2): smaller kappa -> closer to 4
  expect_gt(out_k1, out_k2)
  expect_gt(out_k2, out_k8)
  expect_true(all(c(out_k1, out_k2, out_k8) <= 4 + 1e-12))
  expect_true(all(c(out_k1, out_k2, out_k8) >= 2 - 1e-12))
})

test_that("S6: cci_meta_dt wrapper matches positional call", {
  d <- data.table::data.table(
    s1_min   = c(0, 1, 3),
    s1_max   = c(4, 5, 3),
    s1_mid   = c(2, 3, 3),
    s2_ecci  = c(2.5, 4, 3),
    s3_mi    = c(2.7, 4, 3),
    s4_bayes = c(2.9, 4, 3)
  )
  expect_equal(cci_meta_dt(d),
               cci_meta(d$s1_min, d$s1_max, d$s1_mid,
                        d$s2_ecci, d$s3_mi, d$s4_bayes))
})

test_that("S6: missing column triggers informative error", {
  d <- data.table::data.table(s1_min = 0, s1_max = 4, s1_mid = 2,
                              s2_ecci = 2, s3_mi = 2)   # no s4_bayes
  expect_error(cci_meta_dt(d), "s4_bayes")
})

test_that("S6: returns float, never coerced to integer", {
  out <- cci_meta(0, 4, 2, 2.7, 2.3, 2.5)
  expect_type(out, "double")
})

test_that("S6: handles zero-length input", {
  expect_equal(length(cci_meta(numeric(0), numeric(0), numeric(0),
                               numeric(0), numeric(0), numeric(0))), 0L)
})
