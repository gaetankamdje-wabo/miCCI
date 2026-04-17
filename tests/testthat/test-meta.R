# =============================================================================
# miCCI v0.5.0 — Unit tests for S6 Super Learner meta-estimator
# =============================================================================
library(testthat)
library(data.table)

skip_if_no_sl <- function() {
  testthat::skip_if_not_installed("SuperLearner")
}

# Helper: generate a plausible synthetic cohort where S3 is the best strategy.
.make_synth <- function(n = 300L, seed = 1L) {
  set.seed(seed)
  gold <- sample(0:12, n, replace = TRUE,
                 prob = dpois(0:12, lambda = 2))
  s1_min   <- pmax(gold - sample(0:3, n, TRUE, prob = c(0.4, 0.3, 0.2, 0.1)), 0)
  s1_max   <- gold + sample(0:3, n, TRUE, prob = c(0.4, 0.3, 0.2, 0.1))
  s1_mid   <- (s1_min + s1_max) / 2
  s2_ecci  <- gold + rnorm(n, 0,    0.9)
  s3_mi    <- gold + rnorm(n, 0,    0.4)   # best
  s4_bayes <- gold + rnorm(n, 0.05, 0.7)
  data.table(cci_gold = gold,
             s1_min = s1_min, s1_max = s1_max, s1_mid = s1_mid,
             s2_ecci = s2_ecci, s3_mi = s3_mi, s4_bayes = s4_bayes)
}

test_that("S6: cci_meta_fit returns the expected list structure", {
  skip_if_no_sl()
  d <- .make_synth(200)
  res <- cci_meta_fit(d, V = 5L, seed = 42L)
  expect_named(res, c("predictions", "weights", "cv_risk", "sl_fit"))
  expect_length(res$predictions, nrow(d))
  expect_length(res$weights, 4L)
  expect_named(res$weights, c("S1_mid", "S2_ecci", "S3_mi", "S4_bayes"))
  expect_length(res$cv_risk, 4L)
})

test_that("S6: Super Learner weights are non-negative and finite", {
  skip_if_no_sl()
  d <- .make_synth(200)
  res <- cci_meta_fit(d, V = 5L, seed = 42L)
  expect_true(all(res$weights >= -1e-8))
  expect_true(all(is.finite(res$weights)))
})

test_that("S6: predictions are numerically sane", {
  skip_if_no_sl()
  d <- .make_synth(300)
  res <- cci_meta_fit(d, V = 5L, seed = 42L)
  # All predictions finite
  expect_true(all(is.finite(res$predictions)))
  # All predictions within a reasonable range (Quan CCI max is 29)
  expect_true(all(res$predictions >= -1))
  expect_true(all(res$predictions <= 35))
  # Mean prediction should track the gold mean within 1 CCI point
  expect_lt(abs(mean(res$predictions) - mean(d$cci_gold)), 1.0)
})

test_that("S6: the Super Learner prefers the best base learner", {
  # With S3 being clearly the best (sd=0.4 vs 0.9 and 0.7 for S2 and S4)
  # the NNLS meta-weight on S3 should dominate.
  skip_if_no_sl()
  d <- .make_synth(500, seed = 7L)
  res <- cci_meta_fit(d, V = 5L, seed = 42L)
  expect_gt(res$weights["S3_mi"], res$weights["S2_ecci"])
  expect_gt(res$weights["S3_mi"], res$weights["S4_bayes"])
})

test_that("S6: S6 MAE is at most as bad as the best single base learner", {
  # Oracle property: asymptotically the Super Learner is at least as good
  # as the best single learner in the library. With n=500 we expect
  # MAE(S6) to be within a small tolerance of MAE(S3).
  skip_if_no_sl()
  d <- .make_synth(500, seed = 7L)
  res <- cci_meta_fit(d, V = 10L, seed = 42L)
  mae <- function(x) mean(abs(x - d$cci_gold))
  mae_s6 <- mae(res$predictions)
  mae_s3 <- mae(d$s3_mi)
  # S6 must not be drastically worse than S3 — allow 20% slack for finite-n
  expect_lt(mae_s6, 1.20 * mae_s3)
})

test_that("S6: missing column triggers informative error", {
  skip_if_no_sl()
  d <- .make_synth(50)
  d[, s4_bayes := NULL]
  expect_error(cci_meta_fit(d, V = 5L), "s4_bayes")
})

test_that("S6: degenerate base learners collapse to gold", {
  skip_if_no_sl()
  d <- .make_synth(200, seed = 3L)
  # Collapse all four base learners onto the gold value. Any convex
  # non-negative combination must then also equal the gold, so S6
  # must reproduce it exactly.
  d[, s1_min   := cci_gold]
  d[, s1_max   := cci_gold]
  d[, s1_mid   := cci_gold]
  d[, s2_ecci  := cci_gold]
  d[, s3_mi    := cci_gold]
  d[, s4_bayes := cci_gold]
  res <- cci_meta_fit(d, V = 5L, seed = 42L)
  expect_true(all(abs(res$predictions - d$cci_gold) < 1e-8))
})
