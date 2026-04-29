library(testthat)
library(data.table)

skip_if_no_sl <- function() testthat::skip_if_not_installed("SuperLearner")

# Synthetic cohort with independent residuals (S3 best by design).
.make_synth <- function(n = 300L, seed = 1L) {
  set.seed(seed)
  gold <- sample(0:12, n, replace = TRUE, prob = stats::dpois(0:12, lambda = 2))
  s1_min <- pmax(gold - sample(0:3, n, TRUE, prob = c(0.4, 0.3, 0.2, 0.1)), 0)
  s1_max <- gold + sample(0:3, n, TRUE, prob = c(0.4, 0.3, 0.2, 0.1))
  s1_mid <- (s1_min + s1_max) / 2
  s2_ecci  <- gold + rnorm(n, 0, 0.9)
  s3_mi    <- gold + rnorm(n, 0, 0.4)   # best
  s4_bayes <- gold + rnorm(n, 0.05, 0.7)
  data.table(cci_gold = gold,
             s1_min = s1_min, s1_max = s1_max, s1_mid = s1_mid,
             s2_ecci = s2_ecci, s3_mi = s3_mi, s4_bayes = s4_bayes)
}

# Synthetic cohort with CORRELATED residuals: S3 and S4 share noise.
# This is the realistic situation - S3/S4 both draw from the same
# Destatis prior so their residuals are not independent. NNLS behaviour
# in this regime is the interesting test case.
.make_synth_correlated <- function(n = 300L, seed = 11L) {
  set.seed(seed)
  gold <- sample(0:12, n, replace = TRUE, prob = stats::dpois(0:12, lambda = 2))
  shared <- rnorm(n, 0, 0.4)            # shared component
  s1_min <- pmax(gold - 1, 0); s1_max <- gold + 1
  s1_mid <- (s1_min + s1_max) / 2
  s2_ecci  <- gold + rnorm(n, 0, 0.9)   # independent
  s3_mi    <- gold + shared + rnorm(n, 0, 0.2)   # correlated cluster
  s4_bayes <- gold + shared + rnorm(n, 0.05, 0.2)
  data.table(cci_gold = gold,
             s1_min = s1_min, s1_max = s1_max, s1_mid = s1_mid,
             s2_ecci = s2_ecci, s3_mi = s3_mi, s4_bayes = s4_bayes)
}

test_that("cci_meta_fit returns the expected list structure", {
  skip_if_no_sl()
  d <- .make_synth(200L)
  res <- cci_meta_fit(d, V = 5L, seed = 42L)
  expect_named(res, c("predictions", "weights", "cv_risk", "sl_fit"))
  expect_length(res$predictions, nrow(d))
  expect_length(res$weights, 4L)
  expect_named(res$weights, c("S1_mid", "S2_ecci", "S3_mi", "S4_bayes"))
  expect_length(res$cv_risk, 4L)
})

test_that("Meta weights are non-negative and finite", {
  skip_if_no_sl()
  d <- .make_synth(200L)
  res <- cci_meta_fit(d, V = 5L, seed = 42L)
  expect_true(all(res$weights >= -1e-8))
  expect_true(all(is.finite(res$weights)))
})

test_that("Predictions are numerically sane", {
  skip_if_no_sl()
  d <- .make_synth(300L)
  res <- cci_meta_fit(d, V = 5L, seed = 42L)
  expect_true(all(is.finite(res$predictions)))
  expect_true(all(res$predictions >= -1))
  expect_true(all(res$predictions <= 35))
  expect_lt(abs(mean(res$predictions) - mean(d$cci_gold)), 1.0)
})

test_that("Independent residuals: meta prefers the best learner", {
  skip_if_no_sl()
  d <- .make_synth(500L, seed = 7L)
  res <- cci_meta_fit(d, V = 5L, seed = 42L)
  expect_gt(res$weights["S3_mi"], res$weights["S2_ecci"])
  expect_gt(res$weights["S3_mi"], res$weights["S4_bayes"])
})

test_that("Meta MAE not worse than S3 by more than 20% (oracle property)", {
  skip_if_no_sl()
  d <- .make_synth(500L, seed = 7L)
  res <- cci_meta_fit(d, V = 10L, seed = 42L)
  mae <- function(x) mean(abs(x - d$cci_gold))
  expect_lt(mae(res$predictions), 1.20 * mae(d$s3_mi))
})

test_that("Correlated residuals: weights sum close to one and finite", {
  # S3 and S4 share a noise component - the test is whether NNLS still
  # behaves sanely (no negatives, predictions track gold). It does not
  # require the weight on a specific learner because under correlated
  # residuals the meta regression has multiple near-equivalent solutions.
  skip_if_no_sl()
  d <- .make_synth_correlated(500L, seed = 11L)
  res <- cci_meta_fit(d, V = 5L, seed = 42L)
  expect_true(all(res$weights >= -1e-8))
  # Sum-to-1 (not exact under NNLS, but very close in practice)
  expect_lt(abs(sum(res$weights) - 1), 0.05)
  # Beats the worst single learner
  mae <- function(x) mean(abs(x - d$cci_gold))
  expect_lt(mae(res$predictions),
            max(mae(d$s2_ecci), mae(d$s3_mi), mae(d$s4_bayes), mae(d$s1_mid)))
})

test_that("Missing column raises an informative error", {
  skip_if_no_sl()
  d <- .make_synth(50L); d[, s4_bayes := NULL]
  expect_error(cci_meta_fit(d, V = 5L), "s4_bayes")
})

test_that("Degenerate base learners collapse onto gold", {
  skip_if_no_sl()
  d <- .make_synth(200L, seed = 3L)
  for (col in c("s1_min","s1_max","s1_mid","s2_ecci","s3_mi","s4_bayes"))
    set(d, j = col, value = d$cci_gold)
  res <- cci_meta_fit(d, V = 5L, seed = 42L)
  expect_true(all(abs(res$predictions - d$cci_gold) < 1e-8))
})
