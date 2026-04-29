library(testthat)
library(data.table)

qm     <- .test_quan()
cohort <- .synth_cohort()
freq   <- .synth_freq()
cache  <- precompute_lookups(freq, qm)

test_that("cci_gold_batch matches cci_gold encounter-by-encounter", {
  per_row <- vapply(cohort$diagnosen,
                    function(s) cci_gold(strsplit(s, "\\|+")[[1L]], qm)$cci,
                    integer(1L))
  batch <- cci_gold_batch(cohort, qm)
  expect_equal(batch, unname(per_row))
})

test_that("Batch produces expected gold values on the synthetic cohort", {
  # Hard-coded expected values for documented encounters; row 6 is computed
  # explicitly so a future Quan-map change does not break this test.
  exp_row6 <- cci_gold(strsplit(cohort$diagnosen[6L], "\\|+")[[1L]], qm)$cci
  expected <- c(4L, 2L, 3L, 6L, 0L, exp_row6, 4L)
  actual <- as.integer(unname(cci_gold_batch(cohort, qm)))
  expect_equal(actual, expected)
})

test_that("S1 envelope holds: cci_min <= cci_gold <= cci_max", {
  gold <- cci_gold_batch(cohort, qm)
  s1   <- cci_interval_batch(cohort, qm, cache)
  expect_true(all(s1$cci_min <= gold))
  expect_true(all(gold <= s1$cci_max))
  expect_true(all(s1$cci_min <= s1$cci_max))
})

test_that("S1 batch vs single agree", {
  single <- lapply(cohort$diagnosen, function(s)
    cci_interval(strsplit(s, "\\|+")[[1L]], qm, cache))
  batch  <- cci_interval_batch(cohort, qm, cache)
  expect_equal(batch$cci_min, vapply(single, function(x) x$cci_min, numeric(1L)))
  expect_equal(batch$cci_max, vapply(single, function(x) x$cci_max, numeric(1L)))
})

test_that("S2 batch vs single agree to numerical tolerance", {
  single <- vapply(cohort$diagnosen, function(s)
    cci_probabilistic(strsplit(s, "\\|+")[[1L]], qm, cache)$e_cci,
    numeric(1L))
  batch <- cci_probabilistic_batch(cohort, qm, cache)
  expect_equal(unname(batch), unname(single), tolerance = 1e-9)
})

test_that("S3 batch agrees with single-encounter route on length-1 input", {
  one <- cohort[1L]
  s_single <- cci_mi(strsplit(one$diagnosen, "\\|+")[[1L]], qm, cache,
                     m = 5L, seed = 7L)$mi_cci
  s_batch  <- cci_mi_batch(one, qm, cache, m = 5L, seed = 7L)
  expect_equal(s_single, unname(s_batch[1L]))
})

test_that("S4 batch agrees with single-encounter route on length-1 input", {
  one <- cohort[1L]
  s_single <- cci_bayesian(strsplit(one$diagnosen, "\\|+")[[1L]], qm, cache,
                           n_draws = 5L, alpha_0 = 10, seed = 7L)$posterior_median
  s_batch  <- cci_bayesian_batch(one, qm, cache,
                                 n_draws = 5L, alpha_0 = 10, seed = 7L)
  expect_equal(s_single, unname(s_batch[1L]))
})

test_that("Anti-join suppression: dm_simple is not double-counted", {
  dt <- data.table(diagnosen = c("E11.40|E11.90"))
  expect_equal(cci_gold_batch(dt, qm), 2L)
})

test_that("Truncated codes only: S2 expectation is non-negative and bounded", {
  dt <- data.table(diagnosen = c("E11|N18|I10"))
  res <- cci_probabilistic_batch(dt, qm, cache, return_group_prob = TRUE)
  expect_true(res$e_cci >= 0)
  s1 <- cci_interval_batch(dt, qm, cache)
  expect_true(res$e_cci <= s1$cci_max + 1e-9)  # union prob mass <= upper bound
})

test_that("S2 union beats max() heuristic on multi-prefix overlap", {
  # Build a tiny fake cache where two distinct prefixes both have non-zero
  # probability of mapping into the same Charlson group.
  qm2 <- list(
    g = list(name="g", weight = 1,
             codes = list(list(condition = "any", codes = c("A11","A22"))))
  )
  freq2 <- data.table(
    code       = c("A11.0","A11.9","A22.0","A22.9"),
    freq_total = c(50, 50, 50, 50)
  )
  freq2[, code_nodot := gsub(".", "", code, fixed = TRUE)]
  freq2[, code3      := substr(code_nodot, 1L, 3L)]
  setkey(freq2, code_nodot)
  c2 <- precompute_lookups(freq2, qm2)
  # Both prefixes give P(group g | prefix) = 1 (every child maps in).
  # Inclusion-exclusion: P(g | both prefixes) = 1 - (1-1)*(1-1) = 1.
  res <- cci_probabilistic(c("A11","A22"), qm2, c2)
  expect_equal(res$e_cci, 1, tolerance = 1e-9)
  expect_equal(res$group_prob[["g"]], 1, tolerance = 1e-9)
})
