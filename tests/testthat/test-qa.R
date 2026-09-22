library(testthat)
library(data.table)

qm    <- .test_quan()
freq  <- .synth_freq()
cache <- precompute_lookups(freq, qm)

test_that("qa_group_coverage returns one row per Charlson group", {
  cohort <- .synth_cohort()
  cohort[, cci_gold := cci_gold_batch(cohort, qm)]
  res2 <- cci_probabilistic_batch(cohort, qm, cache, return_group_prob = TRUE)
  cohort[, s2_ecci := res2$e_cci]
  setattr(cohort, "s2_group_prob", res2$group_prob)

  cov <- qa_group_coverage(cohort, qm)
  expect_equal(nrow(cov), length(qm))
  expect_setequal(cov$group, names(qm))
})

test_that("Mass conservation = 1 when strategy probabilities equal gold activations", {
  # If we feed qa_group_coverage a synthetic group_prob table where p = 1
  # for every (idx, gk) that is gold-active and absent otherwise, the
  # ratio should be exactly 1 for every active group.
  cohort <- .synth_cohort()
  cohort[, cci_gold := cci_gold_batch(cohort, qm)]
  cohort[, s2_ecci := 0]  # placeholder, irrelevant for the ratio test

  active <- gold_active_long(cohort, qm)
  active[, p := 1.0]
  setattr(cohort, "s2_group_prob", active)

  cov <- qa_group_coverage(cohort, qm)
  active_groups <- cov[n_active > 0]
  expect_true(all(abs(active_groups$ratio_s2 - 1) < 1e-9))
})

test_that("qa_group_coverage falls back to recomputing S2 when attribute is missing", {
  cohort <- .synth_cohort()
  cohort[, cci_gold := cci_gold_batch(cohort, qm)]
  cohort[, s2_ecci := cci_probabilistic_batch(cohort, qm, cache)]
  # NO setattr - the function should reconstruct via cache.

  cov <- qa_group_coverage(cohort, qm, cache = cache)
  active_groups <- cov[n_active > 0]
  # Mass values should be finite and non-negative.
  expect_true(all(is.finite(active_groups$s2_mass)))
  expect_true(all(active_groups$s2_mass >= 0))
  # S3/S4 should be NA since they have no fallback.
  expect_true(all(is.na(cov$s3_mass)))
  expect_true(all(is.na(cov$s4_mass)))
})

test_that("S2 mass stays bounded above gold mass for full-codes input", {
  # When the cohort uses full-length codes (no truncation), S2 should
  # behave like S2 on a known answer: the expected mass per group should
  # not exceed the gold mass dramatically.
  cohort <- data.table(diagnosen = c(
    "E11.4|N18.4|I10.0",
    "K70.4|K74.6",
    "C34.1|C78.0"
  ))
  cohort[, cci_gold := cci_gold_batch(cohort, qm)]
  res2 <- cci_probabilistic_batch(cohort, qm, cache, return_group_prob = TRUE)
  cohort[, s2_ecci := res2$e_cci]
  setattr(cohort, "s2_group_prob", res2$group_prob)

  cov <- qa_group_coverage(cohort, qm)
  active <- cov[n_active > 0]
  # All ratios should be in [0, 2] for sane synthetic data
  expect_true(all(active$ratio_s2 >= 0 & active$ratio_s2 <= 2))
})

test_that("qa_score_distribution returns frequency and KS tables of expected shape", {
  cohort <- .synth_cohort()
  cohort[, cci_gold := cci_gold_batch(cohort, qm)]
  cohort[, s2_ecci := cci_probabilistic_batch(cohort, qm, cache)]
  cohort[, s3_mi := cohort$s2_ecci + 0.1]    # fake values; just shape test
  cohort[, s4_bayes := cohort$s2_ecci - 0.05]
  cohort[, meta := cohort$s2_ecci]

  out <- qa_score_distribution(cohort)
  expect_true(all(c("frequency", "ks") %in% names(out)))
  expect_true(nrow(out$frequency) > 0)
  expect_true(all(c("strategy","ks_d","ks_p","mean_pred","mean_gold","mean_diff")
                  %in% names(out$ks)))
  expect_true(all(out$ks$ks_d >= 0 & out$ks$ks_d <= 1))
})

test_that("qa_posterior_coverage is monotone: pct_within_0_5 <= pct_within_1 <= pct_within_2", {
  cohort <- .synth_cohort()
  cohort[, cci_gold := cci_gold_batch(cohort, qm)]
  cohort[, s2_ecci  := cci_probabilistic_batch(cohort, qm, cache)]
  cohort[, s3_mi    := cohort$s2_ecci]
  cohort[, s4_bayes := cohort$s2_ecci]
  cohort[, meta     := cohort$s2_ecci]

  pc <- qa_posterior_coverage(cohort)
  expect_true(all(pc$pct_within_0_5 <= pc$pct_within_1 + 1e-9))
  expect_true(all(pc$pct_within_1   <= pc$pct_within_2 + 1e-9))
})
