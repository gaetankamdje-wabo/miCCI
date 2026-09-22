# =============================================================================
# test-revision.R
# Covers the behaviour introduced for the second revision, so `R CMD check`
# exercises it as well as `tools/verify.R` does. The harness is the fuller
# instrument (64 assertions, including a brute-force expectation check); this
# file carries the invariants that must never regress silently.
# =============================================================================

library(testthat)
library(data.table)

# Synthetic reference table, independent of Destatis.
#   E11.2 / E11.7 -> dm_complicated, E11.9 -> dm_simple  (parent + child, one prefix)
#   K70.3 -> liver_mild, K70.4 -> liver_severe           (parent + child, one prefix)
#   I10.0 -> single child, the case that loses dimnames on column subset
#   B20   -> no four-character children at all           (degenerate)
.rev_fixture <- function() {
  bands <- miCCI:::.MICCI_AGE_BANDS
  mk <- function(code, total, young, old) {
    r <- as.list(setNames(rep(0, length(bands)), bands))
    r[[bands[2L]]]  <- young
    r[[bands[20L]]] <- old
    c(list(code = code, freq_total = total), r)
  }
  freq <- rbindlist(list(
    mk("E11.2", 300,  10, 260), mk("E11.7", 100,  5,  80),
    mk("E11.9", 600, 200,  60), mk("K70.3", 200, 50,  50),
    mk("K70.4", 100,  25,  25), mk("I10.0", 500, 100, 100)))
  freq[, code_nodot := gsub(".", "", code, fixed = TRUE)]
  freq[, code3      := substr(code_nodot, 1L, 3L)]
  setkey(freq, code_nodot)
  # Reuse the suite-wide synthetic map: self-contained, no installed
  # package and no Destatis download needed.
  qm <- .test_quan()
  list(freq = freq, quan_map = qm, cache = precompute_lookups(freq, qm))
}

test_that("a prefix with no children in the reference table stands for itself", {
  f <- .rev_fixture()
  pool <- prefix_pool("B20", f$cache)
  expect_equal(attr(pool, "status"), "degenerate")
  expect_equal(pool$code_nodot, "B20")
  expect_equal(pool$prob, 1)
  expect_equal(attr(prefix_pool("E11", f$cache), "status"), "resolved")
})

test_that("AIDS scores identically under all four strategies", {
  f <- .rev_fixture()
  dt <- data.table(diagnosen = "B20")
  expect_equal(cci_gold_batch(dt, f$quan_map)[1L], 6L)
  expect_equal(cci_interval_batch(dt, f$quan_map, f$cache)$cci_mid[1L], 6)
  expect_equal(cci_probabilistic_batch(dt, f$quan_map, f$cache)[1L], 6)
  expect_equal(cci_mi_batch(dt, f$quan_map, f$cache, m = 5L, seed = 1L)[1L], 6)
  expect_equal(cci_bayesian_batch(dt, f$quan_map, f$cache,
                                  n_draws = 5L, seed = 1L)[1L], 6)
})

test_that("the S1 upper bound respects hierarchical exclusion", {
  f <- .rev_fixture()
  # K70 can reach liver_mild (1) and liver_severe (3). No single subcode scores
  # 4, so the reachable maximum is 3.
  expect_equal(cci_interval_batch(data.table(diagnosen = "K70"),
                                  f$quan_map, f$cache)$cci_max[1L], 3)
  # E11 can reach dm_simple (1) and dm_complicated (2). Maximum is 2, not 3.
  expect_equal(cci_interval_batch(data.table(diagnosen = "E11"),
                                  f$quan_map, f$cache)$cci_max[1L], 2)
})

test_that("the S1 envelope still contains the gold score", {
  f <- .rev_fixture()
  for (cc in list("E11.2", "E11.9", "K70.3", "K70.4", "I10.0",
                  c("E11.9", "K70.3"), c("E11.2", "K70.4"))) {
    g  <- cci_gold(cc, f$quan_map)$cci
    iv <- cci_interval(cc, f$quan_map, f$cache)
    expect_lte(iv$cci_min, g)
    expect_gte(iv$cci_max, g)
  }
})

test_that("S2 returns the true expectation over the subcode support", {
  f <- .rev_fixture()
  # Brute force: enumerate every subcode of E11, score it, weight by probability.
  pool <- prefix_pool("E11", f$cache)
  truth <- sum(vapply(seq_len(nrow(pool)), function(i)
    pool$prob[i] * cci_gold(miCCI:::add_dot4(pool$code_nodot[i]),
                            f$quan_map)$cci, numeric(1L)))
  got <- cci_probabilistic_batch(data.table(diagnosen = "E11"),
                                 f$quan_map, f$cache)[1L]
  expect_equal(got, truth, tolerance = 1e-9)
})

test_that("S2 group probabilities decompose the reported score", {
  f <- .rev_fixture()
  wm  <- get(".wt_map", envir = f$cache)
  res <- cci_probabilistic_batch(data.table(diagnosen = "E11|K70"),
                                 f$quan_map, f$cache, return_group_prob = TRUE)
  expect_equal(sum(res$group_prob$p * wm[res$group_prob$gk]),
               res$e_cci[1L], tolerance = 1e-9)
})

test_that("the age-conditioned prior shifts subcode probabilities", {
  f <- .rev_fixture()
  young <- age_to_bin_index(3); old <- age_to_bin_index(87)
  py <- prefix_pool("E11", f$cache, age_idx = young)
  po <- prefix_pool("E11", f$cache, age_idx = old)
  expect_gt(py$prob[py$code_nodot == "E119"], 0.9)   # uncomplicated dominates
  expect_gt(po$prob[po$code_nodot == "E112"], 0.6)   # complicated dominates
  expect_equal(sum(py$prob), 1, tolerance = 1e-9)
  expect_equal(sum(po$prob), 1, tolerance = 1e-9)
})

test_that("a single-child prefix keeps its code under the age prior", {
  f <- .rev_fixture()
  # Subsetting a column of a one-row matrix drops the dimnames, which silently
  # produced an empty code vector before this was guarded.
  p <- prefix_pool("I10", f$cache, age_idx = age_to_bin_index(3))
  expect_equal(nrow(p), 1L)
  expect_equal(p$code_nodot, "I100")
  expect_equal(p$prob, 1, tolerance = 1e-9)
})

test_that("an empty age band falls back to the marginal", {
  f <- .rev_fixture()
  p <- prefix_pool("E11", f$cache, age_idx = age_to_bin_index(30))
  expect_true(attr(p, "age_fallback"))
  expect_true(all(is.finite(p$prob)))
  expect_equal(sum(p$prob), 1, tolerance = 1e-9)
})

test_that("code multiplicity raises the expected score", {
  f <- .rev_fixture()
  one <- cci_probabilistic_batch(data.table(diagnosen = "E11"),
                                 f$quan_map, f$cache)[1L]
  two <- cci_probabilistic_batch(data.table(diagnosen = "E11|E11"),
                                 f$quan_map, f$cache)[1L]
  collapsed <- cci_probabilistic_batch(data.table(diagnosen = "E11|E11"),
                                       f$quan_map, f$cache,
                                       preserve_multiplicity = FALSE)[1L]
  expect_gt(two, one)
  expect_equal(collapsed, one, tolerance = 1e-12)
})

test_that("the cluster bootstrap widens intervals under within-cluster correlation", {
  set.seed(11)
  pid  <- rep(seq_len(200L), each = 5L)
  off  <- rep(stats::rnorm(200L, 0, 1), each = 5L)
  gold <- stats::rpois(1000L, 2)
  pred <- gold + off + stats::rnorm(1000L, 0, 0.1)
  naive <- bootstrap_metrics(pred, gold, B = 300L, seed = 5L)
  clust <- bootstrap_metrics(pred, gold, B = 300L, seed = 5L, cluster = pid)
  expect_gt(clust$mae_hi - clust$mae_lo, naive$mae_hi - naive$mae_lo)
  expect_equal(clust$n_clusters, 200L)
  expect_equal(naive$mae, clust$mae, tolerance = 1e-12)
})

test_that("the paired bootstrap honours a Bonferroni level", {
  set.seed(3)
  gold <- stats::rpois(500L, 2)
  a <- gold + stats::rnorm(500L, 0, 0.30)
  b <- gold + stats::rnorm(500L, 0, 0.34)
  p95  <- paired_bootstrap_diff(a, b, gold, B = 300L, seed = 2L)
  pbon <- paired_bootstrap_diff(a, b, gold, B = 300L, seed = 2L,
                                conf = 1 - 0.05 / 20)
  expect_gt(pbon$diff_hi - pbon$diff_lo, p95$diff_hi - p95$diff_lo)
  expect_false(paired_bootstrap_diff(a, a, gold, B = 200L, seed = 2L)$decisive)
})

test_that("halves round away from zero, not to even", {
  expect_equal(miCCI:::.round_half_up(c(0.5, 1.5, 2.5, -0.5)), c(1, 2, 3, -1))
  expect_equal(round(c(0.5, 2.5)), c(0, 2))   # the artefact being avoided
})

test_that("exact agreement counts half-integers instead of rounding them", {
  fake <- data.table(cci_gold = c(0L, 1L, 1L, 2L, 3L),
                     s1_mid   = c(0, 1.5, 1, 2, 2.5),
                     s4_bayes = c(0, 1,   1, 2, 3))
  ea <- qa_exact_agreement(fake, c("s1_mid", "s4_bayes"), max_value = 5L)
  expect_equal(ea$summary[strategy == "s4_bayes", pct_exact_agreement], 100)
  expect_equal(ea$summary[strategy == "s1_mid",   pct_half_integer], 40)
  expect_equal(ea$summary[strategy == "s1_mid",   n_exact_agreement], 3L)
  v1 <- ea$by_value[strategy == "s1_mid" & cci_value == "1"]
  expect_equal(v1$n_gold, 2L)
  expect_equal(v1$n_exact_match, 1L)
})
