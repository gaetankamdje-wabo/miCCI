library(testthat)

test_that("normalize_icd is idempotent and uppercases", {
  expect_equal(miCCI:::normalize_icd("e11.40"), "E11.40")
  expect_equal(miCCI:::normalize_icd(" E11.4 "), "E11.4")
  expect_equal(miCCI:::normalize_icd("e11,40"), "E11.40")
  expect_equal(miCCI:::normalize_icd("E11.4."), "E11.4")
})

test_that("like_match is bidirectional", {
  expect_true(miCCI:::like_match("E11", "E11.4"))
  expect_true(miCCI:::like_match("E11.4", "E11"))
  expect_false(miCCI:::like_match("E11", "E12"))
})

test_that("split_icd splits and drops short codes", {
  expect_equal(miCCI:::split_icd("C16.3|C16.9|D63.0|E46"),
               c("C16.3","C16.9","D63.0","E46"))
  expect_equal(miCCI:::split_icd(NA), character(0L))
  expect_equal(miCCI:::split_icd("NA"), character(0L))
  expect_equal(miCCI:::split_icd("AB"), character(0L))
})

test_that("truncate_icd reduces to three-char prefixes uniquely", {
  expect_equal(sort(miCCI:::truncate_icd(c("E11.40","E11.90","N18.4"))),
               sort(c("E11","N18")))
})

test_that("age_to_bin matches Destatis bins", {
  expect_equal(miCCI:::age_to_bin(0),   "unter 1")
  expect_equal(miCCI:::age_to_bin(0.5), "unter 1")
  expect_equal(miCCI:::age_to_bin(1),   "1 - 5")
  expect_equal(miCCI:::age_to_bin(4.9), "1 - 5")
  expect_equal(miCCI:::age_to_bin(5),   "5 - 10")
  expect_equal(miCCI:::age_to_bin(95),  "95 u. \u00e4lter")
  expect_equal(miCCI:::age_to_bin(120), "95 u. \u00e4lter")
  expect_true(is.na(miCCI:::age_to_bin(NA)))
  # Vectorised
  expect_equal(miCCI:::age_to_bin(c(0, 5, 95)),
               c("unter 1", "5 - 10", "95 u. \u00e4lter"))
})

test_that(".with_local_seed restores the caller RNG state", {
  set.seed(11)
  before <- runif(1L)
  set.seed(11)
  miCCI:::.with_local_seed(99, runif(1L))   # different inner stream
  after  <- runif(1L)
  # If global state were polluted, after != before.
  expect_equal(before, after)
})
