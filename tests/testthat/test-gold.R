# Gold-standard CCI semantics.
library(testthat)
library(data.table)

qm <- .test_quan()

test_that("Full codes: DMc + kidney = 4", {
  expect_equal(cci_gold(c("E11.40","N18.4","I10.00","E78.0"), qm)$cci, 4L)
})

test_that("Dependency: dm_simple suppressed by dm_complicated", {
  res <- cci_gold(c("E11.40","E11.90"), qm)
  expect_equal(res$cci, 2L)
  expect_false("dm_simple" %in% names(res$active))
})

test_that("Dependency: liver_severe suppresses liver_mild", {
  expect_equal(cci_gold(c("K70.4","K74.6"), qm)$cci, 3L)
})

test_that("Dependency: malignancy_meta suppresses malignancy_nonmeta", {
  expect_equal(cci_gold(c("C34.1","C78.0"), qm)$cci, 6L)
})

test_that("Healthy patient = 0", {
  expect_equal(cci_gold(c("I10.00","E78.0","Z96.1"), qm)$cci, 0L)
})

test_that("Complex multimorbid = 10", {
  expect_equal(
    cci_gold(c("I21.0","I50.0","J44.10","E11.40","N18.4","I63.3","G81.1"), qm)$cci,
    10L
  )
})

test_that("Truncated codes still trigger groups", {
  res <- cci_gold(c("E11","N18","I10","E78"), qm)
  expect_equal(res$cci, 4L)
  expect_true("dm_complicated" %in% names(res$active))
  expect_false("dm_simple"      %in% names(res$active))
})

test_that("Single mild liver code = 1", {
  expect_equal(cci_gold("K74.6", qm)$cci, 1L)
})

test_that("Z49.1 maps to kidney (after manual JSON correction)", {
  expect_equal(cci_gold("Z49.1", qm)$cci, 2L)
})
