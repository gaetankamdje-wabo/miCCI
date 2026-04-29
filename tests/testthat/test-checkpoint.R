library(testthat)
library(data.table)

test_that("Cohort fingerprint changes on row swap", {
  dt1 <- data.table(falnr = 1:5, diagnosen = c("E11","N18","I10","K70","B18"))
  dt2 <- data.table(falnr = 1:5, diagnosen = c("E11","N18","I10","K70","B19"))
  expect_false(identical(miCCI:::.cohort_fingerprint(dt1),
                         miCCI:::.cohort_fingerprint(dt2)))
})

test_that("Cohort fingerprint stable on identical input", {
  dt <- data.table(falnr = 1:100, diagnosen = paste0("E11.", 0:99))
  expect_equal(miCCI:::.cohort_fingerprint(dt),
               miCCI:::.cohort_fingerprint(copy(dt)))
})

test_that("Checkpoint round-trip via save/load matches original", {
  tmp <- tempfile("micci_cp_", fileext = "")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  fp <- "test_fingerprint"
  payload <- list(a = 1:10, b = data.table(x = letters, y = 1:26))
  miCCI:::.cp_save(tmp, "demo", fp, payload)

  loaded <- miCCI:::.cp_load(tmp, "demo", fp)
  expect_equal(loaded, payload)
})

test_that("Checkpoint invalidates on fingerprint mismatch", {
  tmp <- tempfile("micci_cp_", fileext = "")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  miCCI:::.cp_save(tmp, "demo", "fp_a", list(value = 42))
  expect_null(miCCI:::.cp_load(tmp, "demo", "fp_b"))
  expect_equal(miCCI:::.cp_load(tmp, "demo", "fp_a"), list(value = 42))
})

test_that("Missing checkpoint file returns NULL", {
  tmp <- tempfile("micci_cp_", fileext = "")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  expect_null(miCCI:::.cp_load(tmp, "nonexistent_stage", "any_fp"))
})

test_that("Corrupt checkpoint file returns NULL (no crash)", {
  tmp <- tempfile("micci_cp_", fileext = "")
  dir.create(file.path(tmp, "checkpoints"), recursive = TRUE)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  writeLines("not an RDS file", file.path(tmp, "checkpoints", "broken.rds"))
  expect_null(miCCI:::.cp_load(tmp, "broken", "any_fp"))
})
