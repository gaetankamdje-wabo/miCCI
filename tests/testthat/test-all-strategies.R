library(testthat)
library(miCCI)
library(data.table)

# =============================================================================
# Reduced Quan-style map for end-to-end testing.
#
# This local fixture covers the Charlson groups exercised by the assertions
# below: dm_simple/complicated, kidney, liver_mild/severe, malignancy_meta/
# nonmeta, MI, HF, cerebrovascular, plegia, COPD, AIDS, rheumatic, dementia,
# peptic_ulcer, peripheral_vascular. Every public test in this file uses
# `test_quan` so the synthetic Charlson definitions are deterministic and
# do not depend on the package's shipped JSON.
# =============================================================================

test_quan <- list(
  dm_complicated = list(name = "DM complicated", weight = 2,
    codes = list(list(condition = "any", codes = c(
      "E10.2","E10.3","E10.4","E10.5","E10.7",
      "E11.2","E11.3","E11.4","E11.5","E11.7",
      "E12.2","E12.3","E12.4","E12.5","E12.7",
      "E13.2","E13.3","E13.4","E13.5","E13.7",
      "E14.2","E14.3","E14.4","E14.5","E14.7")))),
  dm_simple = list(name = "DM simple", weight = 1,
    codes = list(list(condition = "any", codes = c(
      "E10.0","E10.1","E10.6","E10.8","E10.9",
      "E11.0","E11.1","E11.6","E11.8","E11.9",
      "E12.0","E12.1","E12.6","E12.8","E12.9",
      "E13.0","E13.1","E13.6","E13.8","E13.9",
      "E14.0","E14.1","E14.6","E14.8","E14.9"))),
    depends_on = list("dm_complicated")),
  kidney = list(name = "Kidney disease", weight = 2,
    codes = list(list(condition = "any", codes = c(
      "I12.0","I13.1","N03.2","N03.3","N03.4","N03.5","N03.6","N03.7",
      "N05.2","N05.3","N05.4","N05.5","N05.6","N05.7",
      "N18","N19","N25.0","Z49.0","Z49.1","Z49.2","Z94.0","Z99.2")))),
  liver_mild = list(name = "Mild liver", weight = 1,
    codes = list(list(condition = "any", codes = c(
      "B18","K70.0","K70.1","K70.2","K70.3","K70.9",
      "K71.3","K71.4","K71.5","K71.6","K71.7",
      "K73","K74","K76.0","K76.2","K76.3","K76.4","K76.8","K76.9","Z94.4"))),
    depends_on = list("liver_severe")),
  liver_severe = list(name = "Severe liver", weight = 3,
    codes = list(list(condition = "any", codes = c(
      "I85.0","I85.9","I86.4","I98.2","K70.4","K71.1","K72.1","K72.9",
      "K76.5","K76.6","K76.7")))),
  malignancy_nonmeta = list(name = "Nonmetastatic malignancy", weight = 2,
    codes = list(list(condition = "any", codes = c(
      "C00","C01","C02","C03","C04","C16","C18","C20","C25",
      "C32","C33","C34","C50","C61","C71"))),
    depends_on = list("malignancy_meta")),
  malignancy_meta = list(name = "Metastatic malignancy", weight = 6,
    codes = list(list(condition = "any", codes = c("C77","C78","C79","C80")))),
  mi = list(name = "MI", weight = 1,
    codes = list(list(condition = "any", codes = c("I21","I22","I25.2")))),
  hf = list(name = "HF", weight = 1,
    codes = list(list(condition = "any", codes = c(
      "I09.9","I11.0","I13.0","I13.2","I25.5",
      "I42.0","I42.5","I42.6","I42.7","I42.8","I42.9","I43","I50","P29.0")))),
  cerebrovascular = list(name = "Cerebrovascular", weight = 1,
    codes = list(list(condition = "any", codes = c(
      "G45","G46","H34.0","I60","I61","I62","I63","I64",
      "I65","I66","I67","I68","I69")))),
  plegia = list(name = "Plegia", weight = 2,
    codes = list(list(condition = "any", codes = c(
      "G04.1","G11.4","G80.1","G80.2","G81","G82",
      "G83.0","G83.1","G83.2","G83.3","G83.4","G83.9")))),
  pulmo = list(name = "COPD", weight = 1,
    codes = list(list(condition = "any", codes = c(
      "I27.8","I27.9","J40","J41","J42","J43","J44","J45","J46","J47")))),
  aids = list(name = "AIDS", weight = 6,
    codes = list(list(condition = "any", codes = c("B20","B21","B22","B24")))),
  rheumatic = list(name = "Rheumatic", weight = 1,
    codes = list(list(condition = "any", codes = c(
      "M05","M06","M31.5","M32","M33","M34","M35.1","M35.3","M36.0")))),
  dementia = list(name = "Dementia", weight = 1,
    codes = list(list(condition = "any", codes = c(
      "F00","F01","F02","F03","F05.1","G30","G31.1")))),
  peptic_ulcer = list(name = "Peptic ulcer", weight = 1,
    codes = list(list(condition = "any", codes = c("K25","K26","K27","K28")))),
  peripheral_vascular = list(name = "PVD", weight = 1,
    codes = list(list(condition = "any", codes = c(
      "I70","I71","I73.1","I73.8","I73.9","I77.1","I79.0","I79.2"))))
)

# =============================================================================
# Gold-standard semantics
# =============================================================================

test_that("Example patient full codes = 4", {
  expect_equal(cci_gold(c("E11.40","N18.4","I10.00","E78.0"), test_quan)$cci, 4L)
})

test_that("depends_on: dm_simple suppressed when dm_complicated present", {
  res <- cci_gold(c("E11.40","E11.90"), test_quan)
  expect_equal(res$cci, 2L)
  expect_false("dm_simple" %in% names(res$active))
})

test_that("depends_on: liver severe suppresses mild", {
  expect_equal(cci_gold(c("K70.4","K74.6"), test_quan)$cci, 3L)
})

test_that("depends_on: meta suppresses nonmeta", {
  expect_equal(cci_gold(c("C34.1","C78.0"), test_quan)$cci, 6L)
})

test_that("Healthy patient = 0", {
  expect_equal(cci_gold(c("I10.00","E78.0","Z96.1"), test_quan)$cci, 0L)
})

test_that("Complex multimorbid = 10", {
  expect_equal(cci_gold(c("I21.0","I50.0","J44.10","E11.40","N18.4","I63.3","G81.1"),
                        test_quan)$cci, 10L)
})

test_that("Max CCI: every group active simultaneously sums their weights", {
  # One canonical code per group. We omit dm_simple, liver_mild and
  # malignancy_nonmeta because they are suppressed by their parents.
  codes <- c("E11.40",   # dm_complicated (2)
             "N18.4",    # kidney (2)
             "K70.4",    # liver_severe (3)
             "C78.0",    # malignancy_meta (6)
             "I21.0",    # mi (1)
             "I50.0",    # hf (1)
             "I63.3",    # cerebrovascular (1)
             "G81.1",    # plegia (2)
             "J44.10",   # pulmo (1)
             "B20.0",    # aids (6)
             "M05.0",    # rheumatic (1)
             "F00.1",    # dementia (1)
             "K25.0",    # peptic_ulcer (1)
             "I70.0")    # peripheral_vascular (1)
  expect_equal(cci_gold(codes, test_quan)$cci, 29L)
})

test_that("Only mild liver = 1", {
  expect_equal(cci_gold("K74.6", test_quan)$cci, 1L)
})

test_that("Only simple DM = 1", {
  expect_equal(cci_gold("E11.90", test_quan)$cci, 1L)
})

test_that("Only nonmeta cancer = 2", {
  expect_equal(cci_gold("C34.1", test_quan)$cci, 2L)
})

test_that("Truncated E11+N18+I10+E78 = 4", {
  expect_equal(cci_gold(c("E11","N18","I10","E78"), test_quan)$cci, 4L)
})

# =============================================================================
# Internal helpers - addressed via the `miCCI:::` operator because they are
# not exported by NAMESPACE (intentionally: they are implementation details).
# =============================================================================

test_that("normalize_icd", {
  expect_equal(miCCI:::normalize_icd("e11.40"), "E11.40")
  expect_equal(miCCI:::normalize_icd(" E11.4 "), "E11.4")
  expect_equal(miCCI:::normalize_icd("e11,40"), "E11.40")
})

test_that("like_match bidirectional", {
  expect_true(miCCI:::like_match("E11",   "E11.4"))
  expect_true(miCCI:::like_match("E11.4", "E11"))
  expect_false(miCCI:::like_match("E11", "E12"))
})

test_that("split_icd with single pipe", {
  expect_equal(miCCI:::split_icd("C16.3|C16.9|D63.0|E46"),
               c("C16.3","C16.9","D63.0","E46"))
})

test_that("truncate_icd", {
  expect_equal(sort(miCCI:::truncate_icd(c("E11.40","E11.90","N18.4"))),
               sort(c("E11","N18")))
})
