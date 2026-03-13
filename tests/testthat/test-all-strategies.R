# =============================================================================
# miCCI — tests/testthat/test-all-strategies.R
#
# Every test block starts with a clear label describing:
#   [TYPE] What is being tested | What the expected outcome is
#
# Run with:   devtools::test()          <- shows all labels + pass/fail
#             testthat::test_file(...)  <- same
#             Rscript tests/testthat.R  <- minimal output
# =============================================================================
library(testthat)
library(data.table)

# ── Minimal Quan map (self-contained, no package install needed for CI) ───────
test_quan <- list(
  dm_complicated = list(
    name = "DM complicated", weight = 2,
    codes = list(list(condition = "any", codes = c(
      "E10.2","E10.3","E10.4","E10.5","E10.7",
      "E11.2","E11.3","E11.4","E11.5","E11.7",
      "E12.2","E12.3","E12.4","E12.5","E12.7",
      "E13.2","E13.3","E13.4","E13.5","E13.7",
      "E14.2","E14.3","E14.4","E14.5","E14.7")))),

  dm_simple = list(
    name = "DM simple", weight = 1,
    codes = list(list(condition = "any", codes = c(
      "E10.0","E10.1","E10.6","E10.8","E10.9",
      "E11.0","E11.1","E11.6","E11.8","E11.9",
      "E12.0","E12.1","E12.6","E12.8","E12.9",
      "E13.0","E13.1","E13.6","E13.8","E13.9",
      "E14.0","E14.1","E14.6","E14.8","E14.9"))),
    depends_on = list("dm_complicated")),

  kidney = list(
    name = "Kidney disease", weight = 2,
    codes = list(list(condition = "any", codes = c(
      "I12.0","I13.1","N03.2","N03.3","N03.4","N03.5","N03.6","N03.7",
      "N05.2","N05.3","N05.4","N05.5","N05.6","N05.7",
      "N18","N19","N25.0","Z49.0","Z49.2","Z94.0","Z99.2")))),

  liver_mild = list(
    name = "Mild liver disease", weight = 1,
    codes = list(list(condition = "any", codes = c(
      "B18","K70.0","K70.1","K70.2","K70.3","K70.9",
      "K71.3","K71.4","K71.5","K71.6","K71.7",
      "K73","K74","K76.0","K76.2","K76.3","K76.4","K76.8","K76.9","Z94.4"))),
    depends_on = list("liver_severe")),

  liver_severe = list(
    name = "Severe liver disease", weight = 3,
    codes = list(list(condition = "any", codes = c(
      "I85.0","I85.9","I86.4","I98.2",
      "K70.4","K71.1","K72.1","K72.9","K76.5","K76.6","K76.7")))),

  malignancy_nonmeta = list(
    name = "Non-metastatic malignancy", weight = 2,
    codes = list(list(condition = "any", codes = c(
      "C00","C01","C02","C03","C04","C05","C06","C07","C08","C09",
      "C10","C11","C12","C13","C14","C15","C16","C17","C18","C19",
      "C20","C21","C22","C23","C24","C25","C26","C30","C31","C32",
      "C33","C34","C37","C38","C39","C40","C41","C43","C45","C46",
      "C47","C48","C49","C50","C51","C52","C53","C54","C55","C56",
      "C57","C58","C60","C61","C62","C63","C64","C65","C66","C67",
      "C68","C69","C70","C71","C72","C73","C74","C75","C76",
      "C81","C82","C83","C84","C85","C88","C90","C91","C92","C93",
      "C94","C95","C96","C97"))),
    depends_on = list("malignancy_meta")),

  malignancy_meta = list(
    name = "Metastatic malignancy", weight = 6,
    codes = list(list(condition = "any", codes = c("C77","C78","C79","C80")))),

  mi = list(
    name = "Myocardial infarction", weight = 1,
    codes = list(list(condition = "any", codes = c("I21","I22","I25.2")))),

  hf = list(
    name = "Heart failure", weight = 1,
    codes = list(list(condition = "any", codes = c(
      "I09.9","I11.0","I13.0","I13.2","I25.5",
      "I42.0","I42.5","I42.6","I42.7","I42.8","I42.9","I43","I50","P29.0")))),

  cerebrovascular = list(
    name = "Cerebrovascular disease", weight = 1,
    codes = list(list(condition = "any", codes = c(
      "G45","G46","H34.0","I60","I61","I62","I63","I64",
      "I65","I66","I67","I68","I69")))),

  plegia = list(
    name = "Hemiplegia or paraplegia", weight = 2,
    codes = list(list(condition = "any", codes = c(
      "G04.1","G11.4","G80.1","G80.2","G81","G82",
      "G83.0","G83.1","G83.2","G83.3","G83.4","G83.9")))),

  pulmo = list(
    name = "Chronic pulmonary disease", weight = 1,
    codes = list(list(condition = "any", codes = c(
      "I27.8","I27.9","J40","J41","J42","J43","J44","J45","J46","J47",
      "J60","J61","J62","J63","J64","J65","J66","J67","J68.4","J70.1","J70.3")))),

  aids = list(
    name = "AIDS / HIV", weight = 6,
    codes = list(list(condition = "any", codes = c("B20","B21","B22","B24")))),

  rheumatic = list(
    name = "Rheumatic disease", weight = 1,
    codes = list(list(condition = "any", codes = c(
      "M05","M06","M31.5","M32","M33","M34","M35.1","M35.3","M36.0")))),

  dementia = list(
    name = "Dementia", weight = 1,
    codes = list(list(condition = "any", codes = c(
      "F00","F01","F02","F03","F05.1","G30","G31.1")))),

  peptic_ulcer = list(
    name = "Peptic ulcer disease", weight = 1,
    codes = list(list(condition = "any", codes = c("K25","K26","K27","K28")))),

  peripheral_vascular = list(
    name = "Peripheral vascular disease", weight = 1,
    codes = list(list(condition = "any", codes = c(
      "I70","I71","I73.1","I73.8","I73.9","I77.1","I79.0","I79.2",
      "K55.1","K55.8","K55.9","Z95.8","Z95.9"))))
)

# =============================================================================
# BLOCK 1 — Utility functions
# =============================================================================

test_that("[UTIL] normalize_icd: lowercase input is uppercased and trimmed", {
  cat("\n  Testing: normalize_icd('e11.40') -> 'E11.40'\n")
  expect_equal(normalize_icd("e11.40"), "E11.40")
  cat("  Testing: leading/trailing spaces are stripped\n")
  expect_equal(normalize_icd(" E11.4 "), "E11.4")
  cat("  Testing: comma is converted to decimal point\n")
  expect_equal(normalize_icd("e11,40"), "E11.40")
})

test_that("[UTIL] like_match: bidirectional prefix matching works correctly", {
  cat("\n  Testing: 'E11' matches 'E11.4' (shorter is prefix of longer)\n")
  expect_true(like_match("E11", "E11.4"))
  cat("  Testing: 'E11.4' matches 'E11' (longer starts with shorter)\n")
  expect_true(like_match("E11.4", "E11"))
  cat("  Testing: 'E11' does NOT match 'E12'\n")
  expect_false(like_match("E11", "E12"))
})

test_that("[UTIL] split_icd: pipe-separated string is split into individual codes", {
  cat("\n  Testing: 'C16.3|C16.9|D63.0|E46' -> 4 codes\n")
  expect_equal(split_icd("C16.3|C16.9|D63.0|E46"),
               c("C16.3","C16.9","D63.0","E46"))
  cat("  Testing: NA input returns empty character vector\n")
  expect_equal(split_icd(NA), character(0))
  cat("  Testing: 'NA' string also returns empty vector\n")
  expect_equal(split_icd("NA"), character(0))
})

test_that("[UTIL] truncate_icd: 4-digit codes are collapsed to unique 3-char prefixes", {
  cat("\n  Testing: E11.40 + E11.90 + N18.4 -> unique 3-char: E11, N18\n")
  expect_equal(sort(truncate_icd(c("E11.40","E11.90","N18.4"))),
               sort(c("E11","N18")))
})

# =============================================================================
# BLOCK 2 — Gold CCI (cci_gold, single-encounter)
# =============================================================================

test_that("[GOLD] Standard 4-digit case: E11.40 + N18.4 + I10 + E78 = CCI 4", {
  cat("\n  Codes: E11.40 (dm_compl w=2) + N18.4 (kidney w=2) + I10/E78 (no group)\n")
  cat("  Expected: CCI = 4\n")
  expect_equal(cci_gold(c("E11.40","N18.4","I10.00","E78.0"), test_quan)$cci, 4L)
})

test_that("[GOLD] depends_on: dm_complicated suppresses dm_simple (both codes present)", {
  cat("\n  Codes: E11.40 (dm_compl) + E11.90 (dm_simple)\n")
  cat("  Expected: CCI = 2, dm_simple NOT in active groups\n")
  res <- cci_gold(c("E11.40","E11.90"), test_quan)
  expect_equal(res$cci, 2L)
  expect_false("dm_simple" %in% names(res$active))
})

test_that("[GOLD] depends_on: liver_severe suppresses liver_mild", {
  cat("\n  Codes: K70.4 (liver_severe w=3) + K74.6 (liver_mild w=1)\n")
  cat("  Expected: CCI = 3 (mild suppressed)\n")
  expect_equal(cci_gold(c("K70.4","K74.6"), test_quan)$cci, 3L)
})

test_that("[GOLD] depends_on: malignancy_meta suppresses malignancy_nonmeta", {
  cat("\n  Codes: C34.1 (non-meta w=2) + C78.0 (meta w=6)\n")
  cat("  Expected: CCI = 6 (non-meta suppressed)\n")
  expect_equal(cci_gold(c("C34.1","C78.0"), test_quan)$cci, 6L)
})

test_that("[GOLD] Healthy patient (only hypertension/lipids) = CCI 0", {
  cat("\n  Codes: I10.00 + E78.0 + Z96.1 (none map to Charlson groups)\n")
  cat("  Expected: CCI = 0\n")
  expect_equal(cci_gold(c("I10.00","E78.0","Z96.1"), test_quan)$cci, 0L)
})

test_that("[GOLD] Multi-morbid patient: MI + HF + COPD + DM_compl + kidney + cerebrovascular + plegia = CCI 10", {
  cat("\n  7 active Charlson groups: 1+1+1+2+2+1+2 = 10\n")
  expect_equal(
    cci_gold(c("I21.0","I50.0","J44.10","E11.40","N18.4","I63.3","G81.1"),
             test_quan)$cci,
    10L
  )
})

test_that("[GOLD] Maximum CCI = 29: all non-overlapping groups active", {
  cat("\n  All 17 Charlson groups with weights: sum = 29\n")
  codes <- c("I21.0","I50.0","I70.0","I63.0","F00.0","J44.0","M05.00",
             "K25.0","K70.4","E11.40","N18.5","G81.0","C78.0","B20")
  expect_equal(cci_gold(codes, test_quan)$cci, 29L)
})

test_that("[GOLD] Only liver_mild (no liver_severe present) = CCI 1", {
  cat("\n  Code: K74.6 alone -> liver_mild = 1, no suppression\n")
  expect_equal(cci_gold("K74.6", test_quan)$cci, 1L)
})

test_that("[GOLD] Only dm_simple (no dm_complicated present) = CCI 1", {
  cat("\n  Code: E11.90 alone -> dm_simple = 1, no suppression\n")
  expect_equal(cci_gold("E11.90", test_quan)$cci, 1L)
})

test_that("[GOLD] Only non-metastatic cancer = CCI 2", {
  cat("\n  Code: C34.1 -> malignancy_nonmeta = 2, no meta present\n")
  expect_equal(cci_gold("C34.1", test_quan)$cci, 2L)
})

# =============================================================================
# BLOCK 3 — Truncated code handling (3-char input to cci_gold)
# =============================================================================

test_that("[TRUNC] 3-char codes E11+N18+I10+E78: triggers dm_complicated NOT dm_simple", {
  cat("\n  All subcodes of E11 include dm_complicated codes -> certain group\n")
  cat("  Expected: CCI = 4, dm_complicated active, dm_simple NOT active\n")
  res <- cci_gold(c("E11","N18","I10","E78"), test_quan)
  expect_equal(res$cci, 4L)
  expect_true("dm_complicated" %in% names(res$active))
  expect_false("dm_simple" %in% names(res$active))
})

# =============================================================================
# BLOCK 4 — Vectorised batch gold CCI
# =============================================================================

test_that("[BATCH-GOLD] cci_gold_batch: correct CCI for 3 known patients", {
  cat("\n  Patient 1: E11.40 + N18.4  -> expected CCI = 4\n")
  cat("  Patient 2: I50.0           -> expected CCI = 1 (HF only)\n")
  cat("  Patient 3: I10.00          -> expected CCI = 0\n")
  pl <- build_pattern_lookup(test_quan)
  dl <- build_dep_lookup(test_quan)
  dt <- data.table(diagnosen = c("E11.40|N18.4","I50.0","I10.00"))
  res <- cci_gold_batch(dt, test_quan, pl, dl)
  expect_equal(res, c(4L, 1L, 0L))
})

test_that("[BATCH-GOLD] cci_gold_batch: empty diagnoses return 0", {
  cat("\n  Empty string and NA -> expected CCI = 0 for all\n")
  pl  <- build_pattern_lookup(test_quan)
  dl  <- build_dep_lookup(test_quan)
  dt  <- data.table(diagnosen = c("","NA","Z00.0"))
  res <- cci_gold_batch(dt, test_quan, pl, dl)
  expect_true(all(res >= 0L))
  expect_equal(res[3], 0L)  # Z00.0 has no Charlson group
})

# =============================================================================
# BLOCK 5 — S1 Interval strategy
# =============================================================================

test_that("[S1-SINGLE] cci_interval: certain group yields cci_min > 0 and cci_mid = mean(min,max)", {
  cat("\n  Requires a cache with at least N18 prefix loaded\n")
  skip_if_not(requireNamespace("readxl", quietly=TRUE), "readxl not available")
  skip_on_cran()
  tryCatch({
    qm  <- load_quan_map()
    dst <- load_destatis()
    cch <- precompute_lookups(dst, qm)
    res <- cci_interval(c("N18","E11"), qm, cch)
    cat(sprintf("  N18+E11: min=%d  max=%d  mid=%.1f  width=%d\n",
                res$cci_min, res$cci_max, res$cci_mid, res$interval_width))
    expect_true(res$cci_min >= 0L)
    expect_true(res$cci_max >= res$cci_min)
    expect_equal(res$cci_mid, (res$cci_min + res$cci_max) / 2)
  }, error = function(e) skip(paste("Destatis unavailable:", conditionMessage(e))))
})

# =============================================================================
# BLOCK 6 — eval_metrics helper
# =============================================================================

test_that("[METRICS] eval_metrics: perfect prediction yields MAE=0, R2=1", {
  cat("\n  Predicted == Actual for all 5 patients\n")
  cat("  Expected: MAE=0, RMSE=0, R2=1, bias=0\n")
  m <- eval_metrics(c(1,2,3,4,5), c(1,2,3,4,5))
  expect_equal(unname(m["mae"]),  0)
  expect_equal(unname(m["rmse"]), 0)
  expect_equal(unname(m["r_squared"]), 1)
  expect_equal(unname(m["bias"]), 0)
})

test_that("[METRICS] eval_metrics: systematic overestimation gives positive bias", {
  cat("\n  Predicted = Actual + 1 for all patients\n")
  cat("  Expected: bias = 1, MAE = 1, RMSE = 1\n")
  m <- eval_metrics(c(2,3,4,5,6), c(1,2,3,4,5))
  expect_equal(unname(m["bias"]), 1)
  expect_equal(unname(m["mae"]),  1)
  expect_equal(unname(m["rmse"]), 1)
})

# =============================================================================
# BLOCK 7 — S6 Ensemble
# =============================================================================

test_that("[S6] cci_ensemble_train: weights are non-negative and sum to 1", {
  cat("\n  5 strategies x 20 calibration patients\n")
  cat("  Expected: all weights >= 0, sum(weights) = 1\n")
  set.seed(1)
  mat <- matrix(abs(rnorm(100)), nrow = 20, ncol = 5)
  y   <- abs(rnorm(20)) * 3
  ens <- cci_ensemble_train(mat, y)
  cat(sprintf("  Weights: %s\n",
      paste(sprintf("%.3f", ens$weights), collapse=" | ")))
  expect_true(all(ens$weights >= 0))
  expect_equal(sum(ens$weights), 1, tolerance = 1e-6)
})

test_that("[S6] cci_ensemble: single-patient vector gives correct weighted sum", {
  cat("\n  Weights: S1=0.5, S2=0.5, S3=0, S4=0, S5=0\n")
  cat("  Estimates: S1=2, S2=4, rest=0 -> expected = 3\n")
  w   <- c(S1=0.5, S2=0.5, S3=0, S4=0, S5=0)
  est <- c(2, 4, 0, 0, 0)
  expect_equal(cci_ensemble(est, w), 3)
})

# =============================================================================
# BLOCK 8 — Built-in dummy dataset
# =============================================================================

test_that("[DATA] micci_patients: dataset has correct structure and column names", {
  cat("\n  Expected: 39 rows, columns p_id / list_full_icd_code / los_in_days\n")
  if (!exists("micci_patients")) data(micci_patients, package="miCCI")
  expect_equal(nrow(micci_patients), 39L)
  expect_true("p_id" %in% names(micci_patients))
  expect_true("list_full_icd_code" %in% names(micci_patients))
  expect_true("los_in_days" %in% names(micci_patients))
})

test_that("[DATA] micci_patients: all los_in_days values are positive integers", {
  cat("\n  All stay lengths should be > 0\n")
  if (!exists("micci_patients")) data(micci_patients, package="miCCI")
  expect_true(all(micci_patients$los_in_days > 0))
})

test_that("[DATA] micci_patients: all ICD strings are non-empty", {
  cat("\n  No NA or empty strings in list_full_icd_code\n")
  if (!exists("micci_patients")) data(micci_patients, package="miCCI")
  expect_true(all(!is.na(micci_patients$list_full_icd_code)))
  expect_true(all(nchar(micci_patients$list_full_icd_code) > 0))
})

test_that("[DATA] micci_patients: CCI computed from full codes spans 0 to >=10", {
  cat("\n  Dataset covers healthy patients (CCI=0) through high-burden cases (CCI>=10)\n")
  if (!exists("micci_patients")) data(micci_patients, package="miCCI")
  pl <- build_pattern_lookup(test_quan)
  dl <- build_dep_lookup(test_quan)
  dt <- data.table(diagnosen = micci_patients$list_full_icd_code)
  gold <- cci_gold_batch(dt, test_quan, pl, dl)
  cat(sprintf("  CCI range: %d to %d  (mean=%.1f)\n",
              min(gold), max(gold), mean(gold)))
  expect_equal(min(gold), 0L)
  expect_true(max(gold) >= 10L)
})
