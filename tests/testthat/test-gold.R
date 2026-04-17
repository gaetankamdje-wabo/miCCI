# =============================================================================
# miCCI — Unit tests (all strategies)
# =============================================================================
library(testthat)
library(data.table)

# --- Minimal Quan map for testing (no Destatis needed) ---
test_quan <- list(
  dm_complicated = list(name="DM complicated", weight=2,
    codes=list(list(condition="any", codes=c(
      "E10.2","E10.3","E10.4","E10.5","E10.7",
      "E11.2","E11.3","E11.4","E11.5","E11.7",
      "E12.2","E12.3","E12.4","E12.5","E12.7",
      "E13.2","E13.3","E13.4","E13.5","E13.7",
      "E14.2","E14.3","E14.4","E14.5","E14.7")))),
  dm_simple = list(name="DM simple", weight=1,
    codes=list(list(condition="any", codes=c(
      "E10.0","E10.1","E10.6","E10.8","E10.9",
      "E11.0","E11.1","E11.6","E11.8","E11.9",
      "E12.0","E12.1","E12.6","E12.8","E12.9",
      "E13.0","E13.1","E13.6","E13.8","E13.9",
      "E14.0","E14.1","E14.6","E14.8","E14.9"))),
    depends_on=list("dm_complicated")),
  kidney = list(name="Kidney disease", weight=2,
    codes=list(list(condition="any", codes=c(
      "I12.0","I13.1","N03.2","N03.3","N03.4","N03.5","N03.6","N03.7",
      "N05.2","N05.3","N05.4","N05.5","N05.6","N05.7",
      "N18","N19","N25.0","Z49.0","Z49.2","Z94.0","Z99.2")))),
  liver_mild = list(name="Mild liver", weight=1,
    codes=list(list(condition="any", codes=c(
      "B18","K70.0","K70.1","K70.2","K70.3","K70.9",
      "K71.3","K71.4","K71.5","K71.6","K71.7",
      "K73","K74","K76.0","K76.2","K76.3","K76.4","K76.8","K76.9","Z94.4"))),
    depends_on=list("liver_severe")),
  liver_severe = list(name="Severe liver", weight=3,
    codes=list(list(condition="any", codes=c(
      "I85.0","I85.9","I86.4","I98.2","K70.4","K71.1","K72.1","K72.9",
      "K76.5","K76.6","K76.7")))),
  malignancy_nonmeta = list(name="Nonmetastatic malignancy", weight=2,
    codes=list(list(condition="any", codes=c(
      "C00","C01","C02","C03","C04","C05","C06","C07","C08","C09",
      "C10","C11","C12","C13","C14","C15","C16","C17","C18","C19",
      "C20","C21","C22","C23","C24","C25","C26","C30","C31","C32",
      "C33","C34","C37","C38","C39","C40","C41","C43","C45","C46",
      "C47","C48","C49","C50","C51","C52","C53","C54","C55","C56",
      "C57","C58","C60","C61","C62","C63","C64","C65","C66","C67",
      "C68","C69","C70","C71","C72","C73","C74","C75","C76",
      "C81","C82","C83","C84","C85","C88","C90","C91","C92","C93",
      "C94","C95","C96","C97"))),
    depends_on=list("malignancy_meta")),
  malignancy_meta = list(name="Metastatic malignancy", weight=6,
    codes=list(list(condition="any", codes=c("C77","C78","C79","C80")))),
  mi = list(name="MI", weight=1,
    codes=list(list(condition="any", codes=c("I21","I22","I25.2")))),
  hf = list(name="HF", weight=1,
    codes=list(list(condition="any", codes=c(
      "I09.9","I11.0","I13.0","I13.2","I25.5",
      "I42.0","I42.5","I42.6","I42.7","I42.8","I42.9","I43","I50","P29.0")))),
  cerebrovascular = list(name="Cerebrovascular", weight=1,
    codes=list(list(condition="any", codes=c(
      "G45","G46","H34.0","I60","I61","I62","I63","I64","I65","I66","I67","I68","I69")))),
  plegia = list(name="Plegia", weight=2,
    codes=list(list(condition="any", codes=c(
      "G04.1","G11.4","G80.1","G80.2","G81","G82",
      "G83.0","G83.1","G83.2","G83.3","G83.4","G83.9")))),
  pulmo = list(name="COPD", weight=1,
    codes=list(list(condition="any", codes=c(
      "I27.8","I27.9","J40","J41","J42","J43","J44","J45","J46","J47",
      "J60","J61","J62","J63","J64","J65","J66","J67","J68.4","J70.1","J70.3")))),
  aids = list(name="AIDS", weight=6,
    codes=list(list(condition="any", codes=c("B20","B21","B22","B24")))),
  rheumatic = list(name="Rheumatic", weight=1,
    codes=list(list(condition="any", codes=c(
      "M05","M06","M31.5","M32","M33","M34","M35.1","M35.3","M36.0")))),
  dementia = list(name="Dementia", weight=1,
    codes=list(list(condition="any", codes=c("F00","F01","F02","F03","F05.1","G30","G31.1")))),
  peptic_ulcer = list(name="Peptic ulcer", weight=1,
    codes=list(list(condition="any", codes=c("K25","K26","K27","K28")))),
  peripheral_vascular = list(name="PVD", weight=1,
    codes=list(list(condition="any", codes=c(
      "I70","I71","I73.1","I73.8","I73.9","I77.1","I79.0","I79.2",
      "K55.1","K55.8","K55.9","Z95.8","Z95.9"))))
)

# ── Gold CCI ──
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
  expect_equal(cci_gold(c("I21.0","I50.0","J44.10","E11.40","N18.4","I63.3","G81.1"), test_quan)$cci, 10L)
})

test_that("Max CCI = 29", {
  codes <- c("I21.0","I50.0","I70.0","I63.0","F00.0","J44.0","M05.00",
             "K25.0","K70.4","E11.40","N18.5","G81.0","C78.0","B20")
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

# ── Truncated codes ──
test_that("Truncated E11+N18+I10+E78 = 4", {
  res <- cci_gold(c("E11","N18","I10","E78"), test_quan)
  expect_equal(res$cci, 4L)
  expect_true("dm_complicated" %in% names(res$active))
  expect_false("dm_simple" %in% names(res$active))
})

# ── Utilities ──
test_that("normalize_icd", {
  expect_equal(normalize_icd("e11.40"), "E11.40")
  expect_equal(normalize_icd(" E11.4 "), "E11.4")
  expect_equal(normalize_icd("e11,40"), "E11.40")
})

test_that("like_match bidirectional", {
  expect_true(like_match("E11", "E11.4"))
  expect_true(like_match("E11.4", "E11"))
  expect_false(like_match("E11", "E12"))
})

test_that("split_icd with single pipe", {
  expect_equal(split_icd("C16.3|C16.9|D63.0|E46"), c("C16.3","C16.9","D63.0","E46"))
  expect_equal(split_icd(NA), character(0))
  expect_equal(split_icd("NA"), character(0))
})

test_that("truncate_icd", {
  expect_equal(sort(truncate_icd(c("E11.40","E11.90","N18.4"))), sort(c("E11","N18")))
})
