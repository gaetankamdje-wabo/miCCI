# Helpers shared by every test file.

# Compact synthetic Quan-style map. Uses a small but representative
# subset of Charlson groups so the suite runs in <2s.
.test_quan <- function() {
  list(
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
        "N18","N19","N25.0","Z49.0","Z49.1","Z49.2","Z94.0","Z99.2")))),
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
        "C00","C01","C02","C03","C04","C16","C18","C20","C25","C32","C33","C34","C50","C61","C71"))),
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
        "I27.8","I27.9","J40","J41","J42","J43","J44","J45","J46","J47")))),
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
        "I70","I71","I73.1","I73.8","I73.9","I77.1","I79.0","I79.2"))))
  )
}

# Build a synthetic frequency table covering enough codes to exercise
# every Charlson group in the test Quan map. Two children per parent
# prefix so "possible" classification has bite.
.synth_freq <- function() {
  full <- c(
    # DM complicated and simple, two children per E1x
    "E10.4","E10.9","E11.4","E11.9","E12.4","E12.9","E13.4","E13.9","E14.4","E14.9",
    # kidney
    "I12.0","I12.9","I13.1","I13.9","N18.1","N18.4","N18.9","N19.0","Z49.0","Z49.1",
    # liver mild + severe
    "B18.0","B18.9","K70.3","K70.9","K70.4","K70.0","K71.5","K71.1",
    "K74.6","K74.0","K76.5","K76.7",
    # malignancies
    "C16.0","C16.9","C18.4","C18.9","C50.1","C50.9","C71.0","C71.9",
    "C77.0","C77.9","C78.0","C78.9","C79.5","C79.9","C80.0","C80.9",
    # MI / HF
    "I21.0","I21.4","I22.0","I22.9","I25.2","I25.9",
    "I50.0","I50.9","I43.0","I42.9","P29.0","P29.8",
    # CVD / plegia
    "I63.3","I63.9","I60.0","I69.4","G81.0","G81.9","G82.5","G83.4",
    # COPD
    "J44.0","J44.9","J45.0","J45.9","J47.0","J47.9",
    # rheumatic
    "M05.0","M05.9","M06.0","M06.9","M32.1","M32.9",
    # dementia
    "F00.1","F00.9","F03.0","G30.0","G30.9",
    # PUD
    "K25.0","K25.9","K26.0","K27.9","K28.0","K28.9",
    # PVD
    "I70.0","I70.9","I71.3","I71.9","I73.1","I73.9",
    # AIDS
    "B20.0","B20.9","B21.0","B22.0","B24",
    # baseline non-Charlson, exercises "no group" branches
    "Z01.0","Z01.9","I10.0","I10.9","E78.0","E78.5"
  )
  d <- data.table::data.table(
    code       = full,
    freq_total = stats::setNames(seq_along(full), NULL) * 13L  # arbitrary positive
  )
  d[, code_nodot := gsub(".", "", code, fixed = TRUE)]
  d[, code3      := substr(code_nodot, 1L, 3L)]
  data.table::setkey(d, code_nodot)
  d
}

# Once-built shared cache used by multiple test files (no Destatis network call).
.test_cache <- function() {
  qm   <- .test_quan()
  freq <- .synth_freq()
  precompute_lookups(freq, qm)
}

# Tiny synthetic cohort for batch tests. Mix of full-length and truncated codes,
# plus a healthy patient.
.synth_cohort <- function() {
  data.table::data.table(
    falnr     = seq_len(7L),
    diagnosen = c(
      "E11.4|N18.4|I10.0",         # DMc + kidney + (none) -> 4
      "E11.9|E11.4",               # DMs suppressed by DMc -> 2
      "K70.4|K74.6",               # severe liver suppresses mild -> 3
      "C34.1|C78.0",               # meta suppresses nonmeta -> 6
      "I10.0|E78.0",               # no comorbidity -> 0
      "I21.0|I50.0|J44.0|E11.4|N18.4|I63.3|G81.1",  # complex multimorbid
      "E11|N18|I10|E78"            # truncated three-char codes
    )
  )
}
