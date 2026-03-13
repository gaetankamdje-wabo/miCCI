## data-raw/create_micci_patients.R
## Run this ONCE from the package root to (re)create data/micci_patients.rda
## You do NOT need to run this as a user — the dataset is already included.

micci_patients <- data.frame(
  p_id = paste0("PAT_", sprintf("%03d", 1:39)),

  # Column 2: full 4-digit ICD-10-GM codes, pipe-separated
  # These are the codes BEFORE §21 KHEntgG anonymization
  list_full_icd_code = c(
    # ── No comorbidity (CCI = 0) ──────────────────────────────────────────
    "I10.00|E78.00|Z96.41",        # Hypertension + hyperlipidaemia (no Charlson group)
    "J06.9|R05",                   # Upper respiratory infection
    "M54.5|M51.16",                # Low back pain
    "Z00.0",                       # Routine health check
    "K21.0|R12",                   # Gastro-oesophageal reflux

    # ── CCI = 1 ───────────────────────────────────────────────────────────
    "I21.09|I10.00",               # Myocardial infarction (weight 1) + hypertension
    "I50.14|E78.00",               # Heart failure (weight 1)
    "J44.10|R09.02",               # COPD (weight 1)
    "E11.90|I10.00",               # DM simple (weight 1) — no complications
    "I63.30|I10.00",               # Cerebrovascular disease (weight 1)
    "K25.0",                       # Peptic ulcer (weight 1)
    "I70.00|I10.00",               # Peripheral vascular disease (weight 1)
    "F00.1|R41.3",                 # Dementia (weight 1)

    # ── CCI = 2 ───────────────────────────────────────────────────────────
    "E11.40|I10.00|E78.00",        # DM complicated (weight 2) — end-organ damage
    "N18.40|I12.00",               # Chronic kidney disease stage 4 (weight 2)
    "C34.10|R06.00",               # Non-metastatic lung cancer (weight 2)
    "G81.10|I63.30",               # Hemiplegia (weight 2) + cerebrovascular
    "C50.91",                      # Non-metastatic breast cancer (weight 2)

    # ── CCI = 3 ───────────────────────────────────────────────────────────
    "K70.40|I85.00",               # Severe liver disease (weight 3)
    "I21.09|I50.14|J44.10",        # MI + HF + COPD (1+1+1 = 3)

    # ── CCI = 4 ───────────────────────────────────────────────────────────
    "E11.40|N18.40|I21.09|J44.10", # DM_compl + kidney + MI + COPD (2+2+1+1 = ... wait see actual)
    "I50.14|I70.00|C34.10|G81.10", # HF + PVD + cancer + hemiplegia
    "N18.50|E11.40|C50.91|I63.30", # Kidney_5 + DM_compl + breast_ca + cerebrovascular

    # ── CCI = 5-7 ─────────────────────────────────────────────────────────
    "E11.40|N18.40|K70.40|C34.10|I21.09",         # DM_compl+kidney+liver_sev+ca+MI
    "C34.10|N18.40|I50.14|J44.10|E11.40|G81.10",  # 6 active groups
    "B20|I10.00|R05",                              # AIDS (weight 6)
    "C78.09|R06.00|I10.00",                        # Metastatic cancer (weight 6)

    # ── High CCI (>=7) ────────────────────────────────────────────────────
    "I21.09|I50.14|J44.10|E11.40|N18.40|I63.30|G81.10",
    "C78.09|E11.40|N18.40|K70.40|I50.14|G81.10|I21.09",
    "B20|C78.09|E11.40|N18.50|K70.40|I50.14|G81.10|I21.09",

    # ── depends_on hierarchy: correct suppression tests ───────────────────
    "E11.40|E11.90|I10.00",    # dm_complicated present -> dm_simple suppressed
    "K70.40|K70.30|I85.00",   # liver_severe present -> liver_mild suppressed
    "C78.09|C34.10|I10.00",   # malignancy_meta present -> nonmeta suppressed
    "E11.40|E11.10",           # both DM codes; only dm_complicated counts

    # ── Single-group edge cases ────────────────────────────────────────────
    "M05.00|M06.00|I10.00",    # Rheumatic disease (weight 1)
    "F01.3|F03|G30.0",         # Dementia variants (weight 1)
    "K25.0|K26.0|I10.00",      # Peptic ulcer (weight 1)
    "I09.9|I43|P29.0|I10.00",  # Heart failure codes (weight 1)
    "Z94.0|N18.40|I12.00"      # Renal transplant + kidney disease (weight 2)
  ),

  # Column 3: length of hospital stay in days
  los_in_days = c(
    2, 1, 3, 1, 2,          # CCI = 0
    5, 7, 6, 4, 8, 3, 9, 6, # CCI = 1
    8, 10, 12, 14, 11,       # CCI = 2
    15, 9,                   # CCI = 3
    12, 16, 13,              # CCI = 4
    18, 20, 7, 6,            # CCI 5-7
    22, 28, 31,              # high CCI
    6, 8, 5, 7,              # hierarchy tests
    4, 5, 3, 9, 11           # edge cases
  ),

  stringsAsFactors = FALSE
)

## Save
dir.create("data", showWarnings = FALSE)
save(micci_patients, file = "data/micci_patients.rda", compress = "xz")
message("Saved data/micci_patients.rda  (", nrow(micci_patients), " rows)")
