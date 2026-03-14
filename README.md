# miCCI: Charlson Comorbidity Index from Anonymized ICD-10-GM Codes

[![R package](https://img.shields.io/badge/R-%3E%3D%204.1.0-blue)](https://cran.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Version](https://img.shields.io/badge/version-0.4.3-green)](DESCRIPTION)

> **Rationale**: German hospitals routinely strip diagnosis codes down to three characters due to data privacy reasons (applying for example k-anonymity, l-diversity, etc.), before sharing data for research. That truncation makes it impossible to calculate the Charlson Comorbidity Index (CCI) accurately. miCCI provides approaches for fixing that. It offers six computational methods: from a simple rule-based estimate to a machine learning model- each producing a CCI score from those shortened codes with well-characterized accuracy.

---

## Table of Contents

1. [What problem does this solve?](#what-problem-does-this-solve)
2. [The six strategies at a glance](#the-six-strategies-at-a-glance)
3. [Installation](#installation)
4. [Quick start — five minutes to your first result](#quick-start)
5. [Using each strategy independently](#using-each-strategy-independently)
6. [Using all strategies together](#using-all-strategies-together)
7. [Built-in dummy dataset](#built-in-dummy-dataset)
8. [Running the full validation pipeline](#running-the-full-validation-pipeline)
9. [Running tests](#running-tests)
10. [Citation](#citation)
11. [License](#license)

---

## What problem does this solve?

The **Charlson Comorbidity Index** is the most widely used tool for quantifying how sick a patient is based on their diagnoses [1]. It assigns a numeric weight to 17 chronic conditions, so from 1 for myocardial infarction to 6 for AIDS or metastatic cancer; and sums them into a single score predicting mortality and resource use [2].

Computing CCI correctly requires the full-digit ICD code. For example, `E11.40` (diabetic nephropathy) clearly maps to the *complicated diabetes* group (weight 2). The code `E11.90` (simple type-2 diabetes, no complications) maps to *simple diabetes* (weight 1). The three-character prefix `E11` alone is ambiguous.

German law, however, requires hospitals to report full diagnostic detail for billing under § 21 KHEntgG, yet data privacy regulations for external researchers often mandate a reduction to 3-digit ICD codes [3]. For instance, for the 17 million annual inpatient cases, this may create an important granularity gap. This makes the accurate computation of the CCI difficult, and leads often to a systematic underestimation or overestimation of the true comorbid burden.

miCCI can aid resolving this by using advanced computational methods, to model the truncation probabilistically, using for example the national frequency data from the German Federal Statistical Office (Destatis 23131-01) as a prior distribution over which four-digit subcode was actually documented [4].

---

## The six strategies at a glance

| Strategy | Method | Needs Destatis? | Needs training data? | Speed | Best for |
|---|---|:---:|:---:|---|---|
| **S1** Interval | ICD hierarchy logic: This non-statistical approach calculates a CCI range by identifying "certain" versus "possible" Charlson groups for each 3-character prefix and uses the interval midpoint as the point estimate. | No | No | ~5 s / 500k | Uncertainty bounds, no external data |
| **S2** Probabilistic | Expected value: This strategy uses national frequency data to assign continuous probability weights to each prefix, producing a deterministic expected-value estimate in a single vectorized pass. | Yes | No | ~5 s / 500k | Fast, good approximation |
| **S3** MI-CCI | Multiple imputation: Treating truncation as a missing data problem, this method draws 20 independent 4-digit subcodes per prefix from the national distribution and averages the resulting gold-standard CCI values.| Yes | No | ~4 h / 500k | Best point accuracy |
| **S4** Bayesian | Dirichlet posterior: This approach samples subcodes from a posterior Dirichlet distribution that combines national priors with Jeffreys-type regularization, using the posterior median across 25 draws as the final estimate.| Yes | No | ~8 h / 500k | Downstream Bayesian analysis |
| **S5** ML-CCI | XGBoost: An XGBoost model predicts the gold-standard CCI by training on 3-character indicators and clinical features, followed by isotonic regression for post-hoc calibration. | Yes | **Yes** | ~75 s train | Orthogonal signal in ensemble |
| **S6** Ensemble | NNLS meta-learner: This strategy combines the predictions of S1 through S5 using a Non-Negative Least Squares meta-learner to assign optimal, non-negative weights that minimize estimation error.| — | **Yes (S1-S5)** | Instant | Maximum R² |

> **Validated results (n=144,660 prospective):** S3 MAE=0.136, R²=0.973 · S4 MAE=0.135 · S6 R²=0.974 · S1 interval coverage=95.1%

---

## Installation

### From GitHub (recommended)

```r
# Step 1 — install the devtools package if you do not have it
install.packages("devtools")

# Step 2 — install miCCI directly from GitHub
devtools::install_github("https://github.com/gaetankamdje-wabo/miCCI")
```

### From a local clone

```bash
# Step 1 — clone the repository to your machine
git clone https://github.com/gaetankamdje-wabo/miCCI.git

# Step 2 — open R and install from the local folder
```

```r
# Step 3 — install from the cloned folder (adjust the path if needed)
devtools::install("path/to/miCCI")
```

### Dependencies

All dependencies install automatically. The core ones are:

```r
# These are installed automatically with miCCI
install.packages(c("data.table", "jsonlite", "readxl",
                   "xgboost", "nnls", "ggplot2"))

# Optional — only needed for reading Parquet files in run_pipeline()
install.packages("arrow")
```

### Python users

miCCI is an R package. Python users can call it via the `rpy2` bridge:

```python

# ║  Cell 1 — run once to install everything         

# !apt-get install -y -q r-base
# !pip install -q rpy2 pandas matplotlib
# !git clone https://github.com/gaetankamdje-wabo/miCCI.git
# !Rscript -e "install.packages(c('remotes','data.table','readxl','nnls'), repos='https://cloud.r-project.org', quiet=TRUE)"
# !Rscript -e "remotes::install_local('miCCI', upgrade='never')"


# ║  Cell 2 — run the session                        


import os, tempfile
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri, default_converter
from rpy2.robjects.conversion import localconverter
import rpy2.rinterface_lib.callbacks as cb
import pandas as pd, numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

cb.consolewrite_print     = lambda x: None
cb.consolewrite_warnerror = lambda x: print(x, end="")

# ── setup (library + cache) ───────────────────────────────────────────────────
ro.r("library(data.table)")          # must come before importr("miCCI")
micci = importr("miCCI")
qm    = micci.load_quan_map()
cache = micci.precompute_lookups(micci.load_destatis(), qm)
ro.globalenv["qm_r"]  = qm
ro.globalenv["cch_r"] = cache
print("Ready.")

# ── helpers ───────────────────────────────────────────────────────────────────
def rval(v):
    x = v[0]; return x.item() if hasattr(x, "item") else x

def to_df(r):
    with localconverter(default_converter + pandas2ri.converter):
        return ro.conversion.rpy2py(r)

def run_r(code):
    """source() keeps data.table objects in .GlobalEnv across Python calls."""
    f = tempfile.NamedTemporaryFile(mode="w", suffix=".R", delete=False)
    f.write(code); f.close()
    try:    ro.r(f'source("{f.name}")')
    finally: os.unlink(f.name)

# ── single patient (the README example — now working) ─────────────────────────
result = micci.cci_gold(ro.StrVector(["E11.40", "N18.4", "I21.0"]), qm)
print("Gold CCI:", rval(result.rx2("cci")),
      " | active groups:", list(result.rx2("active").names))
# Gold CCI: 5  | active groups: ['kidney', 'mi', 'dm_complicated']

# ── batch S1–S6 on built-in 39-patient test dataset ───────────────────────────
run_r("""
data(micci_patients, package="miCCI")
dt <- as.data.table(micci_patients)
dt[, diagnosen:=list_full_icd_code]; dt[, stay_in_days:=los_in_days]; dt[, age:=60]

qm2 <- load_quan_map()
dt[, gold := cci_gold_batch(dt, qm2, build_pattern_lookup(qm2), build_dep_lookup(qm2))]

dt_anon <- copy(dt)
dt_anon[, diagnosen := sapply(list_full_icd_code, function(x) {
    paste(unique(substr(toupper(gsub("[^A-Z0-9]","",unlist(strsplit(x,"[|]+")))),1,3)), collapse="|")
})]

s1b <- cci_interval_batch(dt_anon, qm_r, cch_r)
dt[, S1 := s1b$cci_mid]; dt[, S1_min := s1b$cci_min]; dt[, S1_max := s1b$cci_max]
dt[, S2 := cci_probabilistic_batch(dt_anon, qm_r, cch_r)]
dt[, S3 := cci_mi_batch(dt_anon, qm_r, cch_r, m=5L)]
dt[, S4 := cci_bayesian_batch(dt_anon, qm_r, cch_r, n_draws=5L)]

if (requireNamespace("xgboost", quietly=TRUE)) {
    ml <- cci_ml_train(dt, qm_r)
    dt[, S5 := cci_ml_predict(ml, dt_anon, qm_r)]
} else {
    message("xgboost not found – S5 falls back to S3. Run: install.packages('xgboost')")
    dt[, S5 := S3]
}

mat <- as.matrix(dt[, .(S1, S2, S3, S4, S5)])
ens <- cci_ensemble_train(mat, dt$gold)
ens_b <<- ens
dt[, S6 := cci_ensemble(mat, ens$weights)]
""")

df = to_df(ro.r("dt[,.(p_id,gold,S1,S1_min,S1_max,S2,S3,S4,S5,S6)]")).round(2)
print(df[["p_id","gold","S1","S2","S3","S4","S5","S6"]].head(15).to_string(index=False))

# ── coverage metrics ──────────────────────────────────────────────────────────
strategies = ["S1","S2","S3","S4","S5","S6"]
gold = df["gold"].values
metrics = {}
for s in strategies:
    p = df[s].values
    metrics[s] = {
        "exact":   (np.abs(gold - p) < 0.5).mean() * 100,
        "within1": (np.abs(gold - p) < 1.5).mean() * 100,
        "within2": (np.abs(gold - p) < 2.5).mean() * 100,
        "mae":     np.abs(gold - p).mean(),
    }
s1_interval_cov = ((df["gold"] >= df["S1_min"]) & (df["gold"] <= df["S1_max"])).mean() * 100

# ── plots ─────────────────────────────────────────────────────────────────────
C = {"S1":"#4C72B0","S2":"#DD8452","S3":"#55A868","S4":"#C44E52","S5":"#8172B2","S6":"#937860"}
BG, DARK = "#f0f2f5", "#222222"
C_EXACT, C_W1, C_W2 = "#2ecc71", "#f39c12", "#e74c3c"

fig = plt.figure(figsize=(16, 13))
fig.patch.set_facecolor("#f8f9fa")
fig.suptitle("miCCI  ·  Strategy Coverage vs Gold CCI  (n = 39 patients)",
             fontsize=14, fontweight="bold", color=DARK, y=0.99)
gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.48, wspace=0.34, top=0.94, bottom=0.06)

# top: stacked coverage bars
ax1 = fig.add_subplot(gs[0, :])
x, bw = np.arange(len(strategies)), 0.55
exact_v = [metrics[s]["exact"]                         for s in strategies]
w1_v    = [metrics[s]["within1"] - metrics[s]["exact"]  for s in strategies]
w2_v    = [metrics[s]["within2"] - metrics[s]["within1"] for s in strategies]
b1      = [e + w for e, w in zip(exact_v, w1_v)]

ax1.bar(x, exact_v, bw, color=C_EXACT, label="Exact (±0)",  zorder=3)
ax1.bar(x, w1_v,    bw, bottom=exact_v, color=C_W1, label="Within ±1", zorder=3)
ax1.bar(x, w2_v,    bw, bottom=b1,      color=C_W2, label="Within ±2", zorder=3)
for xi, v in zip(x, exact_v):
    if v > 5:
        ax1.text(xi, v/2, f"{v:.0f}%", ha="center", va="center",
                 fontsize=10, fontweight="bold", color="white")
ax1.axhline(s1_interval_cov, color=C["S1"], lw=1.8, ls="--", zorder=4)
ax1.text(len(strategies)-0.28, s1_interval_cov+1.5,
         f"S1 interval coverage {s1_interval_cov:.0f}%",
         color=C["S1"], fontsize=8.5, ha="right")
ax1.set_xticks(x); ax1.set_xticklabels(strategies, fontsize=12)
ax1.set_ylabel("Patients (%)"); ax1.set_ylim(0, 118)
ax1.set_title("Coverage: % patients within tolerance of Gold CCI", fontsize=11, pad=6)
ax1.legend(loc="upper left", fontsize=9, framealpha=0.9, ncol=3)
ax1.set_facecolor(BG); ax1.grid(axis="y", color="white", lw=1.3, zorder=2)
ax1.spines[["top","right","left"]].set_visible(False)

# bottom: scatter per strategy
lim = max(gold.max(), max(df[s].max() for s in strategies)) + 1.5
for i, s in enumerate(strategies):
    r, c = [(1,0),(1,1),(1,2),(2,0),(2,1),(2,2)][i]
    ax = fig.add_subplot(gs[r, c])
    pred = df[s].values
    em = np.abs(gold-pred) < 0.5
    w1m = (np.abs(gold-pred) >= 0.5) & (np.abs(gold-pred) < 1.5)
    om  = np.abs(gold-pred) >= 1.5
    ax.scatter(gold[em],  pred[em],  color=C_EXACT, s=60, zorder=4, label="Exact")
    ax.scatter(gold[w1m], pred[w1m], color=C_W1,   s=60, zorder=4, label="±1")
    ax.scatter(gold[om],  pred[om],  color=C_W2,   s=60, zorder=4, label="±2+")
    if s == "S1":
        for _, row in df.iterrows():
            ax.plot([row["gold"]]*2, [row["S1_min"], row["S1_max"]],
                    color=C["S1"], alpha=0.4, lw=1.3, zorder=3)
    ax.plot([0,lim],[0,lim], color=DARK, lw=1, ls="--", alpha=0.4, zorder=3)
    ax.set_xlim(-0.5,lim); ax.set_ylim(-0.5,lim)
    ax.set_xlabel("Gold CCI", fontsize=8.5); ax.set_ylabel("Predicted CCI", fontsize=8.5)
    ax.set_title(f"{s}    exact {metrics[s]['exact']:.0f}%  MAE {metrics[s]['mae']:.2f}",
                 fontsize=10, fontweight="bold", color=C[s], pad=5)
    ax.set_facecolor(BG); ax.grid(color="white", lw=1.1, zorder=2)
    ax.spines[["top","right"]].set_visible(False)
    if i == 0:
        ax.legend(fontsize=8, loc="upper left", framealpha=0.85)

plt.savefig("micci_coverage.png", dpi=150, bbox_inches="tight")
plt.show()
print("Saved → micci_coverage.png")
```

A standalone Python example script is included at `scripts/micci_python_example.py`.

---

## Quick start

The fastest way to verify the installation and see all six strategies in action:

```r
library(miCCI)

# Run the demo on the built-in 39-patient dataset
# m=5 and n_draws=5 are reduced for speed; production defaults are 20 and 25
result <- run_demo(m = 5, n_draws = 5)

# Inspect the output table
print(result[, .(p_id, cci_gold, s1_mid, s2_ecci, s3_mi, s4_bayes)])
```

Expected output (truncated):

```
=== miCCI v0.4.3  --  Quick Demo
 Built-in dataset: micci_patients (39 patients)
...
--- Demo Results ---
  S1 Interval          MAE=0.xxx  RMSE=0.xxx  R²=0.xxx
  S2 Probabilistic     MAE=0.xxx  RMSE=0.xxx  R²=0.xxx
  S3 MI-CCI            MAE=0.xxx  RMSE=0.xxx  R²=0.xxx
  S4 Bayesian          MAE=0.xxx  RMSE=0.xxx  R²=0.xxx
  S5 ML-CCI            MAE=0.xxx  RMSE=0.xxx  R²=0.xxx
  S6 Ensemble          MAE=0.xxx  RMSE=0.xxx  R²=0.xxx
  S1 interval coverage: xx.x%
```

---

## Using each strategy independently

Every strategy is a self-contained function. You can call any one of them without running the others. The setup steps (loading the Quan map, downloading Destatis, building the cache) are shared and need to run only once per session.

### One-time setup

```r
library(miCCI)
library(data.table)

# Load the Quan et al. (2005) ICD-to-Charlson mapping
# This reads from inst/extdata/codes_quan_orig.json — no internet needed
quan_map <- load_quan_map()

# Download Destatis 23131-01 national frequency table
# Required for S2, S3, S4 (and the S6 ensemble)
# Internet connection needed the first time; result is cached in memory
destatis <- load_destatis()

# Build the prefix lookup cache — called once, reused by all strategies
cache <- precompute_lookups(destatis, quan_map)
```

### Gold-standard CCI (from full 4-digit codes)

```r
# Single patient — full 4-digit codes as a character vector
gold <- cci_gold(c("E11.40", "N18.4", "I21.0"), quan_map)
cat("CCI:", gold$cci)          # 5  (dm_compl=2 + kidney=2 + MI=1)
cat("Active groups:", names(gold$active))

# Batch of patients — data.table with a 'diagnosen' column (pipe-separated)
patients <- data.table(diagnosen = c(
  "E11.40|N18.4|I21.0",   # CCI = 5
  "J44.10|I50.0",          # CCI = 2
  "I10.00|E78.0"           # CCI = 0
))
pl   <- build_pattern_lookup(quan_map)
dl   <- build_dep_lookup(quan_map)
gold_vec <- cci_gold_batch(patients, quan_map, pl, dl)
print(gold_vec)  # 5 2 0
```

### S1 — Interval estimator (no external data needed)

```r
# Single patient — input is 3-character anonymized codes
s1 <- cci_interval(c("E11", "N18", "I21"), quan_map, cache)
cat("Min:", s1$cci_min, "  Max:", s1$cci_max,
    "  Mid:", s1$cci_mid, "  Width:", s1$interval_width)

# Batch
patients_anon <- data.table(diagnosen = c("E11|N18|I21", "J44|I50", "I10"))
s1_batch <- cci_interval_batch(patients_anon, quan_map, cache)
print(s1_batch)
#    cci_min cci_max cci_mid interval_width
# 1:       4       5     4.5              1
# 2:       2       2     2.0              0
# 3:       0       0     0.0              0
```

### S2 — Probabilistic expected value

```r
# Single patient
s2 <- cci_probabilistic(c("E11", "N18", "I21"), quan_map, cache)
cat("E-CCI:", s2$e_cci)

# Batch — returns a numeric vector
s2_vec <- cci_probabilistic_batch(patients_anon, quan_map, cache)
print(s2_vec)
```

### S3 — Multiple Imputation (best accuracy, slowest)

```r
# Single patient — m=20 for production, m=5 for quick testing
s3 <- cci_mi(c("E11", "N18", "I21"), quan_map, cache, m = 5)
cat("MI-CCI:", s3$mi_cci)

# Batch — m=20 takes ~4 hours for 500,000 patients
s3_vec <- cci_mi_batch(patients_anon, quan_map, cache, m = 5)
print(s3_vec)
```

### S4 — Bayesian Dirichlet estimator

```r
# Single patient
s4 <- cci_bayesian(c("E11", "N18", "I21"), quan_map, cache, n_draws = 5)
cat("Posterior median:", s4$posterior_median)

# Batch
s4_vec <- cci_bayesian_batch(patients_anon, quan_map, cache, n_draws = 5)
print(s4_vec)
```

### S5 — XGBoost machine learning

```r
# S5 requires a training set with gold-standard CCI labels
# Here we use the built-in dummy data as an illustration
data(micci_patients)
dt_train <- as.data.table(micci_patients)
dt_train[, diagnosen    := list_full_icd_code]
dt_train[, stay_in_days := los_in_days]
dt_train[, age          := 60]

# Train the model (also runs calibration if dt_calibrate is supplied)
ml_model <- cci_ml_train(dt_train, quan_map)

# Predict on anonymized codes
dt_anon <- copy(dt_train)
dt_anon[, diagnosen := sapply(list_full_icd_code, function(x)
  paste(unique(substr(gsub("[^A-Z0-9]","",toupper(
    unlist(strsplit(x,"\\|+")))),1,3)), collapse="|"))]
s5_vec <- cci_ml_predict(ml_model, dt_anon, quan_map)
print(s5_vec)
```

### S6 — NNLS Ensemble (requires S1-S5 outputs)

```r
# IMPORTANT: S6 takes the predictions of S1-S5 as its input,
# not raw ICD codes. You must run all five strategies first.

# Assemble the five predictions into a matrix (one column per strategy)
est_matrix <- cbind(
  S1_mid   = s1_batch$cci_mid,
  S2_ecci  = s2_vec,
  S3_mi    = s3_vec,
  S4_bayes = s4_vec,
  S5_ml    = s5_vec
)

# Train weights on a calibration set (here using the same data for illustration)
ens <- cci_ensemble_train(est_matrix, gold_vec)
cat("Weights:", sprintf("%s=%.3f", names(ens$weights), ens$weights), "\n")

# Apply to new patients
s6_vec <- cci_ensemble(est_matrix, ens$weights)
print(s6_vec)
```

---

## Using all strategies together

A complete workflow from setup through all six strategies with accuracy evaluation:

```r
library(miCCI)
library(data.table)

# ── 1. Load reference data ────────────────────────────────────────────────────
quan_map <- load_quan_map()
destatis  <- load_destatis()
cache     <- precompute_lookups(destatis, quan_map)
pl        <- build_pattern_lookup(quan_map)
dl        <- build_dep_lookup(quan_map)

# ── 2. Prepare patient data ───────────────────────────────────────────────────
data(micci_patients)
dt <- as.data.table(micci_patients)
dt[, diagnosen    := list_full_icd_code]
dt[, stay_in_days := los_in_days]
dt[, age          := 60]

# Compute gold CCI (from full codes — this is the "truth" for validation)
dt[, cci_gold := cci_gold_batch(dt, quan_map, pl, dl)]

# Simulate §21 anonymization: truncate to 3 characters
dt[, diagnosen_anon := sapply(list_full_icd_code, function(x) {
  codes <- unlist(strsplit(x, "\\|+"))
  paste(unique(substr(toupper(gsub("[^A-Z0-9]","",codes)),1,3)), collapse="|")
})]
dt_anon        <- copy(dt)
dt_anon[, diagnosen := diagnosen_anon]

# ── 3. Run all six strategies ─────────────────────────────────────────────────
s1  <- cci_interval_batch(dt_anon, quan_map, cache)
dt[, c("s1_min","s1_max","s1_mid") := .(s1$cci_min, s1$cci_max, s1$cci_mid)]
dt[, s2 := cci_probabilistic_batch(dt_anon, quan_map, cache)]
dt[, s3 := cci_mi_batch(dt_anon, quan_map, cache, m = 5L)]
dt[, s4 := cci_bayesian_batch(dt_anon, quan_map, cache, n_draws = 5L)]
ml   <- cci_ml_train(dt, quan_map)
dt[, s5 := cci_ml_predict(ml, dt_anon, quan_map)]

# S6: assemble S1-S5 predictions, train ensemble weights
mat  <- as.matrix(dt[, .(s1_mid, s2, s3, s4, s5)])
ens  <- cci_ensemble_train(mat, dt$cci_gold)
dt[, s6 := cci_ensemble(mat, ens$weights)]

# ── 4. Evaluate accuracy ──────────────────────────────────────────────────────
for (col in c("s1_mid","s2","s3","s4","s5","s6")) {
  m <- eval_metrics(dt[[col]], dt$cci_gold)
  cat(sprintf("%-8s  MAE=%.3f  R2=%.3f\n", col, m["mae"], m["r_squared"]))
}
```

---

## Built-in dummy dataset

```r
data(micci_patients)
str(micci_patients)
# 'data.frame': 39 obs. of 3 variables:
#  $ p_id               : chr  "PAT_001" "PAT_002" ...
#  $ list_full_icd_code : chr  "I10.00|E78.00|Z96.41" "J06.9|R05" ...
#  $ los_in_days        : num  2 1 3 1 2 5 7 6 4 8 ...
```

The dataset covers:

- **Rows 1-5**: Healthy patients (CCI = 0) — hypertension, respiratory infections, back pain
- **Rows 6-13**: Single comorbidity (CCI = 1) — MI, HF, COPD, simple DM, cerebrovascular disease
- **Rows 14-18**: Moderate burden (CCI = 2) — complicated DM, kidney disease, cancer, hemiplegia
- **Rows 19-20**: CCI = 3 — severe liver disease, triple comorbidity
- **Rows 21-30**: High burden (CCI = 4 to 17)
- **Rows 31-34**: `depends_on` hierarchy test cases — verify suppression logic
- **Rows 35-39**: Edge cases — rheumatic disease, dementia variants, heart failure subtypes

---

## Running the full validation pipeline

If you have access to the UMM Parquet extract (requires institutional data access):

```r
library(miCCI)

run_pipeline(
  data_path   = "A:/HLZ/Projekt/2025/comorbidity_gaetan/comorbidity_diagnoses_only_pseudonym.parquet",
  output_dir  = "results/",
  sample_size = NULL    # Set to e.g. 10000 for a quick test run
)
```

Output files written to `results/`:

| File | Description |
|---|---|
| `cci_predictions_validation.csv` | Per-encounter predictions for all 6 strategies |
| `cci_metrics_summary.csv` | MAE, RMSE, R² across all splits |
| `ensemble_weights.csv` | NNLS weights for S6 |
| `fig_01_mae_comparison.pdf` | Bar chart of MAE per strategy |
| `fig_02_predicted_vs_gold.pdf` | Hexbin scatter (predicted vs. gold) |
| `fig_03_temporal_stability.pdf` | MAE by year |
| `fig_04_density_overlay.pdf` | CCI distribution: gold vs. all strategies |
| `fig_05_bland_altman.pdf` | Bland-Altman agreement plots |
| `tab_01_descriptive_statistics.csv` | Summary statistics |
| `tab_02_agreement_metrics.csv` | Bias, LoA, coverage |

---

## Running tests

```r
# From within RStudio or an R console in the package directory
devtools::test()

# Alternatively, from the terminal
Rscript -e "devtools::test('path/to/miCCI')"
```

Tests print a label for every check, describing **what** is being tested and **what** the expected outcome is. For example:

```
[ FAIL 0 | WARN 0 | SKIP 0 | PASS 3 ]
── Test results ─────────────────────────────────────────────────────────
  Testing: normalize_icd('e11.40') -> 'E11.40'
  Testing: leading/trailing spaces are stripped
  Testing: comma is converted to decimal point
✓ [UTIL] normalize_icd: lowercase input is uppercased and trimmed
...
```

---

## Citation

If you use miCCI in published research, please cite:

> Kamdje Wabo G, Sokolowski P, Siegel F, Ganslandt T. Estimating the Charlson
> Comorbidity Index from anonymized diagnosis codes while preserving data quality:
> development and validation of six computational strategies. *JMIR Medical Informatics*.
> 2025 (in preparation).

### References

1. Charlson ME, Pompei P, Ales KL, MacKenzie CR. A new method of classifying prognostic
   comorbidity in longitudinal studies: development and validation. *J Chronic Dis*.
   1987;40(5):373-383. doi:10.1016/0021-9681(87)90171-8

2. Quan H, Sundararajan V, Halfon P, et al. Coding algorithms for defining comorbidities
   in ICD-9-CM and ICD-10 administrative data. *Med Care*. 2005;43(11):1130-1139.
   doi:10.1097/01.mlr.0000182534.19832.83

3. §21 Krankenhausentgeltgesetz (KHEntgG). Übermittlungspflicht der Krankenhäuser.
   Bundesministerium der Justiz. 2004.

4. Statistisches Bundesamt (Destatis). Diagnosedaten der Krankenhäuser
   (Krankenhausdiagnosestatistik) 23131-01. Wiesbaden: Destatis; 2024.
   https://www.destatis.de/static/DE/dokumente/5231301237015_SB.xlsx

---

## License

MIT © 2026 Gaetan Kamdje Wabo et al.
See [LICENSE](LICENSE) for full terms.
