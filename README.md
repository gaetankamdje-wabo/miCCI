# miCCI

**Estimate the Charlson Comorbidity Index from anonymized ICD-10-GM codes.**

When hospital data integration centers truncate ICD-10-GM codes to three characters for privacy, standard CCI algorithms break. miCCI reconstructs the CCI by using national subcode frequencies from the German Federal Statistical Office (Destatis) as a population prior.

Five strategies, one package. From a fast deterministic bound to a full Bayesian model.

| Strategy | What it does | Speed |
|----------|-------------|-------|
| **S1** Interval | Reports the [min, max] range; point estimate = midpoint | instant |
| **S2** Probabilistic | Weights each Charlson group by its Destatis probability | instant |
| **S3** MI-CCI | Draws subcodes 20x from Destatis, averages the CCI | seconds |
| **S4** Bayesian | Draws subcodes 25x from a Dirichlet prior, takes the median | seconds |
| **S5** Meta Learner | Combines S1 to S4 via cross-validated Super Learner (NNLS) | minutes |

---

## Install

```r
# Option 1: from GitHub
devtools::install_github("gaetankamdje-wabo/miCCI")

# Option 2: from a local clone
git clone https://github.com/gaetankamdje-wabo/miCCI.git
devtools::install("miCCI")
```

### Required R packages (installed automatically)

```
data.table, jsonlite, readxl, arrow, ggplot2, SuperLearner
```

### Optional

```
future.apply   (parallel bootstrap)
scales         (plot label formatting)
testthat       (unit tests)
```

---

## Your data must look like this

One row per encounter. One column `diagnosen` with ICD-10-GM codes separated by `|`.

```
falnr   | age | diagnosen
--------|-----|------------------------------
P00001  | 72  | I25.1|E11.9|N18.3|J44.1
P00002  | 45  | K80.1|K85.9
P00003  | 68  | C34.1|C78.0|E11.4|I10.0
```

Codes can be full length (`E11.9`) or already truncated to three characters (`E11`). miCCI handles both.

---

## Quick start: simulate and run

This self-contained example generates synthetic encounters and runs all five strategies. Copy-paste it into R.

```r
library(miCCI)
library(data.table)

# ── 1. Create toy encounters ──────────────────────────────────────────
toy <- data.table(
  falnr    = paste0("P", 1:6),
  age      = c(72, 45, 68, 55, 80, 33),
  diagnosen = c(
    "I25.1|E11.9|N18.3|J44.1",
    "K80.1|K85.9",
    "C34.1|C78.0|E11.4|I10.0",
    "E11.9|I10.0",
    "C61|C78.5|G30.1|F00.1",
    "S72.0|M80.0"
  )
)

# ── 2. Load the Quan mapping and Destatis frequencies ─────────────────
quan_map <- load_quan_map()
dt_dest  <- load_destatis()

# ── 3. Build the lookup cache (run ONCE, reuse for all strategies) ────
cache <- precompute_lookups(dt_dest, quan_map)

# ── 4. Gold standard (requires full-length codes) ────────────────────
gold <- cci_gold(c("I25.1", "E11.9", "N18.3", "J44.1"), quan_map)
gold$cci      # 4
gold$active   # shows which Charlson groups fired
```

### Run each strategy on a single encounter

```r
codes <- c("I25.1", "E11.9", "N18.3", "J44.1")

# S1: Deterministic interval
cci_interval(codes, quan_map, cache)
# => list(cci_min, cci_max, cci_mid, interval_width)

# S2: Probabilistic heuristic
cci_probabilistic(codes, quan_map, cache)
# => list(e_cci)

# S3: Multiple imputation (20 draws, mean)
cci_mi(codes, quan_map, cache, m = 20L, seed = 42L)
# => list(mi_cci)

# S4: Bayesian Dirichlet prior (25 draws, median)
cci_bayesian(codes, quan_map, cache, n_draws = 25L, seed = 42L)
# => list(posterior_median)
```

### Run all strategies at once on a full cohort

```r
# compute_predictions() runs S1 to S5 in sequence and writes
# a self-contained artefact to disk.
#
# This is designed for the real Mannheim cohort (~16h on 500k encounters).
# For the toy example:

result <- compute_predictions(
  data_path   = "path/to/your/data.parquet",
  output_dir  = "results/",
  sample_size = NULL,     # set 5000L for a smoke test
  mi_m        = 20L,      # S3: imputation rounds
  bayes_draws = 25L,      # S4: Dirichlet draws
  sl_cv_folds = 10L,      # S5: Super Learner CV folds
  seed        = 42L
)

# Output: results/predictions.parquet + results/run_manifest.json
```

### Generate all tables and figures from a completed run

```r
build_report(
  input_dir  = "results/",
  output_dir = "results/",
  B          = 1000L,     # bootstrap resamples (global)
  B_strata   = 500L       # bootstrap resamples (per stratum)
)

# Output: tab_01..tab_12 CSV files + fig_01..fig_09 PDF/PNG
```

---

## Function reference

### Data loading

| Function | Purpose |
|----------|---------|
| `load_quan_map()` | Load the Quan et al. ICD-10 to Charlson group mapping (bundled JSON) |
| `load_destatis()` | Download and parse Destatis 23131-01 subcode frequencies |
| `precompute_lookups(dt_destatis, quan_map)` | Build the internal cache. Run once, pass to all strategies |
| `load_cohort(path)` | Load a parquet file in the expected schema (falnr, age, diagnosen, ...) |

### Gold standard

| Function | Purpose |
|----------|---------|
| `cci_gold(icd_codes, quan_map)` | CCI from full-length codes (single encounter) |

### Strategy functions (single encounter)

| Function | Strategy | Returns |
|----------|----------|---------|
| `cci_interval(codes, quan_map, cache)` | S1 | `list(cci_min, cci_max, cci_mid, interval_width)` |
| `cci_probabilistic(codes, quan_map, cache)` | S2 | `list(e_cci)` |
| `cci_mi(codes, quan_map, cache, m, seed)` | S3 | `list(mi_cci)` |
| `cci_bayesian(codes, quan_map, cache, n_draws, seed)` | S4 | `list(posterior_median)` |

### Pipeline functions (full cohort)

| Function | Purpose |
|----------|---------|
| `compute_predictions(data_path, output_dir, ...)` | Stage 1: run S1 to S5 on every encounter, write predictions.parquet |
| `load_predictions(output_dir)` | Read a previously computed predictions artefact |
| `build_report(input_dir, output_dir, ...)` | Stage 2: generate all tables and figures from predictions.parquet |
| `cci_meta_fit(dt, V, seed)` | Fit the S5 Super Learner (called internally by compute_predictions) |

### Evaluation and QA

| Function | Purpose |
|----------|---------|
| `bootstrap_metrics(gold, predicted, B)` | MAE, RMSE, R², Bias with bootstrap 95% CI |
| `bootstrap_strategies(preds, strategies, B)` | Run bootstrap_metrics for all strategies at once |
| `qa_group_coverage(preds, quan_map)` | Per-Charlson-group mass conservation check |
| `qa_score_distribution(preds, strategy_cols)` | CCI bin frequencies + KS distance vs gold |
| `qa_posterior_coverage(preds, strategy_cols)` | % of encounters within ±0.5 / ±1.0 / ±2.0 of gold |
| `stratified_evaluation(preds, strategies, ...)` | MAE by ICD chapter, Charlson group, multimorbid pattern |
| `stratum_winners(strat_dt)` | Pick the best strategy per stratum with decisive/tied verdict |

---

## Two-stage pipeline

```
Stage 1: compute_predictions()          Stage 2: build_report()
(S1 to S5, ~16h on 500k encounters)    (tables + figures, ~minutes)

   your_data.parquet                     predictions.parquet
         |                                      |
         v                                      v
  predictions.parquet  ──────────────>   tab_01..12.csv
  run_manifest.json                      fig_01..09.pdf/png
```

**Stage 1 runs once.** Stage 2 runs as often as needed for cosmetic edits.

---

## Output columns (predictions.parquet)

| Column | Type | Description |
|--------|------|-------------|
| `cci_gold` | int | Gold-standard CCI from full-length codes |
| `s1_min` | int | S1 lower bound |
| `s1_max` | int | S1 upper bound |
| `s1_mid` | float | S1 midpoint estimate |
| `s2_ecci` | float | S2 expected CCI |
| `s3_mi` | float | S3 multiple imputation mean |
| `s4_bayes` | float | S4 Bayesian prior median |
| `s6_meta` | float | S5 Meta Learner (column name kept for backward compatibility) |

To convert any float to integer: `as.integer(round(x))`

---

## S5 Meta Learner

The Meta Learner is a cross-validated Super Learner (van der Laan, Polley & Hubbard 2007). It takes the outputs of S1 to S4, learns non-negative weights via NNLS over 10-fold CV, and returns their convex combination. It does not look at ICD codes or Destatis frequencies.

The internal column is named `s6_meta` for backward compatibility with earlier development versions. All display labels and figures use "Meta".

---

## Citation

If you use miCCI, please cite:

> Kamdje Wabo G, Santhanam N, Jannesari Ladani M, Hagmann M, Sokolowski PP,
> Ganslandt T, Siegel F. Estimating the Charlson Comorbidity Index From
> Anonymized Diagnosis Codes While Preserving Data Quality: Methodological
> Study Comparing Five Computational Strategies. *JMIR* (submitted).

---

## License

MIT. See [LICENSE](LICENSE).
