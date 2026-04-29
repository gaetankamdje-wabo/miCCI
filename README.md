# miCCI v0.5.0

Quality-preserving Charlson Comorbidity Index (CCI) estimation under
ICD-code anonymisation.

## What is new in v0.5.0

* **Single consistent cohort 2010-2024.** No train / calibrate / validate
  split. All strategies are training-free and can be applied to the entire
  cohort directly.
* **S5 (XGBoost ML) removed.** Conceptually outside the paper.
* **S6 reworked into a deterministic meta-estimator** that combines S1-S4
  on every encounter without any learned parameters
  (`R/60_cci_meta.R`, formula in the manuscript).
* **All probabilistic outputs are floats.** Rounding is left to the user.
* **Two-stage pipeline.** A long compute stage that produces a self-
  contained predictions artefact, and a fast reporting stage that
  regenerates every table and plot from that artefact in minutes.
* **Bootstrap 95% confidence intervals** for MAE / RMSE / R^2 / Bias on
  every strategy and on every stratum.
* **Quality-assurance module** comparing the JSON Quan map, the gold
  standard, and the predictions on coverage / mass conservation /
  distributional fidelity (Plot 6, Tables 3-6).
* **Stratified diagnostic-cluster benchmark** identifying the best
  strategy per ICD chapter, per Charlson group and per multimorbid
  pattern, with a categorical recommendation per stratum.
* **Plot 02 (Predicted vs Gold) rewritten** at base font 16, axis text 14,
  bold facet labels, axes 0-20, per-facet MAE+CI annotation, both PDF
  (cairo vector) and PNG (300 dpi).

## Two-stage workflow

```
+----------------------+         +------------------------+
| 01_run_compute.R     |  -->    | predictions.parquet    |
| (S1-S4 + S6, ~16h)   |         | run_manifest.json      |
+----------------------+         +------------------------+
                                              |
                                              v
                                 +------------------------+
                                 | 02_run_report.R        |
                                 | (CIs, QA, plots, ~min) |
                                 +------------------------+
                                              |
                                              v
                                 tab_*.csv  +  fig_*.{pdf,png}
```

## Run instructions

### 1. Install (one time)

```r
# from the unpacked miCCI/ directory
devtools::install(".")
# or, without installing:
# devtools::load_all(".")
```

Required packages: `data.table`, `jsonlite`, `readxl`, `arrow`, `ggplot2`.
Optional: `future.apply` (parallel bootstrap), `scales`.

### 2. Stage 1 — compute predictions (~16 h on the full Mannheim cohort)

Edit `scripts/01_run_compute.R`:

```r
DATA_PATH   <- "<path to comorbidity_diagnoses_only_pseudonym.parquet>"
OUTPUT_DIR  <- "<directory for results>"
SAMPLE_SIZE <- NULL      # set 5000L for a smoke test
```

Then run:

```bash
Rscript scripts/01_run_compute.R
```

Output of this stage:

```
<OUTPUT_DIR>/predictions.parquet
<OUTPUT_DIR>/run_manifest.json
```

After it has finished, **never run it again unless the cohort changes**.

### 3. Stage 2 — build the report (~minutes)

Edit `scripts/02_run_report.R` if you want a different `INPUT_DIR`,
bootstrap budget or pattern count, then:

```bash
Rscript scripts/02_run_report.R
```

Output of this stage in the same directory:

```
tab_01_strategy_metrics_ci.csv          Global per-strategy MAE/RMSE/R^2/Bias + 95% CI
tab_02_descriptive_statistics.csv       Gold vs each strategy
tab_03_qa_group_coverage.csv            Per Charlson group: gold-active rate + mass ratios
tab_04_qa_score_frequency.csv           Frequency of CCI bins per source
tab_05_qa_ks_distance.csv               KS distance of each strategy vs gold
tab_06_qa_posterior_coverage.csv        % within +/-0.5 / 1 / 2 CCI points
tab_07_stratified_chapter.csv           Per ICD chapter
tab_08_stratified_charlson_group.csv    Per gold-active Charlson group
tab_09_stratified_pattern.csv           Per multimorbid pattern (top 30)
tab_10_winners_chapter.csv              Best strategy per chapter (with verdict)
tab_11_winners_group.csv                Best strategy per Charlson group
tab_12_winners_pattern.csv              Best strategy per pattern

fig_01_mae_with_ci.{pdf,png}            MAE comparison with bootstrap CIs
fig_02_predicted_vs_gold.{pdf,png}      HIGH DPI, big fonts, per-facet captions
fig_03_temporal_stability.{pdf,png}     MAE by year
fig_04_density_overlay.{pdf,png}        Distributional fidelity
fig_05_bland_altman.{pdf,png}           Per-strategy agreement
fig_06_qa_mass_conservation.{pdf,png}   QA: per-group mass ratios
fig_07_stratified_chapter.{pdf,png}     Stratified MAE by chapter
fig_08_stratified_charlson_group.{pdf,png}  Stratified MAE by Charlson group
fig_09_stratified_pattern.{pdf,png}     Stratified MAE by multimorbid pattern
```

### 4. Reviewer-driven cosmetic edits

If a reviewer asks for a different colour, size, font or label on any plot,
edit the corresponding block in `R/95_report.R` and **only re-run Stage 2**:

```bash
Rscript scripts/02_run_report.R
```

The 16-hour compute stage stays cached in `predictions.parquet`.

## S6 deterministic meta-estimator

```
m_dist  = mean(s2_ecci, s3_mi, s4_bayes)
w       = s1_width / (s1_width + kappa)        with kappa = 2
raw     = (1 - w) * s1_mid + w * m_dist
s6_meta = clamp(raw, s1_min, s1_max)
```

* No training, no calibration, no learned parameters.
* Bounded by S1's certainty interval -> inherits S1's coverage.
* Reduces to `s1_mid` when the interval collapses.
* Reduces to the consensus of (S2, S3, S4) when the interval is wide.
* `kappa = 2` is paper-stated as "twice the smallest non-trivial Charlson
  weight" and exposed as `MICCI_META_KAPPA`.

## Output dtypes

| Column     | Type  |
|------------|-------|
| `cci_gold` | int   |
| `s1_min`   | float |
| `s1_max`   | float |
| `s1_mid`   | float |
| `s1_width` | float |
| `s2_ecci`  | float |
| `s3_mi`    | float |
| `s4_bayes` | float |
| `s6_meta`  | float |

Users can apply `as.integer(round(.))` downstream if they need integer scores.
