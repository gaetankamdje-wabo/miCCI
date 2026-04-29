# miCCI — Charlson Comorbidity Index from Anonymized ICD-10-GM Codes

`miCCI` reconstructs the Charlson Comorbidity Index (CCI) from truncated
(anonymised) three-character ICD-10-GM codes using five complementary
strategies and a cross-validated meta learner that combines them.

| Strategy | Idea | Input |
|---|---|---|
| **S1** Interval | Lower / upper / midpoint of the CCI given what the truncated codes alone imply | Frequency table (only as taxonomy) |
| **S2** Probabilistic | Inclusion–exclusion expected CCI: `P(group active) = 1 − Π(1 − qᵢ)` | Frequency table |
| **S3** Multiple Imputation | Draw full-length codes from the empirical subcode distribution; recompute gold; mean across `m` draws | Frequency table |
| **S4** Bayesian | Dirichlet posterior over subcode probabilities; posterior median across `n_draws` | Frequency table |
| **Meta** | Cross-validated NNLS Super Learner over S1..S4 | Gold and the four strategy outputs |

## Installation

```r
# install.packages("remotes")
remotes::install_local("miCCI_1.0.0.tar.gz")
```

## Quick start (Destatis as the default reference)

```r
library(miCCI)

# 1) Charlson mapping (Quan et al. 2005, enhanced ICD-10).
quan_map <- load_quan_map()

# 2) Population frequencies. Defaults to Destatis 23131-01.
freq <- load_destatis()

# 3) Precompute per-prefix lookups (one-time, ~30s).
cache <- precompute_lookups(freq, quan_map)

# 4) One encounter, all strategies.
codes_anon <- c("E11", "N18", "I10")
cci_interval(codes_anon, quan_map, cache)
cci_probabilistic(codes_anon, quan_map, cache)
cci_mi(codes_anon, quan_map, cache, m = 20L)
cci_bayesian(codes_anon, quan_map, cache, n_draws = 25L, alpha_0 = 10)
```

## Using your own population frequencies

Anywhere you would call `load_destatis()`, pass your own `data.table`
with the same schema instead. Required columns:

| Column        | Type      | Meaning |
|---------------|-----------|---------|
| `code`        | character | Dotted four-character ICD-10-GM code (e.g. `"E11.4"`). |
| `code_nodot`  | character | The same code without the dot (e.g. `"E114"`). |
| `code3`       | character | The three-character prefix (e.g. `"E11"`). |
| `freq_total`  | numeric   | Encounter count for that code in your reference population. |
| *age bands*   | numeric   | Optional. One column per Destatis age bin (`"unter 1"`, `"1 - 5"`, …, `"95 u. älter"`). If supplied, miCCI can produce age-stratified subcode probabilities. |

Then:

```r
# 'my_freq' is a data.table with the schema above.
cache  <- precompute_lookups(my_freq, quan_map)

# Pipeline run with your own reference instead of Destatis:
compute_predictions(
  data_path  = "path/to/cohort.parquet",
  output_dir = "out/",
  freq_table = my_freq,                 # <-- your reference goes here
  mi_m        = 20L,
  bayes_draws = 25L,
  alpha_0     = 10,
  sl_cv_folds = 10L,
  seed        = 42L
)
```

If `freq_table` is `NULL`, miCCI falls back to `load_destatis()`.

### Minimal example: building `freq_table` from scratch

```r
library(data.table)
my_freq <- data.table(
  code        = c("E11.4", "N18.4", "I10.0"),
  code_nodot  = c("E114",  "N184",  "I100"),
  code3       = c("E11",   "N18",   "I10"),
  freq_total  = c(123456,   45678,  789012)
)
cache <- precompute_lookups(my_freq, quan_map)
```

Age-stratified columns are optional — if you do not have them, every age
will fall back to the marginal distribution.

## Two-stage pipeline on a cohort

```r
# Stage 1: heavy compute (S1..S4 + Meta on the full cohort). Run once.
compute_predictions(
  data_path  = "cohort.parquet",
  output_dir = "out/"
)

# Stage 2: tables + figures from the persisted artefact. Run as often
# as you like for cosmetic / reviewer-driven tweaks.
build_report(input_dir = "out/")
```

## Reproducibility

* Every internal RNG call is wrapped so the user's `.Random.seed` is
  restored after the call. Running miCCI inside a larger simulation
  never silently resets your RNG stream.
* `seed` is a parameter on every public function that uses randomness.

## Citation

```r
citation("miCCI")
```

## License

MIT.
