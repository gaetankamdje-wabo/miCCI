# miCCI — Charlson Comorbidity Index from Anonymised ICD-10-GM Codes

`miCCI` reconstructs the Charlson Comorbidity Index (CCI) from truncated
(anonymised) three-character ICD-10-GM codes using five complementary
strategies plus a cross-validated meta learner.

| Strategy | Idea | Key parameters |
|---|---|---|
| **S1** Interval | Lower / upper / midpoint of the CCI given what the truncated codes alone imply | none |
| **S2** Probabilistic | Inclusion–exclusion expected CCI: `P(group active) = 1 − Π(1 − qᵢ)` | none |
| **S3** Multiple Imputation | Draw full-length codes from the empirical subcode distribution; recompute gold; mean across `m` draws | `m`, `seed` |
| **S4** Bayesian | Dirichlet posterior over subcode probabilities; posterior median across `n_draws` | `n_draws`, `alpha_0`, `seed` |
| **Meta** | Cross-validated NNLS Super Learner over S1..S4 | `V` (folds), `seed` |

This package provides the **estimators only**. Bootstrap, stratified
evaluation, plotting, and any pipeline orchestration are deliberately
left to the user — `miCCI` produces the per-encounter predictions; how
you analyse them is up to you.

## Installation

```r
# install.packages("remotes")
remotes::install_github("gaetankamdje-wabo/miCCI")

```

`SuperLearner` is needed only for the meta learner (`cci_meta_fit`); the
five base strategies work without it. Install on demand:

```r
install.packages("SuperLearner")
```

## Quick start

```r
library(miCCI)
library(data.table)

# 1) Charlson mapping (Quan et al. 2005, enhanced ICD-10).
quan_map <- load_quan_map()

# 2) Reference frequencies. Defaults to Destatis 23131-01 (Germany).
freq <- load_destatis()

# 3) Precompute per-prefix lookups (one-time, ~30s).
cache <- precompute_lookups(freq, quan_map)

# 4) One encounter, all strategies.
codes_anon <- c("E11", "N18", "I10")
cci_interval(codes_anon, quan_map, cache)
cci_probabilistic(codes_anon, quan_map, cache)
cci_mi(codes_anon, quan_map, cache, m = 20L, seed = 42L)
cci_bayesian(codes_anon, quan_map, cache, n_draws = 25L, alpha_0 = 10, seed = 42L)
```

## Cohort-scale usage

The vectorised batch functions take a `data.table` with a column
`diagnosen` containing pipe-separated ICD-10-GM codes per encounter and
return one row of estimates per encounter:

```r
dt <- data.table(diagnosen = c(
  "E11.4|N18.4|I10.0",
  "K70.4|K74.6",
  "I21.0|I50.0|J44.0"
))

dt[, cci_gold := cci_gold_batch(dt, quan_map)]

s1 <- cci_interval_batch(dt, quan_map, cache)
dt[, c("s1_min","s1_max","s1_mid","s1_width") := .(
  s1$cci_min, s1$cci_max, s1$cci_mid, s1$interval_width
)]

dt[, s2_ecci := cci_probabilistic_batch(dt, quan_map, cache)]
dt[, s3_mi   := cci_mi_batch(dt, quan_map, cache, m = 20L, seed = 42L)]
dt[, s4_bayes := cci_bayesian_batch(dt, quan_map, cache,
                                    n_draws = 25L, alpha_0 = 10, seed = 42L)]

# Cross-validated meta learner over S1..S4
meta <- cci_meta_fit(dt[, .(cci_gold, s1_min, s1_max, s1_mid,
                            s2_ecci, s3_mi, s4_bayes)],
                     V = 10L, seed = 42L)
dt[, meta := meta$predictions]
meta$weights   # NNLS weights over S1..S4
```

By construction, `s1_min <= cci_gold <= s1_max` for any encounter whose
original full-length codes are in the frequency table that built `cache`.

## Using your own population frequencies

`load_destatis()` is just one possible reference. Pass any `data.table`
of the same schema to `precompute_lookups()` instead:

| Column        | Type      | Meaning |
|---------------|-----------|---------|
| `code`        | character | Dotted four-character ICD-10-GM code (e.g. `"E11.4"`). |
| `code_nodot`  | character | The same code without the dot (`"E114"`). |
| `code3`       | character | The three-character prefix (`"E11"`). |
| `freq_total`  | numeric   | Encounter count for that code in your reference population. |

```r
my_freq <- data.table(
  code        = c("E11.4", "N18.4", "I10.0"),
  code_nodot  = c("E114",  "N184",  "I100"),
  code3       = c("E11",   "N18",   "I10"),
  freq_total  = c(123456,   45678,  789012)
)
cache_local <- precompute_lookups(my_freq, quan_map)
```

## Reproducibility

Every randomised call (`cci_mi*`, `cci_bayesian*`, `cci_meta_fit`)
restores the caller's `.Random.seed` after running, so calling these
functions inside a larger simulation never silently mutates your RNG
stream. The `seed` argument controls the inner stream.

## Citation

```r
citation("miCCI")
```

## License

MIT.
