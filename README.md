# miCCI

### Quality-Preserving Charlson Comorbidity Index Estimation Under Diagnosis Truncation

`miCCI` answers one question: *if a patient's diagnosis codes have been cut down to three characters for privacy reasons, how sick were they really?*

It reconstructs the Charlson Comorbidity Index (CCI) from truncated ICD-10-GM codes using five complementary strategies, plus a data-driven meta learner that blends them. Everything below was run and checked against the real package source in this repository, including a live call to the German Federal Statistical Office (Destatis) reference table, so every number you see is real, not illustrative.

---

## Contents

1. [Rationale](#rationale)
2. [The Charlson Comorbidity Index](#the-charlson-comorbidity-index)
3. [Conceptual Overview](#conceptual-overview)
4. [Estimation Strategies](#estimation-strategies)
5. [Configuration Parameters](#configuration-parameters)
6. [Requirements](#requirements)
7. [Installation](#installation)
8. [Loading](#loading)
9. [Testing](#testing)
10. [Quick Start](#quick-start)
11. [Full-Length Codes](#full-length-codes)
12. [Use Case A: Destatis, No Age](#use-case-a-destatis-no-age)
13. [Use Case B: Destatis, With Age Stratification](#use-case-b-destatis-with-age-stratification)
14. [Use Case C: Custom Reference Population, No Age Bands](#use-case-c-custom-reference-population-no-age-bands)
15. [Use Case D: Custom Reference Population, With Age Bands](#use-case-d-custom-reference-population-with-age-bands)
16. [Function Reference](#function-reference)
17. [Reproducibility](#reproducibility)
18. [Implementation Notes](#implementation-notes)
19. [Citation](#citation)
20. [License](#license)

---

## Rationale

Hospitals, insurers, and national statistics offices routinely share diagnosis data for research, but they do not share it whole. A common privacy safeguard is to cut every ICD-10-GM code down to its three-character category. `E11.42` (type 2 diabetes with severe kidney complications) becomes `E11`. `N18.4` (stage 4 chronic kidney disease) becomes `N18`. The fourth and fifth characters, the ones that actually carry the clinical severity, are gone.

That is a problem for the Charlson Comorbidity Index, because the CCI does not treat every diabetes code the same way. `E11.9` (uncomplicated diabetes) contributes one point. `E11.4` (diabetes with neurological complications) contributes two, and it silently cancels the one point from the uncomplicated version if both are somehow present, because the more severe condition subsumes the milder one. Once the fourth character is gone, that distinction is gone too, and a naive researcher is left guessing whether a truncated `E11` was mild or severe.

`miCCI` does not guess blindly. It uses the fact that most of the codes that used to hide behind a given three-character prefix, before truncation, are not equally likely. Type 2 diabetes without complications is far more common than diabetes with severe end-organ damage, and that imbalance shifts with age, because comorbidity burden rises as people get older. National statistics on how often each four-character code actually occurs, stratified by age, give miCCI enough information to reconstruct a defensible CCI estimate from the truncated fragment alone, without ever seeing the original patient-level code.

The package offers five ways to do this reconstruction, from a simple deterministic bound to a fully Bayesian posterior, plus a meta learner that works out which of the five to trust in which situation. The whole approach, the five strategies, and their behaviour on 507,949 real inpatient encounters at Universitätsmedizin Mannheim are described in the accompanying methodological study (see [Citation](#citation)). The ICD-10-GM to Charlson mapping bundled with the package, the file `inst/extdata/codes_quan_orig.json`, follows the German-specific adaptation of the Quan et al. (2005) enhanced ICD-10 algorithm published in Sokołowski et al. (2026), which is also cited below.

## The Charlson Comorbidity Index

The CCI is a single number that summarizes how much chronic illness a patient is carrying, built from 17 weighted comorbidity groups (diabetes, heart failure, kidney disease, cancer, and so on). Each group contributes a fixed weight, from 1 to 6, if the patient has a diagnosis that falls into it, and a small number of groups suppress a milder sibling group when the more severe version is also present (severe liver disease outweighs mild liver disease, metastatic cancer outweighs non-metastatic cancer, complicated diabetes outweighs uncomplicated diabetes). The CCI is widely used to risk-adjust mortality studies, to compare patient populations across hospitals, and to control for baseline health in almost any clinical or health-services research design.

`miCCI` currently bundles 17 Charlson groups, matching the enhanced ICD-10-GM mapping described above:

| Group key | Meaning | Weight | Suppressed by |
|---|---|---|---|
| `malignancy_meta` | Metastatic malignancy | 6 | none |
| `aids` | AIDS | 6 | none |
| `liver_severe` | Severe liver disease | 3 | none |
| `dm_complicated` | Diabetes with end-organ complications | 2 | none |
| `kidney` | Kidney disease | 2 | none |
| `malignancy_nonmeta` | Non-metastatic malignancy | 2 | `malignancy_meta` |
| `plegia` | Plegia / paralysis | 2 | none |
| `liver_mild` | Mild liver disease | 1 | `liver_severe` |
| `dm_simple` | Diabetes without complications | 1 | `dm_complicated` |
| `mi` | Myocardial infarction | 1 | none |
| `hf` | Heart failure | 1 | none |
| `cerebrovascular` | Cerebrovascular disease | 1 | none |
| `pulmo` | Chronic pulmonary disease | 1 | none |
| `rheumatic` | Rheumatic disease | 1 | none |
| `dementia` | Dementia | 1 | none |
| `peptic_ulcer` | Peptic ulcer disease | 1 | none |
| `peripheral_vascular` | Peripheral vascular disease | 1 | none |

## Conceptual Overview

![miCCI concept flowchart: truncated ICD-10-GM codes and Destatis subcode frequencies feed four strategies (S1 through S4), which a meta learner combines into a single CCI estimate](man/figures/miCCI_concept_flowchart.svg)

Read it top to bottom. Two things go in at the top: the truncated codes for one encounter, and the national reference frequencies (Destatis by default, or your own table). The truncated codes alone are enough to run S1, the deterministic interval, because S1 only needs to know which four-character children of a prefix *exist*, not how common each one is. S2, S3, and S4 need more: they need to know *how likely* each hidden child was, so they also draw on the Destatis frequencies (the dashed arrow). All four strategies run independently and land on four different point estimates. The meta learner at the bottom does not average them naively. It fits non-negative weights, cross-validated across ten folds, so that whichever strategy has actually been most reliable on your data gets the most say, and the weights are guaranteed to sum to one. What comes out the bottom is the combined CCI estimate.

## Estimation Strategies

If you only remember one table from this document, make it this one.

| # | Name | What it actually does | What it costs | When to reach for it |
|---|---|---|---|---|
| **S1** | Interval | Works out the honest range: the lowest possible CCI (only counting groups that *every* hidden child would trigger) and the highest possible CCI (counting every group that *any* hidden child could trigger). Reports the midpoint too. | Practically free, no randomness | You want a hard, defensible bound with zero assumptions about which code is more likely |
| **S2** | Probabilistic | Weighs every possible hidden child by how often it actually occurs nationally, and works out the *expected* CCI exactly, using inclusion-exclusion so overlapping and mutually-exclusive groups are handled correctly rather than double-counted | One pass, no simulation, exact math | You want the single best point estimate and you trust the reference frequencies |
| **S3** | Multiple Imputation | Draws a plausible full-length code 20 times (by default) from the empirical frequency distribution, scores each imputed encounter with the ordinary CCI algorithm, and averages the 20 scores | 20 passes, needs a seed | You want an estimate that behaves like "what a real coder would have written," with imputation uncertainty you can quantify by looking at the spread across draws |
| **S4** | Bayesian | Same idea as S3, but instead of trusting the empirical frequencies as fixed truth, it places a Dirichlet prior around them and draws from the resulting posterior 25 times (by default), reporting the posterior median | 25 passes, needs a seed and a prior strength (`alpha_0`) | You want built-in smoothing for prefixes with sparse reference data, or you want a genuinely Bayesian answer with a tunable prior |
| **Meta** | Cross-validated ensemble | Fits non-negative least-squares weights over S1 through S4 using 10-fold cross-validation (via the Super Learner framework), so the combination is never overfit to the encounters it was trained on | Needs the `SuperLearner` package and enough encounters for cross-validation to be meaningful | You have a large cohort and want the best achievable point estimate, letting the data decide which of the four base strategies to trust |

On the 507,949-encounter validation cohort described in the accompanying study, the meta learner put almost all of its weight on S3 (multiple imputation), with S1 picking up a small remainder and S2 and S4 contributing almost nothing, because S3 and S4 are highly correlated and the optimizer only needs one of them. Your own cohort may weigh things differently. That is the entire point of fitting the weights rather than hard-coding them.

## Configuration Parameters

Every function in this package sits at the intersection of two independent choices, and mixing them up is the single most common source of confusion, so it gets its own section before anything else.

**Do you know each patient's age?**
If yes, pass it, and miCCI conditions the subcode probabilities on the patient's Destatis age band (22 bands, from "under 1" to "95 and older"). This matters because comorbidity severity is not age-neutral: a 35-year-old with a truncated `E11` diabetes code is statistically far more likely to have the uncomplicated form than a 78-year-old with the same truncated code. If you do not have age, or choose not to use it, miCCI falls back to the marginal (age-aggregated) distribution automatically. Nothing breaks either way.

**Do you want the German national reference, or your own?**
The default is Destatis table 23131-01, the Federal Statistical Office's published national inpatient frequency table, downloaded once and cached for the session. If you are working outside Germany, or you have a better institutional reference (your own hospital's historical diagnosis frequencies, for instance), you can hand miCCI any table with the same four required columns instead. Every single-encounter function and every batch function accepts a `freq_table` argument for exactly this purpose.

These two choices are independent and combine into four concrete recipes, which is why the four use cases below walk through all four combinations with real, run numbers for each.

## Requirements

`miCCI` needs:

- **R ≥ 4.1.0** (declared in `DESCRIPTION`). Development and the examples in this document that use real live data were run against **R 4.6.0 (2026-04-24, "ucrt")**, matching `R.version.string`. Everything was additionally re-verified in this session against R 4.3.3 to confirm nothing in the package relies on anything newer than the declared minimum.
- **Hard dependencies (`Imports`):** `data.table`, `jsonlite`.
- **Optional dependencies (`Suggests`), each needed only for the feature it powers:**

  | Package | Unlocks |
  |---|---|
  | `readxl` | `load_destatis()`: parsing the Destatis Excel workbook |
  | `SuperLearner` | `cci_meta_fit()`: the cross-validated meta learner |
  | `testthat (>= 3.0.0)` | running the test suite |

None of the five core strategies (S1 through S4 and gold-standard scoring) require anything beyond `data.table` and `jsonlite`. You can compute a CCI reconstruction with nothing else installed.

## Installation

```r
# install.packages("remotes")
remotes::install_github("gaetankamdje-wabo/miCCI")

# or, from a local clone or extracted archive:
remotes::install_local("path/to/miCCI")
```

This was run end to end in this session: `remotes::install_local()` against the exact source in this repository built cleanly and installed as `miCCI 1.1.0` with no errors.

Install the optional pieces you actually plan to use:

```r
install.packages(c("readxl", "SuperLearner"))
```

## Loading

```r
library(miCCI)

quan_map <- load_quan_map()
length(quan_map)
#> [1] 17
```

If that returns `17`, the package is installed correctly and the bundled Charlson mapping loaded without a hitch.

## Testing

```r
library(testthat)
library(miCCI)
test_dir(system.file("tests", "testthat", package = "miCCI"))

# or, from a package source checkout, either of:
devtools::test()
testthat::test_dir("tests/testthat")
```

This was run in this session against the real test suite: **60 `test_that()` blocks, 142 individual expectations, 0 failures.**

One thing worth knowing before you run it yourself: on a system whose default locale is not UTF-8 (some minimal Linux containers, for instance), you may see warnings like `unable to translate '95 u. älter' to native encoding`. That is a locale artifact from the German "ä" in one of the 22 Destatis age-band labels, not a bug in miCCI, and it does not cause a test failure. If it bothers you, set a UTF-8 locale before loading the package:

```r
Sys.setlocale("LC_ALL", "en_US.UTF-8")   # or "C.UTF-8" on Linux
```

## Quick Start

This block was run in this session against the **real, live Destatis 23131-01 table** (7,684 four-character ICD-10-GM codes across 1,324 distinct three-character prefixes, downloaded fresh from `destatis.de`) and the full bundled 17-group Charlson mapping. Every number below is exactly what came back.

```r
library(miCCI)

# 1) Charlson mapping (Quan et al. 2005, German ICD-10-GM adaptation).
quan_map <- load_quan_map()

# 2) National reference frequencies. First call downloads and caches.
freq <- load_destatis()
#> Destatis loaded: 7684 codes (1324 distinct three-character prefixes)

# 3) Precompute per-prefix lookups. One-time cost per session.
cache <- precompute_lookups(freq, quan_map)
#> Precomputed lookups for 1324 ICD-3 prefixes   (took ~6.5s)

# 4) One encounter, all four strategies.
codes_anon <- c("E11", "N18", "I10")   # truncated diabetes, kidney, hypertension

cci_interval(codes_anon, quan_map, cache)
#> $cci_min          2
#> $cci_max          4
#> $cci_mid          3
#> $interval_width   2

cci_probabilistic(codes_anon, quan_map, cache)$e_cci
#> 3.586574

cci_mi(codes_anon, quan_map, cache, m = 20L, seed = 42L)$mi_cci
#> 3.6

cci_bayesian(codes_anon, quan_map, cache, n_draws = 25L, alpha_0 = 10, seed = 42L)$posterior_median
#> 4
```

Read that as: this encounter's true CCI, before truncation, is guaranteed to sit somewhere between 2 and 4 (S1). The single best point estimate under the national frequency distribution is 3.59 (S2). Twenty plausible reconstructions of the original codes average out to 3.6 (S3). A Bayesian posterior over twenty-five draws lands on 4 (S4). All four are consistent with each other and with the S1 bound, which is exactly the sanity check S1 exists to provide.

## Full-Length Codes

Everything above solves the harder problem: reconstructing a CCI from codes that have already been cut down. If your diagnosis codes are already complete, four or five characters, for example `E11.42` or `N18.4`, none of that machinery is needed. The Charlson score is then a fixed calculation with a single, exact answer. No reference population, no age, no randomness.

```r
library(miCCI)

quan_map <- load_quan_map()

cci_gold(c("E11.40", "N18.4", "I10.00"), quan_map)
#> $cci      4
#> $active   dm_complicated=2, kidney=2
```

For a full cohort, `cci_gold_batch()` does the same thing across every row of a `data.table` with a pipe-separated `diagnosen` column, and is two to three orders of magnitude faster than calling `cci_gold()` in a loop:

```r
cohort <- data.table::data.table(
  falnr     = 1:3,
  diagnosen = c("E11.4|N18.4|I10.0", "K70.4|K74.6", "C34.1|C78.0")
)
cci_gold_batch(cohort, quan_map)
#> 4  3  6
```

`cci_gold()` and `cci_gold_batch()` are also what every strategy above calls internally once it has settled on an imputed or bounded code set. The Charlson logic itself, including the hierarchical suppression rules, is identical whether you are scoring a complete code or a hypothesis about what a truncated one might have been.

## Use Case A: Destatis, No Age

This is what you get if you call every function with its defaults: the national reference population, aggregated across all ages.

**Real anecdote from the live data:** the reason the package treats a prefix with zero matching children in the reference table as "certain, not zero" rather than silently dropping it is AIDS. Confirmed live in this session:

```r
freq <- load_destatis()
quan_map <- load_quan_map()
cache <- precompute_lookups(freq, quan_map)

freq[code3 %in% c("B20", "B21", "B22", "B24")]
#> Empty data.table (0 rows)   -- Destatis genuinely has NO 4-character
#>                                 children for any AIDS code

prefix_pool("B20", cache)
#>   code_nodot  prob
#> 1:        B20     1
attr(prefix_pool("B20", cache), "status")
#> [1] "degenerate"

dt <- data.table::data.table(diagnosen = c("B20", "B21", "B22", "B24"))
cci_gold_batch(dt, quan_map)             #> 6 6 6 6
cci_interval_batch(dt, quan_map, cache)$cci_mid   #> 6 6 6 6
cci_probabilistic_batch(dt, quan_map, cache)      #> 6 6 6 6
```

Every AIDS code correctly scores 6 (the full AIDS weight) under every strategy, because a "degenerate" prefix (one Destatis never subdivides) stands for itself with probability 1, rather than contributing nothing just because it happens to have no recorded children.

**Cohort-level walkthrough**, using a small reference table and mapping bundled with the package's own test suite, so you can copy this block verbatim and get identical numbers on your own machine without a network call:

```r
qm    <- miCCI:::.test_quan()      # 17-group synthetic mapping, same shape as the real one
freq  <- miCCI:::.synth_freq()     # tiny reference table, same 4 required columns as Destatis
cache <- precompute_lookups(freq, qm)

cohort <- data.table::data.table(
  falnr = 1:7,
  diagnosen = c(
    "E11.4|N18.4|I10.0",                          # diabetes + kidney + hypertension
    "E11.9|E11.4",                                 # uncomplicated + complicated diabetes together
    "K70.4|K74.6",                                 # severe + mild liver disease together
    "C34.1|C78.0",                                 # primary + metastatic cancer together
    "I10.0|E78.0",                                 # hypertension + lipid disorder, no CCI hit
    "I21.0|I50.0|J44.0|E11.4|N18.4|I63.3|G81.1",    # complex multimorbid patient
    "E11|N18|I10|E78"                              # already-truncated codes
  ))

cohort[, cci_gold := cci_gold_batch(cohort, qm)]
#>  4  2  3  6  0  10  4

s1 <- cci_interval_batch(cohort, qm, cache)
s1$cci_mid
#>  3  1  2  6  0   9  3

cci_probabilistic_batch(cohort, qm, cache)
#>  3.428571  1.673469  1.510204  6.000000  0.000000  9.428571  3.428571

cci_mi_batch(cohort, qm, cache, m = 20L, seed = 42L)
#>  3.50  1.60  1.40  6.00  0.00  9.20  3.45

cci_bayesian_batch(cohort, qm, cache, n_draws = 25L, alpha_0 = 10, seed = 42L)
#>  3  2  1  6  0  9  3
```

Two things to notice, both by design and both checked by the test suite. Row 2 (`E11.9|E11.4`) scores 2, not 3, because the complicated-diabetes group suppresses the uncomplicated one once both are triggered on the same encounter. Row 3 (`K70.4|K74.6`) scores 3, not 4, for the identical reason applied to liver disease. And in every row, `s1_min ≤ cci_gold ≤ s1_max`: the S1 interval never fails to contain the true score, because that containment guarantee is exactly what S1 is built to provide.

## Use Case B: Destatis, With Age Stratification

Pass `age` (single-encounter functions) or a precomputed `age_idx` vector (batch functions, for speed) and every probability-based strategy conditions on the patient's Destatis age band instead of the national marginal. Run live against the real table, isolating one prefix (`E11`, diabetes) to make the effect easy to see:

```r
freq <- load_destatis()
quan_map <- load_quan_map()
cache <- precompute_lookups(freq, quan_map)

cci_probabilistic("E11", quan_map, cache)$e_cci              # no age: marginal
#> 1.586574

cci_probabilistic("E11", quan_map, cache, age = 35)$e_cci     # 35-year-old
#> 1.232791

cci_probabilistic("E11", quan_map, cache, age = 78)$e_cci     # 78-year-old
#> 1.647351
```

The reconstructed diabetes severity climbs with age, exactly as clinical experience would predict, because the underlying donor pools genuinely shift:

```r
prefix_pool("E11", cache, age_idx = age_to_bin_index(35))
#>   code_nodot       prob
#> 1:       E119 0.565    <- uncomplicated diabetes dominates in a 35-year-old
#> 2:       E117 0.188
#> 3:       E116 0.125

prefix_pool("E11", cache, age_idx = age_to_bin_index(78))
#>   code_nodot       prob
#> 1:       E117 0.530    <- complicated diabetes dominates in a 78-year-old
#> 2:       E119 0.190
#> 3:       E116 0.129
```

**Single-encounter functions take `age` directly.** `cci_probabilistic()`, `cci_mi()`, and `cci_bayesian()` each accept a raw numeric `age` and convert it internally. **Batch functions take a precomputed `age_idx` instead**, because converting age to a Destatis band index is vectorised and only needs doing once per cohort, not once per row inside a loop:

```r
dt <- data.table::data.table(
  diagnosen = c("E11|N18|I10", "K70|K74"),
  age       = c(35, 78)
)
age_idx <- age_to_bin_index(dt$age)     # do this once, outside the batch call

cci_probabilistic_batch(dt[, .(diagnosen)], quan_map, cache, age_idx = age_idx)
```

If age is missing for a given encounter, that encounter's `age_idx` is `NA`, and miCCI silently falls back to the marginal distribution for that row only.

## Use Case C: Custom Reference Population, No Age Bands

Anywhere the examples above call `load_destatis()`, you can substitute your own `data.table` instead, as long as it carries the four required columns: `code` (dotted four-character ICD-10-GM), `code_nodot` (the same without the dot), `code3` (the three-character prefix), and `freq_total` (the encounter count for that code in your population). Age-band columns are entirely optional.

```r
library(data.table)

my_freq <- data.table(
  code        = c("E11.4", "N18.4", "I10.0"),
  code_nodot  = c("E114",  "N184",  "I100"),
  code3       = c("E11",   "N18",   "I10"),
  freq_total  = c(123456,   45678,  789012)
)

quan_map <- load_quan_map()
cache    <- precompute_lookups(my_freq, quan_map)

cci_interval(c("E11", "N18", "I10"), quan_map, cache)
#> $cci_min          4
#> $cci_max          4
#> $cci_mid          4
#> $interval_width   0

cci_probabilistic(c("E11", "N18", "I10"), quan_map, cache)$e_cci
#> 4
```

The interval collapses to a single point (min = max = 4, width 0) for an honest reason: each of the three prefixes has exactly one child listed in this tiny reference table, so there is no ambiguity left to bound. That is the correct, verifiable behaviour, not a coincidence. With a richer table, more than one child per prefix, the interval widens again exactly as in Use Case A.

## Use Case D: Custom Reference Population, With Age Bands

Add columns named after the 22 Destatis age-band labels (`"unter 1"`, `"1 - 5"`, `"5 - 10"`, ... `"90 - 95"`, `"95 u. älter"`) and age-conditioning switches on automatically, exactly as it does for Destatis itself, because the loader and the age-conditioned prior both read from the same single set of band labels.

```r
library(data.table)

age_bands <- miCCI:::.MICCI_AGE_BANDS   # the 22 official Destatis band names, in order

mk_row <- function(code, code_nodot, code3, freq_total, band_vals) {
  row <- as.list(setNames(rep(0, length(age_bands)), age_bands))
  row[names(band_vals)] <- band_vals
  c(list(code = code, code_nodot = code_nodot, code3 = code3,
         freq_total = freq_total), row)
}

freq_age <- rbindlist(list(
  mk_row("E11.4", "E114", "E11", 1000, list("45 - 50" = 10,  "70 - 75" = 990)),  # complicated: mostly elderly
  mk_row("E11.9", "E119", "E11", 1000, list("45 - 50" = 990, "70 - 75" = 10))    # uncomplicated: mostly middle-aged
))

cache <- precompute_lookups(freq_age, quan_map)

cci_probabilistic("E11", quan_map, cache)$e_cci               # no age passed: marginal
#> 1.5

cci_probabilistic("E11", quan_map, cache, age = 47)$e_cci      # falls in the "45 - 50" band
#> 1.01

cci_probabilistic("E11", quan_map, cache, age = 72)$e_cci      # falls in the "70 - 75" band
#> 1.99
```

With the marginal prior, complicated and uncomplicated diabetes are each 50% likely, so the expected CCI sits exactly halfway between weight 2 and weight 1, at 1.5. Conditioning on age 47 pulls the estimate down toward the uncomplicated weight of 1, because in this reference table that age band is 99% uncomplicated. Conditioning on age 72 pulls it up toward the complicated weight of 2, for the mirror-image reason. This is precisely the mechanism that produced the real Destatis result in Use Case B, just built by hand here so you can see every input and verify the arithmetic yourself.

If an age band exists in your table's column headers but happens to carry zero recorded cases for a specific prefix (a genuinely empty cell, not a missing column), miCCI detects that and falls back to the marginal for that one (prefix, age band) combination rather than dividing by zero. This is checked directly in the test suite (`"an empty age band falls back to the marginal"`).

## Function Reference

Every exported function, organised the way the source files are, with a short description, its signature, and how to call it. Arguments shared across a whole family (`quan_map`, `cache`) are only spelled out once per family to keep this readable.

### Charlson Mapping

#### `load_quan_map(path = NULL)`
Loads the bundled ICD-10-GM to Charlson mapping from `inst/extdata/codes_quan_orig.json`. Pass `path` to point at your own curated mapping file with the same JSON schema instead.
```r
quan_map <- load_quan_map()
```
Returns a named list, one entry per Charlson group, each holding `name`, `weight`, `codes` (the matching ICD-10-GM patterns), and optionally `depends_on` (the parent group(s) that suppress this one).

#### `cci_gold(icd_codes, quan_map)`
The reference (gold-standard) CCI for one encounter's **full-length** codes, applying hierarchical suppression. Use this for exploration and unit tests. For a whole cohort, use `cci_gold_batch()` instead, which is two to three orders of magnitude faster and numerically identical.
```r
cci_gold(c("E11.40", "N18.4", "I10.00"), quan_map)
#> $cci      4
#> $active   dm_complicated=2, kidney=2
```

#### `cci_gold_batch(dt, quan_map, pattern_lookup = NULL, dep_lookup = NULL)`
Vectorised gold CCI for a full cohort. `dt` needs one character column, `diagnosen`, with pipe-separated codes. Pass a precomputed `pattern_lookup` / `dep_lookup` (see below) if you are calling this many times in a loop, to skip rebuilding them each call.
```r
cci_gold_batch(cohort, quan_map)
```
Returns an integer vector, one CCI per row of `dt`, in row order.

### Reference Data and Lookup Engine

#### `load_destatis(url = "<official Destatis URL>", cache = TRUE, force = FALSE, quiet = TRUE)`
Downloads and parses Destatis table 23131-01. The first call in a session hits the network. Every subsequent call in the same session returns the cached table for free, unless you pass `force = TRUE` (useful the day Destatis publishes an updated annual release).
```r
freq <- load_destatis()
```
Returns a `data.table` with columns `code`, `code_nodot`, `code3`, `freq_total`, and 22 age-band columns.

#### `precompute_lookups(freq_table, quan_map)`
Builds the per-prefix cache every strategy reads from: for each three-character prefix, which Charlson groups it certainly/possibly triggers, and the (optionally age-stratified) probability distribution over its children. This is the "one-time, roughly 30 seconds on the full Destatis table" step (6.5 seconds in this session's run). Call it once per `freq_table`, reuse the resulting `cache` for every subsequent strategy call.
```r
cache <- precompute_lookups(freq, quan_map)
```

#### `prefix_pool(prefix_3char, cache, age_idx = NULL)`
The lower-level building block behind every probabilistic strategy: "what could this three-character prefix have been, and with what probability?" Returns a `data.table` with `code_nodot` and `prob`, tagged with a `status` attribute: `"resolved"` (the reference table lists children, draw from them), `"degenerate"` (no children listed, so the prefix stands for itself with probability 1, as with AIDS above), or `"empty"` (blank input). Mostly useful for debugging or building your own strategy on top of the cache.
```r
prefix_pool("E11", cache, age_idx = age_to_bin_index(72))
```

#### `age_to_bin_index(age)`
Vectorised: converts a numeric age (or vector of ages) into the 1-based index of its Destatis age band, `NA` for missing/unparseable input. Every batch strategy's `age_idx` argument expects exactly this.
```r
age_to_bin_index(c(0.5, 47, 95, NA))
#> 1  11  22  NA
```

### Core Strategies

Each strategy has a single-encounter form (for exploration) and a `_batch` form (for a cohort). The single-encounter form always defers internally to the batch form on a length-1 input, so the two can never numerically disagree.

#### `cci_interval(icd_anon, quan_map, cache)` / `cci_interval_batch(dt, quan_map, cache, return_group_prob = FALSE)`
S1. Returns `cci_min`, `cci_max`, `cci_mid`, `interval_width` (batch: one row per encounter). Set `return_group_prob = TRUE` to also get a per-(encounter, group) probability table, useful for auditing which Charlson groups drove a given estimate.

#### `cci_probabilistic(icd_anon, quan_map, cache, age = NULL, preserve_multiplicity = TRUE)` / `cci_probabilistic_batch(dt, quan_map, cache, return_group_prob = FALSE, age_idx = NULL, preserve_multiplicity = TRUE)`
S2. Returns a single expected CCI (`e_cci`, single-encounter form returns it inside a list alongside `group_prob`). `preserve_multiplicity = TRUE` (the default) means a prefix coded twice on the same encounter is treated as two independent diagnosis positions, which correctly raises the expected score versus collapsing them into one. This is verified directly in the test suite.

#### `cci_mi(icd_anon, quan_map, cache, age = NULL, m = 20L, seed = 42L, preserve_multiplicity = TRUE)` / `cci_mi_batch(dt, quan_map, cache, m = 20L, seed = 42L, return_group_count = FALSE, age_idx = NULL, preserve_multiplicity = TRUE, return_rounds = FALSE)`
S3. `m` is the number of imputations, averaged for the final score. `return_rounds = TRUE` (batch only) returns the full `n × m` matrix of per-round scores rather than the average, useful for a sensitivity analysis across different values of `m`, because round `r` consumes the same random draws regardless of `m`, so the first `k` columns of an `m`-round run are exactly the `k`-round run.

#### `cci_bayesian(icd_anon, quan_map, cache, age = NULL, n_draws = 25L, alpha_0 = 10, seed = 42L, preserve_multiplicity = TRUE, aggregator = c("median","mean"))` / `cci_bayesian_batch(dt, quan_map, cache, n_draws = 25L, alpha_0 = 10, seed = 42L, return_group_count = FALSE, age_idx = NULL, preserve_multiplicity = TRUE, aggregator = c("median","mean"), return_draws = FALSE)`
S4. `alpha_0` controls how strongly the Dirichlet prior is pulled toward the reference frequencies: large values make S4 behave like S3, small values make it more dispersed. With the default `median` aggregator and an odd `n_draws`, S4's output is always an exact integer, sitting on the same support as the gold-standard score.

### Meta Learner

#### `cci_meta_fit(dt, V = 10L, seed = 42L, verbose = FALSE)`
Fits the cross-validated NNLS ensemble over S1 through S4. `dt` needs columns `cci_gold`, `s1_min`, `s1_max`, `s1_mid`, `s2_ecci`, `s3_mi`, `s4_bayes`. Requires the `SuperLearner` package (`stop()`s with a clear message if it is not installed).
```r
meta <- cci_meta_fit(preds, V = 10L, seed = 42L)
meta$weights    # named NNLS weight per base strategy, sums to 1
meta$cv_risk    # cross-validated risk per base strategy, lower is better
meta$predictions
```
On a small (56-row, heavily repeated) demonstration cohort run in this session, the weights came out `S1_mid = 0.36, S2_ecci = 0.64, S3_mi = 0, S4_bayes = 0`, which is not a meaningful result on its own. NNLS weights need a real cohort with genuine variation to be interpretable, and that is exactly what the 507,949-encounter validation cohort in the accompanying study provides (`S3_mi ≈ 0.92, S1_mid ≈ 0.06, S2_ecci ≈ 0.02, S4_bayes ≈ 0.00` there).

## Reproducibility

Every function in this package that uses randomness (`cci_mi*`, `cci_bayesian*`, `cci_meta_fit`) takes an explicit `seed` argument and wraps its random draws so that your global `.Random.seed` is captured before the call and restored afterward, whatever happens inside. Calling `cci_mi()` inside a larger simulation you are running never silently perturbs the rest of that simulation's random stream. Two consequences worth knowing: the same `seed` always reproduces the same output, and, for both S3 and S4, round `r` of an `m`-round (or `n_draws`-round) run consumes exactly the same random numbers regardless of what `m` is set to, so the first `k` rounds of a large run are numerically identical to a standalone `k`-round run. This is what makes a sensitivity sweep across different values of `m` or `n_draws` cheap: run the largest value once with `return_rounds = TRUE` / `return_draws = TRUE`, then slice the matrix.

## Implementation Notes

A few behaviours are easy to mistake for defects the first time you see them. All are verified directly in the package's own test suite, and each exists for a specific, checkable reason.

**A prefix absent from the reference table scores as certain, not zero.** AIDS (`B20`–`B24`) genuinely has no four-character children in the real Destatis table, confirmed live above. Treating "no children on record" as "this diagnosis contributes nothing" would silently discard a real, codeable diagnosis. Instead, a childless (`"degenerate"`) prefix stands for itself with probability 1.

**S1's upper bound respects hierarchical suppression too, not just the lower bound.** A truncated `K70` could reach either mild liver disease (weight 1) or severe liver disease (weight 3), but never both scored independently, because no single hidden subcode belongs to both groups at once and the severe form always suppresses the mild one when both are active. The reachable maximum is therefore 3, not 1 + 3 = 4.

**S2's math accounts for mutual exclusivity exactly, not by treating "child active" and "parent active" as independent events.** When a suppressing parent group and its child both live under the same truncated prefix, a single hidden subcode can only ever belong to one of them. The exact inclusion-exclusion formula for "child active and no parent active" is used rather than an approximation that would understate the suppression.

**Repeating the same truncated prefix on one encounter raises the expected score.** `E11|E11` (diabetes coded twice) scores higher under S2 than a single `E11`, because two independent diagnosis positions genuinely carry more evidence of disease than one. This "multiplicity" behaviour can be switched off with `preserve_multiplicity = FALSE`, which restores the older behaviour of collapsing repeated prefixes to one.

**Half-integer outputs round away from zero, never to even.** R's built-in `round()` uses banker's rounding (round-half-to-even), which would push S1's half-integer midpoints disproportionately onto even CCI values and create a fake-looking pattern in any downstream summary. Every binning operation in this package uses a "round half away from zero" convention instead.

## Citation

If you use `miCCI`, please cite the methodological study describing the five reconstruction strategies and their validation:

```
Kamdje Wabo G, Sokolowski PP, Jannesari Ladani M, Santhanam N, Hagmann M,
Ganslandt T, Siegel F.
Estimating the Charlson Comorbidity Index From Privacy-Truncated Diagnosis
Codes: Methodological Study Evaluating Computational Strategies in 507,949
Inpatient Encounters.
JMIR Preprints. 14/07/2026:107041.
DOI: 10.2196/preprints.107041
URL: https://preprints.jmir.org/preprint/107041
```

The ICD-10-GM to Charlson group mapping bundled with the package (`inst/extdata/codes_quan_orig.json`) follows the German-specific adaptation described in:

```
Sokołowski PP, Hagmann M, Maros ME, Kamdje Wabo G, Meerjanssen JM, Siegel F.
Developing Country-Specific Charlson Comorbidity Index Mappings for Use With
German Administrative Data: Methodological Comparative Study.
JMIR Med Inform. 2026;14:e93923.
DOI: 10.2196/93923
PMID: 42593349
PMCID: 13586706
```

Both citations are also available from within R once the package is installed:

```r
citation("miCCI")
```

## License

MIT. See `LICENSE` for the full text.
