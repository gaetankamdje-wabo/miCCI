#!/usr/bin/env python3
"""
miCCI — Python Usage Example
=============================
This script calls the miCCI R package from Python using rpy2.
It demonstrates single-patient and batch usage for each strategy.

Requirements
------------
    pip install rpy2
    R package miCCI installed (see README.md)

Run
---
    python scripts/micci_python_example.py
"""

import os
import sys

# ── 1. Check rpy2 availability ────────────────────────────────────────────────
try:
    import rpy2.robjects as ro
    from rpy2.robjects.packages import importr
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.vectors import StrVector, IntVector
    pandas2ri.activate()
except ImportError:
    sys.exit(
        "rpy2 not found. Install it with:\n"
        "    pip install rpy2\n"
        "and make sure R is on your PATH."
    )

print("=" * 60)
print(" miCCI Python Example")
print("=" * 60)

# ── 2. Load the R package ─────────────────────────────────────────────────────
try:
    micci = importr("miCCI")
except Exception as e:
    sys.exit(
        f"Could not load miCCI: {e}\n"
        "Install the package in R first:\n"
        "    devtools::install_github('PLACEHOLDER/miCCI')"
    )

dt = importr("data.table")

# ── 3. Load reference objects (needed by all strategies) ─────────────────────
print("\nLoading Quan mapping and Destatis frequencies...")
quan_map = micci.load_quan_map()
destatis = micci.load_destatis()
cache    = micci.precompute_lookups(destatis, quan_map)
pl       = micci.build_pattern_lookup(quan_map)
dl       = micci.build_dep_lookup(quan_map)
print("Ready.\n")

# ── 4. Gold CCI — single patient ──────────────────────────────────────────────
print("--- Gold CCI (full 4-digit codes) ---")
icd_full = StrVector(["E11.40", "N18.4", "I21.0"])
gold = micci.cci_gold(icd_full, quan_map)
print(f"  CCI = {list(gold.rx2('cci'))[0]}")         # 5
print(f"  Active groups: {list(gold.rx2('active').names)}")

# ── 5. S1 — Interval (single patient) ────────────────────────────────────────
print("\n--- S1: Interval ---")
icd_anon = StrVector(["E11", "N18", "I21"])
s1 = micci.cci_interval(icd_anon, quan_map, cache)
print(f"  min={list(s1.rx2('cci_min'))[0]}  "
      f"max={list(s1.rx2('cci_max'))[0]}  "
      f"mid={list(s1.rx2('cci_mid'))[0]}")

# ── 6. S2 — Probabilistic (single patient) ────────────────────────────────────
print("\n--- S2: Probabilistic ---")
s2 = micci.cci_probabilistic(icd_anon, quan_map, cache)
print(f"  E-CCI = {list(s2.rx2('e_cci'))[0]:.3f}")

# ── 7. S3 — Multiple Imputation (single patient, m=5 for speed) ──────────────
print("\n--- S3: Multiple Imputation (m=5) ---")
s3 = micci.cci_mi(icd_anon, quan_map, cache, m=5)
print(f"  MI-CCI = {list(s3.rx2('mi_cci'))[0]:.3f}")

# ── 8. S4 — Bayesian (single patient, n_draws=5 for speed) ───────────────────
print("\n--- S4: Bayesian Dirichlet (n_draws=5) ---")
s4 = micci.cci_bayesian(icd_anon, quan_map, cache, n_draws=5)
print(f"  Posterior median = {list(s4.rx2('posterior_median'))[0]:.3f}")

# ── 9. Batch mode via run_demo() ──────────────────────────────────────────────
print("\n--- run_demo() on built-in 39-patient dataset ---")
result = micci.run_demo(m=5, n_draws=5)
# Convert to pandas
import pandas as pd
df = pandas2ri.rpy2py(result)
print(df[["p_id", "cci_gold", "s1_mid", "s2_ecci", "s3_mi"]].head(10).to_string(index=False))

print("\nDone. For the full pipeline, use run_pipeline() with your Parquet file.")
