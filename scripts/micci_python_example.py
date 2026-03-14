#!/usr/bin/env python3
"""
miCCI — Python Usage Example
=============================
This script calls the miCCI R package from Python using rpy2.
It demonstrates an implmentation using synthetic data and batch usage for each strategy.
It can takes around 45 min to be completed (on Jupyter Notebook or Google Colab)
"""

# ── 1. packages Installation ────────────────────────────────────────────────
!apt-get install -y -q r-base
!pip install -q rpy2 pandas
!git clone https://github.com/gaetankamdje-wabo/miCCI.git
!Rscript -e "install.packages(c('remotes','data.table','readxl','nnls'), repos='https://cloud.r-project.org', quiet=TRUE)"
!Rscript -e "remotes::install_local('miCCI', upgrade='never')"

# ── 2. Run a local demo implementation ───────────────────────────────────────────────

import os, tempfile
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri, default_converter
from rpy2.robjects.conversion import localconverter
import rpy2.rinterface_lib.callbacks as cb
import pandas as pd
import matplotlib
matplotlib.use("Agg")          # remove this line inside a Colab/Jupyter notebook
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

cb.consolewrite_print     = lambda x: None
cb.consolewrite_warnerror = lambda x: print(x, end="")

# ── setup ─────────────────────────────────────────────────────────────────────
ro.r("library(data.table)")
micci = importr("miCCI")
qm    = micci.load_quan_map()
cache = micci.precompute_lookups(micci.load_destatis(), qm)
ro.globalenv["qm_r"]  = qm
ro.globalenv["cch_r"] = cache
print("Ready.")

def to_df(r):
    with localconverter(default_converter + pandas2ri.converter):
        return ro.conversion.rpy2py(r)

def run_r(code):
    f = tempfile.NamedTemporaryFile(mode="w", suffix=".R", delete=False)
    f.write(code); f.close()
    try:    ro.r(f'source("{f.name}")')
    finally: os.unlink(f.name)

# ── S1–S6 on built-in 39-patient test dataset ─────────────────────────────────
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

# S1 – interval (keep min/max for coverage)
s1b <- cci_interval_batch(dt_anon, qm_r, cch_r)
dt[, S1 := s1b$cci_mid]; dt[, S1_min := s1b$cci_min]; dt[, S1_max := s1b$cci_max]

# S2 – probabilistic expected value
dt[, S2 := cci_probabilistic_batch(dt_anon, qm_r, cch_r)]

# S3 – multiple imputation  (m=20 for publication)
dt[, S3 := cci_mi_batch(dt_anon, qm_r, cch_r, m=5L)]

# S4 – bayesian dirichlet   (n_draws=25 for publication)
dt[, S4 := cci_bayesian_batch(dt_anon, qm_r, cch_r, n_draws=5L)]

# S5 – XGBoost ML
if (requireNamespace("xgboost", quietly=TRUE)) {
    ml <- cci_ml_train(dt, qm_r)
    dt[, S5 := cci_ml_predict(ml, dt_anon, qm_r)]
} else {
    message("xgboost not found – S5 falls back to S3. Run: install.packages('xgboost')")
    dt[, S5 := S3]
}

# S6 – NNLS ensemble over S1-S5
mat <- as.matrix(dt[, .(S1, S2, S3, S4, S5)])
ens <- cci_ensemble_train(mat, dt$gold)
dt[, S6 := cci_ensemble(mat, ens$weights)]
""")

# Show the origninal table with truncated diagnosis list
to_df(ro.r("dt_anon[, .(p_id, list_full_icd_code, diagnosen, stay_in_days, age, gold)]")).round(2).head(10)

# Results table -  
df = to_df(ro.r("dt[,.(p_id,gold,S1,S1_min,S1_max,S2,S3,S4,S5,S6)]")).round(2)
print(df[["p_id","gold","S1","S2","S3","S4","S5","S6"]].head(15).to_string(index=False))



# ── 3. Visualization of the coverage ────────────────────────────────────────────────

# ── coverage metrics 
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

# ── style constants 
C = {"S1":"#4C72B0","S2":"#DD8452","S3":"#55A868","S4":"#C44E52","S5":"#8172B2","S6":"#937860"}
BG, DARK = "#f0f2f5", "#222222"
C_EXACT, C_W1, C_W2 = "#2ecc71", "#f39c12", "#e74c3c"

fig = plt.figure(figsize=(16, 13))
fig.patch.set_facecolor("#f8f9fa")
fig.suptitle("miCCI  ·  Strategy Coverage vs Gold CCI  (n = 39 patients)",
             fontsize=14, fontweight="bold", color=DARK, y=0.99)
gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.48, wspace=0.34,
                       top=0.94, bottom=0.06)

# ── TOP: stacked coverage bar chart
ax1 = fig.add_subplot(gs[0, :])
x   = np.arange(len(strategies))
bw  = 0.55

exact_v   = [metrics[s]["exact"]                        for s in strategies]
w1_v      = [metrics[s]["within1"] - metrics[s]["exact"] for s in strategies]
w2_v      = [metrics[s]["within2"] - metrics[s]["within1"] for s in strategies]
b1        = [e + w for e, w in zip(exact_v, w1_v)]

ax1.bar(x, exact_v, bw, color=C_EXACT, label="Exact (±0)",  zorder=3)
ax1.bar(x, w1_v,    bw, bottom=exact_v, color=C_W1, label="Within ±1", zorder=3)
ax1.bar(x, w2_v,    bw, bottom=b1,      color=C_W2, label="Within ±2", zorder=3)

# label exact % inside bar
for xi, v in zip(x, exact_v):
    if v > 5:
        ax1.text(xi, v / 2, f"{v:.0f}%", ha="center", va="center",
                 fontsize=10, fontweight="bold", color="white")

# S1 interval coverage dashed line
ax1.axhline(s1_interval_cov, color=C["S1"], lw=1.8, ls="--", zorder=4)
ax1.text(len(strategies) - 0.28, s1_interval_cov + 1.5,
         f"S1 interval coverage {s1_interval_cov:.0f}%",
         color=C["S1"], fontsize=8.5, ha="right")

ax1.set_xticks(x); ax1.set_xticklabels(strategies, fontsize=12)
ax1.set_ylabel("Patients (%)", fontsize=10)
ax1.set_ylim(0, 118)
ax1.set_title("Coverage: % patients within tolerance of Gold CCI", fontsize=11, pad=6)
ax1.legend(loc="upper left", fontsize=9, framealpha=0.9, ncol=3)
ax1.set_facecolor(BG)
ax1.grid(axis="y", color="white", lw=1.3, zorder=2)
ax1.spines[["top","right","left"]].set_visible(False)

# ── BOTTOM 6: gold vs predicted scatter per strategy 
positions = [(1,0),(1,1),(1,2),(2,0),(2,1),(2,2)]
lim = max(gold.max(), max(df[s].max() for s in strategies)) + 1.5

for i, s in enumerate(strategies):
    r, c = positions[i]
    ax = fig.add_subplot(gs[r, c])
    pred = df[s].values

    exact_m  = np.abs(gold - pred) < 0.5
    within1_m = (np.abs(gold - pred) >= 0.5) & (np.abs(gold - pred) < 1.5)
    off_m     = np.abs(gold - pred) >= 1.5

    ax.scatter(gold[exact_m],   pred[exact_m],   color=C_EXACT, s=60, zorder=4, label="Exact")
    ax.scatter(gold[within1_m], pred[within1_m], color=C_W1,   s=60, zorder=4, label="±1")
    ax.scatter(gold[off_m],     pred[off_m],     color=C_W2,   s=60, zorder=4, label="±2+")

    # S1 only: vertical error bars showing interval width
    if s == "S1":
        for _, row in df.iterrows():
            ax.plot([row["gold"], row["gold"]], [row["S1_min"], row["S1_max"]],
                    color=C["S1"], alpha=0.4, lw=1.3, zorder=3)

    ax.plot([0, lim], [0, lim], color=DARK, lw=1, ls="--", alpha=0.4, zorder=3)
    ax.set_xlim(-0.5, lim); ax.set_ylim(-0.5, lim)
    ax.set_xlabel("Gold CCI", fontsize=8.5)
    ax.set_ylabel("Predicted CCI", fontsize=8.5)
    ax.set_title(
        f"{s}    exact {metrics[s]['exact']:.0f}%  MAE {metrics[s]['mae']:.2f}",
        fontsize=10, fontweight="bold", color=C[s], pad=5
    )
    ax.set_facecolor(BG)
    ax.grid(color="white", lw=1.1, zorder=2)
    ax.spines[["top","right"]].set_visible(False)
    if i == 0:
        ax.legend(fontsize=8, loc="upper left", framealpha=0.85)

plt.savefig("micci_coverage.png", dpi=150, bbox_inches="tight")
plt.show()
