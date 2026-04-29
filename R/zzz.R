# Suppress R CMD check NOTEs for data.table non-standard evaluation. Every
# name listed below is a data.table column name used inside `dt[, ...]` or
# `dt[, ..col]`. R CMD check has no way to know these are bound at runtime,
# so we explicitly tell it.
utils::globalVariables(c(
  ".", ":=", "..keep", "..out_cols",
  ".chapter", ".pattern",
  "MAE", "N", "Source",
  "age", "best_mae",
  "bias", "cci_bin", "cci_max", "cci_mid", "cci_min", "child_gk",
  "cn", "code", "code3", "code_nodot", "code_raw", "contribution",
  "date_admission", "diagnosen", "diff_val",
  "facet_caption", "freq", "freq_total",
  "gk", "gold_mass", "group",
  "idx", "interval_width",
  "label", "loa_hi", "loa_lo",
  "mae", "mae_hi", "mae_lo", "mean_val", "meta",
  "n", "n_active", "n_diagnoses",
  "p", "p_parent", "panel", "pct", "pct_gold_active",
  "predicted", "prob", "q_clip",
  "ratio", "ratio_s2", "ratio_s3", "ratio_s4", "ratio_type",
  "s1_max", "s1_mid", "s1_min", "s1_width",
  "s2_ecci", "s2_mass", "s3_mass", "s3_mi", "s4_bayes", "s4_mass",
  "sd_diff", "short", "stay_in_days",
  "strategy", "strategy_label", "stratum",
  "w", "weight", "year"
))
