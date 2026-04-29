# Suppress R CMD check NOTEs for data.table non-standard evaluation. Every
# name listed below is a data.table column name used inside `dt[, ...]` or
# `dt[, ..col]`. R CMD check has no way to know these are bound at runtime,
# so we explicitly tell it.
utils::globalVariables(c(
  ".", ":=",
  "child_gk", "cci_max", "cci_mid", "cci_min",
  "cn", "code", "code3", "code_nodot", "code_raw",
  "contribution", "diagnosen", "factor",
  "freq", "freq_total", "gk", "idx", "interval_width",
  "n_active",
  "p", "p_parent", "pref", "prob", "q", "q_clip",
  "st", "w"
))
