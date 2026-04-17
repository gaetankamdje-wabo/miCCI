# =============================================================================
# miCCI v1.0.0 — Cohort Descriptive Analysis (single consistent cohort)
# -----------------------------------------------------------------------------
# Operates on the unsplit 2010-2024 cohort, matching the v1.0.0 pipeline.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(arrow); library(ggplot2)
})

DATA_PATH  <- "A:/HLZ/Projekt/2025/comorbidity_gaetan/comorbidity_diagnoses_only_pseudonym.parquet"
OUTPUT_DIR <- "A:/HLZ/Promotionen/Gaetan Kamdje Wabo/Results CCI Model v050/cohort_descriptive"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("=== miCCI v1.0.0 — Cohort Descriptive Analysis ===\n\n")

cat("[1/4] Loading parquet...\n")
dt <- as.data.table(read_parquet(DATA_PATH))
dt <- dt[, .(falnr, age, date_admission, date_discharge, stay_in_days, diagnosen)]
dt[, date_admission := as.Date(as.character(date_admission))]
dt[, year := as.integer(format(date_admission, "%Y"))]
dt[, diagnosen := as.character(diagnosen)]
dt <- dt[!is.na(date_admission) & !is.na(year)]
dt <- dt[!is.na(diagnosen) & diagnosen != "" & diagnosen != "NA"]
dt[, stay_in_days := as.numeric(stay_in_days)]
dt <- dt[!is.na(stay_in_days) & stay_in_days >= 0]
dt[, age := as.numeric(age)]
dt <- dt[!is.na(age) & age >= 0]
dt <- dt[date_admission >= as.Date("2010-01-01") &
         date_admission <= as.Date("2024-09-30")]
dt[, n_diagnoses := lengths(strsplit(diagnosen, "\\|+"))]
dt <- unique(dt, by = "falnr")
cat(sprintf("  loaded %d encounters (%d-%d)\n\n",
            nrow(dt), min(dt$year), max(dt$year)))

breaks_age <- c(0, 18, 30, 40, 50, 60, 70, 80, 90, Inf)
labels_age <- c("<18","18-29","30-39","40-49","50-59","60-69","70-79","80-89",">=90")
dt[, age_class := cut(age, breaks = breaks_age, labels = labels_age,
                      right = FALSE, include.lowest = TRUE)]

cat("[2/4] Summary tables\n")
.desc <- function(x) list(
  mean = round(mean(x, na.rm=TRUE), 2),
  sd   = round(sd(x, na.rm=TRUE),   2),
  median = round(median(x, na.rm=TRUE), 2),
  p25  = round(quantile(x, 0.25, na.rm=TRUE), 2),
  p75  = round(quantile(x, 0.75, na.rm=TRUE), 2)
)
age_s <- .desc(dt$age); icd_s <- .desc(dt$n_diagnoses); los_s <- .desc(dt$stay_in_days)

tab_overall <- data.table(
  Cohort           = "Mannheim 2010-2024 (consistent)",
  N_encounters     = nrow(dt),
  Year_range       = sprintf("%d - %d", min(dt$year), max(dt$year)),
  Age_mean_sd      = sprintf("%.1f +- %.1f", age_s$mean, age_s$sd),
  Age_median_iqr   = sprintf("%.0f [%.0f-%.0f]", age_s$median, age_s$p25, age_s$p75),
  ICD_mean_sd      = sprintf("%.1f +- %.1f", icd_s$mean, icd_s$sd),
  ICD_median_iqr   = sprintf("%.0f [%.0f-%.0f]", icd_s$median, icd_s$p25, icd_s$p75),
  LOS_mean_sd      = sprintf("%.1f +- %.1f", los_s$mean, los_s$sd),
  LOS_median_iqr   = sprintf("%.0f [%.0f-%.0f]", los_s$median, los_s$p25, los_s$p75)
)
fwrite(tab_overall, file.path(OUTPUT_DIR, "tab_cohort_01_overall.csv"))

tab_age  <- dt[, .N, by = age_class][order(age_class)][, pct := round(100*N/sum(N), 2)]
fwrite(tab_age,  file.path(OUTPUT_DIR, "tab_cohort_02_age_distribution.csv"))

tab_year <- dt[, .N, by = year][order(year)]
fwrite(tab_year, file.path(OUTPUT_DIR, "tab_cohort_03_yearly_volume.csv"))

cat("[3/4] Figures\n")
theme_micci <- function(base = 13) theme_minimal(base_size = base) %+replace% theme(
  panel.grid.minor = element_blank(),
  plot.title = element_text(face = "bold", size = base + 2, hjust = 0),
  plot.subtitle = element_text(color = "grey40", size = base, hjust = 0),
  axis.title = element_text(size = base + 1),
  axis.text  = element_text(size = base),
  legend.position = "bottom"
)
save2 <- function(p, base, w, h) {
  ggsave(file.path(OUTPUT_DIR, paste0(base, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(OUTPUT_DIR, paste0(base, ".png")), p, width = w, height = h, dpi = 300)
}

p1 <- ggplot(dt, aes(x = age)) +
  geom_histogram(bins = 50, fill = "#2C7BB6", color = "white", linewidth = 0.15, alpha = 0.9) +
  geom_vline(xintercept = median(dt$age), linetype = "dashed",
             color = "#1B1B1B", linewidth = 0.8) +
  labs(title = "Age Distribution - Consistent Cohort 2010-2024",
       subtitle = sprintf("n = %s encounters | dashed = median %.0f years",
                          format(nrow(dt), big.mark = ","), median(dt$age)),
       x = "Age (years)", y = "Encounters (n)") +
  theme_micci()
save2(p1, "fig_cohort_01_age_overall", 10, 5.5)

icd_cap <- quantile(dt$n_diagnoses, 0.99, na.rm = TRUE)
p2 <- ggplot(dt[n_diagnoses <= icd_cap], aes(x = n_diagnoses)) +
  geom_histogram(bins = min(icd_cap, 40), fill = "#1A9641", color = "white",
                 linewidth = 0.15, alpha = 0.9) +
  geom_vline(xintercept = median(dt$n_diagnoses), linetype = "dashed",
             color = "#1B1B1B", linewidth = 0.8) +
  labs(title = "Number of ICD-10-GM Codes per Encounter",
       subtitle = sprintf("Up to 99th percentile (<= %d codes) | dashed = median",
                          as.integer(icd_cap)),
       x = "Number of ICD codes", y = "Encounters (n)") +
  theme_micci()
save2(p2, "fig_cohort_02_icd_n", 10, 5.5)

los_cap <- quantile(dt$stay_in_days, 0.99, na.rm = TRUE)
p3 <- ggplot(dt[stay_in_days <= los_cap], aes(x = stay_in_days)) +
  geom_histogram(bins = 40, fill = "#D7191C", color = "white",
                 linewidth = 0.15, alpha = 0.9) +
  geom_vline(xintercept = median(dt$stay_in_days), linetype = "dashed",
             color = "#1B1B1B", linewidth = 0.8) +
  labs(title = "Length of Stay - Consistent Cohort 2010-2024",
       subtitle = sprintf("Up to 99th percentile (<= %d days) | dashed = median",
                          as.integer(los_cap)),
       x = "Length of stay (days)", y = "Encounters (n)") +
  theme_micci()
save2(p3, "fig_cohort_03_los", 10, 5.5)

p4 <- ggplot(tab_year, aes(x = year, y = N)) +
  geom_col(fill = "#2C3E50", width = 0.7, alpha = 0.9) +
  geom_text(aes(label = format(N, big.mark = ",")), vjust = -0.4, size = 3.5) +
  scale_x_continuous(breaks = tab_year$year) +
  labs(title = "Annual Encounter Volume",
       subtitle = "Single consistent cohort - no train/calibration/validation split",
       x = "Year", y = "Encounters (n)") +
  theme_micci() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save2(p4, "fig_cohort_04_year_volume", 12, 6)

p5 <- ggplot(tab_age, aes(x = age_class, y = pct)) +
  geom_col(fill = "#56B4E9", alpha = 0.9, width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", pct)), vjust = -0.4, size = 3.5) +
  labs(title = "Age-Class Distribution - Consistent Cohort 2010-2024",
       x = "Age class (years)", y = "Encounters (%)") +
  theme_micci()
save2(p5, "fig_cohort_05_age_class", 10, 5.5)

cat("[4/4] Done.\n")
cat(sprintf("Output: %s\n", OUTPUT_DIR))
