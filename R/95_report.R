# =============================================================================
# miCCI / 95_report.R
# Stage 2: build_report().
#
# Reads the artefact written by compute_predictions() and rebuilds every
# table, plot, CI and stratified evaluation. Pure consumer: never recomputes
# strategies or touches Destatis. Wall-clock: minutes, not hours.
# =============================================================================

#' Save a ggplot in PDF (cairo) and PNG (300 dpi) at the same base path.
#' @keywords internal
.save_plot <- function(p, base_path, width, height) {
  if (!isTRUE(requireNamespace("ggplot2", quietly = TRUE)))
    stop("build_report() requires the 'ggplot2' package.")
  ggplot2::ggsave(paste0(base_path, ".pdf"), p,
                  width = width, height = height, device = cairo_pdf)
  ggplot2::ggsave(paste0(base_path, ".png"), p,
                  width = width, height = height, dpi = 300)
}

#' Strategy display: column name -> pretty label.
#' @keywords internal
.strategies_default <- function() {
  c(
    s1_mid   = "S1 Interval",
    s2_ecci  = "S2 Probabilistic",
    s3_mi    = "S3 MI-CCI",
    s4_bayes = "S4 Bayesian",
    meta     = "Meta"
  )
}

#' Standard publication theme.
#' @keywords internal
.theme_micci <- function(base = 14) {
  gg <- asNamespace("ggplot2")
  `%+replace%` <- get("%+replace%", envir = gg)
  gg$theme_minimal(base_size = base) %+replace%
    gg$theme(
      panel.grid.minor   = gg$element_blank(),
      strip.text         = gg$element_text(face = "bold", size = base + 2),
      strip.background   = gg$element_rect(fill = "#f0f2f5", color = NA),
      plot.title         = gg$element_text(face = "bold", size = base + 4, hjust = 0),
      plot.subtitle      = gg$element_text(color = "grey35", size = base, hjust = 0),
      plot.caption       = gg$element_text(color = "grey35", size = base - 2, hjust = 0),
      axis.title         = gg$element_text(size = base + 1),
      axis.text          = gg$element_text(size = base),
      legend.position    = "bottom",
      legend.title       = gg$element_blank(),
      legend.text        = gg$element_text(size = base)
    )
}

#' Stage 2: build all report artefacts from a predictions file.
#'
#' @param input_dir    directory containing predictions.parquet (or .rds)
#'   plus run_manifest.json.
#' @param output_dir   destination for tables and figures (default: input_dir).
#' @param B            bootstrap resamples for global per-strategy CIs.
#' @param B_strata     bootstrap resamples per stratum.
#' @param top_patterns number of multimorbid patterns to evaluate.
#' @param parallel     if TRUE, parallelise the bootstrap (needs `future.apply`).
#' @return invisibly, list of tables and a vector of written file paths.
#' @export
build_report <- function(input_dir,
                         output_dir   = input_dir,
                         B            = 1000L,
                         B_strata     = 500L,
                         top_patterns = 30L,
                         parallel     = FALSE) {
  if (!isTRUE(requireNamespace("ggplot2", quietly = TRUE)))
    stop("build_report() requires the 'ggplot2' package.")
  gg <- asNamespace("ggplot2")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  message("=== miCCI Stage 2 - build_report ===")

  message("[1/7] Loading predictions artefact")
  preds <- load_predictions(input_dir)
  message(sprintf("  loaded %d encounters, %d columns",
                  nrow(preds), ncol(preds)))
  quan_map <- load_quan_map()
  strategies <- .strategies_default()

  pal <- c("Gold Standard"    = "#1B1B1B",
           "S1 Interval"      = "#E69F00",
           "S2 Probabilistic" = "#56B4E9",
           "S3 MI-CCI"        = "#009E73",
           "S4 Bayesian"      = "#CC79A7",
           "Meta"             = "#6E6E6E")

  written <- character(0L)
  add <- function(p) { written <<- c(written, p); invisible(p) }
  save_csv <- function(x, name) {
    p <- file.path(output_dir, name); fwrite(x, p); add(p)
  }
  save_fig <- function(p, base, w, h) {
    bp <- file.path(output_dir, base); .save_plot(p, bp, w, h)
    add(paste0(bp, ".pdf")); add(paste0(bp, ".png"))
  }

  message("[2/7] Bootstrap CIs (global per-strategy)")
  global_ci <- bootstrap_strategies(preds, strategies,
                                    B = B, seed = 42L, parallel = parallel)
  save_csv(global_ci, "tab_01_strategy_metrics_ci.csv")

  message("[3/7] Descriptive statistics")
  desc_fn <- function(x, nm) data.table(
    Source = nm, N = sum(!is.na(x)),
    Min  = round(min(x, na.rm = TRUE), 3L),
    Q1   = round(stats::quantile(x, 0.25, na.rm = TRUE), 3L),
    Median = round(stats::median(x, na.rm = TRUE), 3L),
    Mean = round(mean(x, na.rm = TRUE), 3L),
    Q3   = round(stats::quantile(x, 0.75, na.rm = TRUE), 3L),
    Max  = round(max(x, na.rm = TRUE), 3L),
    SD   = round(stats::sd(x, na.rm = TRUE), 3L)
  )
  desc_tab <- rbind(
    desc_fn(preds$cci_gold, "Gold Standard"),
    rbindlist(lapply(names(strategies),
                     function(c) desc_fn(preds[[c]], strategies[c])))
  )
  save_csv(desc_tab, "tab_02_descriptive_statistics.csv")

  message("[4/7] QA (coverage, mass conservation, score distribution)")
  qa_cov <- qa_group_coverage(preds, quan_map)
  save_csv(qa_cov, "tab_03_qa_group_coverage.csv")
  qa_dist <- qa_score_distribution(preds, names(strategies))
  save_csv(qa_dist$frequency, "tab_04_qa_score_frequency.csv")
  save_csv(qa_dist$ks,        "tab_05_qa_ks_distance.csv")
  qa_post <- qa_posterior_coverage(preds, names(strategies))
  save_csv(qa_post, "tab_06_qa_posterior_coverage.csv")

  message("[5/7] Stratified evaluation (chapter / group / pattern)")
  strat <- stratified_evaluation(preds, strategies, quan_map,
                                 B = B_strata, top_patterns = top_patterns,
                                 parallel = parallel)
  save_csv(strat$chapter, "tab_07_stratified_chapter.csv")
  save_csv(strat$group,   "tab_08_stratified_charlson_group.csv")
  save_csv(strat$pattern, "tab_09_stratified_pattern.csv")
  save_csv(stratum_winners(strat$chapter), "tab_10_winners_chapter.csv")
  save_csv(stratum_winners(strat$group),   "tab_11_winners_group.csv")
  save_csv(stratum_winners(strat$pattern), "tab_12_winners_pattern.csv")

  message("[6/7] Generating publication plots")

  ## Plot 01 - MAE comparison with CI error bars
  d1 <- copy(global_ci)
  d1[, label := factor(label, levels = label[order(mae)])]
  p1 <- gg$ggplot(d1, gg$aes(x = label, y = mae, fill = label)) +
    gg$geom_col(width = 0.6, alpha = 0.9) +
    gg$geom_errorbar(gg$aes(ymin = mae_lo, ymax = mae_hi),
                     width = 0.2, linewidth = 0.6, color = "#1B1B1B") +
    gg$geom_text(gg$aes(label = sprintf("%.4f", mae)), vjust = -0.6, size = 5) +
    gg$scale_fill_manual(values = pal[-1L], guide = "none") +
    gg$labs(
      title    = "Mean Absolute Error by Strategy",
      subtitle = sprintf("Bootstrap 95%% CI (B=%d). Lower is better.", B),
      caption  = "Strategies sorted by point MAE. Error bars are tight at large n; numerical CIs in tab_01.",
      x = NULL, y = "MAE"
    ) +
    .theme_micci(14) +
    gg$theme(axis.text.x = gg$element_text(angle = 25, hjust = 1))
  save_fig(p1, "fig_01_mae_with_ci", 10, 6.5)

  ## Plot 02 - Predicted vs Gold (axis_max computed from the data)
  ld <- melt(preds[, c("cci_gold", names(strategies)), with = FALSE],
             id.vars = "cci_gold",
             variable.name = "strategy", value.name = "predicted")
  ld[, strategy_label := factor(strategies[as.character(strategy)],
                                levels = unname(strategies))]
  pred_max <- max(c(preds$cci_gold,
                    sapply(names(strategies), function(c) max(preds[[c]], na.rm = TRUE))),
                  na.rm = TRUE)
  axis_max <- max(20L, ceiling(pred_max / 4) * 4L)

  cap_dt <- copy(global_ci)
  cap_dt[, strategy_label := factor(label, levels = unname(strategies))]
  cap_dt[, facet_caption := sprintf("MAE = %.3f  [%.3f, %.3f]",
                                    mae, mae_lo, mae_hi)]

  p2 <- gg$ggplot(ld, gg$aes(x = cci_gold, y = predicted)) +
    gg$geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                   color = "#D7191C", linewidth = 0.8) +
    gg$geom_hex(bins = 32) +
    gg$scale_fill_viridis_c(
      trans = "log10", name = "Encounters",
      breaks = c(1, 10, 100, 1000, 10000),
      guide = gg$guide_colorbar(barwidth = grid::unit(0.6, "cm"),
                                barheight = grid::unit(5, "cm"))
    ) +
    gg$geom_text(data = cap_dt,
                 gg$aes(x = 0.5, y = axis_max - 1,
                        label = facet_caption),
                 hjust = 0, vjust = 1, size = 5, color = "#1B1B1B",
                 inherit.aes = FALSE, fontface = "bold") +
    gg$facet_wrap(~ strategy_label, ncol = 3) +
    gg$coord_fixed(xlim = c(0, axis_max), ylim = c(0, axis_max)) +
    gg$scale_x_continuous(breaks = seq(0, axis_max, 4)) +
    gg$scale_y_continuous(breaks = seq(0, axis_max, 4)) +
    gg$labs(
      title    = "Predicted vs. Gold-Standard CCI",
      subtitle = "Hex density on log10 scale; red dashed = identity.",
      x = "Gold-Standard CCI", y = "Predicted CCI"
    ) +
    .theme_micci(16) +
    gg$theme(legend.position = "right",
             legend.title = gg$element_text(size = 13, face = "bold"))
  save_fig(p2, "fig_02_predicted_vs_gold", 18, 11)

  ## Plot 03 - Temporal stability
  ym <- preds[, lapply(.SD, function(x) mean(abs(x - cci_gold), na.rm = TRUE)),
              by = year, .SDcols = names(strategies)]
  yl <- melt(ym, id.vars = "year", variable.name = "strategy", value.name = "MAE")
  yl[, label := factor(strategies[as.character(strategy)],
                       levels = unname(strategies))]
  p3 <- gg$ggplot(yl, gg$aes(x = year, y = MAE, color = label)) +
    gg$geom_line(linewidth = 1) + gg$geom_point(size = 2.5) +
    gg$scale_color_manual(values = pal[-1L]) +
    gg$scale_x_continuous(breaks = sort(unique(yl$year))) +
    gg$scale_y_continuous(limits = c(0, NA),
                          expand = gg$expansion(mult = c(0, 0.05))) +
    gg$labs(title = "Temporal Stability: MAE by Calendar Year",
            subtitle = "Lower is better. All strategies are training-free.",
            x = "Calendar year", y = "MAE") +
    .theme_micci(14) +
    gg$theme(axis.text.x = gg$element_text(angle = 30, hjust = 1))
  save_fig(p3, "fig_03_temporal_stability", 12, 6.5)

  ## Plot 04 - Bin histogram
  bin_max <- 5L
  hist_rows <- list()
  for (c in names(strategies)) {
    lbl <- strategies[c]
    g <- preds$cci_gold
    g_int <- pmin(pmax(as.integer(round(g)), 0L), bin_max)
    g_tab <- data.table(cci_bin = g_int)[, .N, by = cci_bin]
    g_tab[, pct := 100 * N / nrow(preds)]
    g_tab[, Source := "Gold Standard"]; g_tab[, panel := lbl]

    s <- preds[[c]]
    s_int <- pmin(pmax(as.integer(round(s)), 0L), bin_max)
    s_tab <- data.table(cci_bin = s_int)[, .N, by = cci_bin]
    s_tab[, pct := 100 * N / nrow(preds)]
    s_tab[, Source := lbl]; s_tab[, panel := lbl]

    hist_rows[[length(hist_rows) + 1L]] <- g_tab
    hist_rows[[length(hist_rows) + 1L]] <- s_tab
  }
  hist_dt <- rbindlist(hist_rows, use.names = TRUE, fill = TRUE)
  all_bins <- CJ(cci_bin = 0:bin_max,
                 Source  = unique(hist_dt$Source),
                 panel   = unique(hist_dt$panel))
  all_bins <- all_bins[Source == "Gold Standard" | Source == panel]
  hist_dt <- merge(all_bins, hist_dt, by = c("cci_bin", "Source", "panel"),
                   all.x = TRUE)
  hist_dt[is.na(pct), pct := 0]
  hist_dt[is.na(N), N := 0L]
  hist_dt[, panel  := factor(panel,  levels = unname(strategies))]
  hist_dt[, Source := factor(Source, levels = c("Gold Standard", unname(strategies)))]
  hist_pal <- c(pal["Gold Standard"], pal[unname(strategies)])

  p4 <- gg$ggplot(hist_dt, gg$aes(x = cci_bin, y = pct, fill = Source)) +
    gg$geom_col(position = gg$position_dodge(width = 0.75),
                width = 0.7, alpha = 0.9) +
    gg$facet_wrap(~ panel, ncol = 3) +
    gg$scale_fill_manual(values = hist_pal, drop = FALSE) +
    gg$scale_x_continuous(breaks = 0:bin_max) +
    gg$scale_y_continuous(limits = c(0, NA),
                          expand = gg$expansion(mult = c(0, 0.1)),
                          labels = function(x) paste0(x, "%")) +
    gg$labs(title = "Distributional Fidelity (CCI 0..5)",
            subtitle = sprintf("Per-bin share of cohort (n=%s).",
                               format(nrow(preds), big.mark = ",")),
            x = "CCI score", y = "Share of encounters") +
    .theme_micci(13) +
    gg$theme(legend.position = "bottom",
             panel.spacing = grid::unit(1.2, "lines"))
  save_fig(p4, "fig_04_bin_histogram", 16, 9)

  ## Plot 05 - Bland-Altman
  ba_rows <- list()
  for (c in names(strategies)) {
    ba_rows[[length(ba_rows) + 1L]] <- data.table(
      mean_val = (preds$cci_gold + preds[[c]]) / 2,
      diff_val = preds[[c]] - preds$cci_gold,
      label    = strategies[c]
    )
  }
  ba_dt <- rbindlist(ba_rows)
  ba_dt[, label := factor(label, levels = unname(strategies))]
  ba_stats <- ba_dt[, .(bias = mean(diff_val, na.rm = TRUE),
                        sd_diff = stats::sd(diff_val, na.rm = TRUE)), by = label]
  ba_stats[, loa_lo := bias - 1.96 * sd_diff]
  ba_stats[, loa_hi := bias + 1.96 * sd_diff]

  p5 <- gg$ggplot(ba_dt, gg$aes(x = mean_val, y = diff_val)) +
    gg$geom_hline(yintercept = 0, color = "grey50") +
    gg$geom_point(alpha = 0.20, size = 0.8, color = "#2C3E50") +
    gg$geom_hline(data = ba_stats, gg$aes(yintercept = bias),
                  linetype = "dashed", color = "#D55E00", linewidth = 0.8) +
    gg$geom_hline(data = ba_stats, gg$aes(yintercept = loa_lo),
                  linetype = "dotted", color = "#0072B2", linewidth = 0.7) +
    gg$geom_hline(data = ba_stats, gg$aes(yintercept = loa_hi),
                  linetype = "dotted", color = "#0072B2", linewidth = 0.7) +
    gg$facet_wrap(~ label, ncol = 3) +
    gg$labs(title = "Bland-Altman Agreement: Predicted - Gold",
            subtitle = "Dashed = bias  |  Dotted = 95% Limits of Agreement",
            x = "Mean of Gold and Predicted CCI",
            y = "Difference (Predicted - Gold)") +
    .theme_micci(14)
  save_fig(p5, "fig_05_bland_altman", 16, 10)

  ## Plot 06 - QA mass conservation (per-strategy ratios, real per-group values)
  qa_long <- melt(qa_cov, id.vars = c("group", "weight", "pct_gold_active"),
                  measure.vars = c("ratio_s2", "ratio_s3", "ratio_s4"),
                  variable.name = "ratio_type", value.name = "ratio")
  qa_long[, ratio_type := factor(ratio_type,
            levels = c("ratio_s2", "ratio_s3", "ratio_s4"),
            labels = c("S2 / Gold", "S3 / Gold", "S4 / Gold"))]
  qa_pal <- c("S2 / Gold" = unname(pal["S2 Probabilistic"]),
              "S3 / Gold" = unname(pal["S3 MI-CCI"]),
              "S4 / Gold" = unname(pal["S4 Bayesian"]))
  p6 <- gg$ggplot(qa_long, gg$aes(x = stats::reorder(group, pct_gold_active),
                                  y = ratio, fill = ratio_type)) +
    gg$geom_col(position = gg$position_dodge(width = 0.75),
                width = 0.7, alpha = 0.9) +
    gg$geom_hline(yintercept = 1, linetype = "dashed",
                  color = "#1B1B1B", linewidth = 0.6) +
    gg$scale_fill_manual(values = qa_pal) +
    gg$coord_flip() +
    gg$labs(title = "QA: Per-Group Mass Conservation",
            subtitle = "Strategy mass divided by gold mass. Ratio 1 means perfect conservation.",
            x = "Charlson group", y = "Mass ratio (Strategy / Gold)") +
    .theme_micci(13)
  save_fig(p6, "fig_06_qa_mass_conservation", 12, 9)

  ## Plot 07 - Stratified MAE by chapter
  if (nrow(strat$chapter) > 0L) {
    p7 <- gg$ggplot(strat$chapter,
                    gg$aes(x = stats::reorder(stratum, mae), y = mae, color = label)) +
      gg$geom_pointrange(gg$aes(ymin = mae_lo, ymax = mae_hi),
                         position = gg$position_dodge(width = 0.55),
                         size = 0.5, linewidth = 0.7) +
      gg$scale_color_manual(values = pal[-1L]) +
      gg$coord_flip() +
      gg$labs(title = "Stratified MAE by ICD-10 Chapter",
              subtitle = sprintf("Bootstrap 95%% CI per stratum (B=%d).", B_strata),
              x = NULL, y = "MAE") +
      .theme_micci(13)
    save_fig(p7, "fig_07_stratified_chapter", 13, 9)
  }

  ## Plot 08 - Stratified MAE by Charlson group
  if (nrow(strat$group) > 0L) {
    p8 <- gg$ggplot(strat$group,
                    gg$aes(x = stats::reorder(stratum, mae), y = mae, color = label)) +
      gg$geom_pointrange(gg$aes(ymin = mae_lo, ymax = mae_hi),
                         position = gg$position_dodge(width = 0.55),
                         size = 0.5, linewidth = 0.7) +
      gg$scale_color_manual(values = pal[-1L]) +
      gg$coord_flip() +
      gg$labs(title = "Stratified MAE by Active Charlson Group",
              subtitle = sprintf("Bootstrap 95%% CI per stratum (B=%d)", B_strata),
              x = "Charlson group", y = "MAE") +
      .theme_micci(13)
    save_fig(p8, "fig_08_stratified_charlson_group", 13, 9)
  }

  ## Plot 09 - Top multimorbid patterns
  if (nrow(strat$pattern) > 0L) {
    pat_dt <- copy(strat$pattern)
    pat_dt[, short := substr(stratum, 1L, 60L)]
    p9 <- gg$ggplot(pat_dt,
                    gg$aes(x = stats::reorder(short, mae), y = mae, color = label)) +
      gg$geom_pointrange(gg$aes(ymin = mae_lo, ymax = mae_hi),
                         position = gg$position_dodge(width = 0.55),
                         size = 0.5, linewidth = 0.7) +
      gg$scale_color_manual(values = pal[-1L]) +
      gg$coord_flip() +
      gg$labs(title = sprintf("Stratified MAE - Top %d Patterns", top_patterns),
              x = "Gold-active Charlson pattern (truncated)", y = "MAE") +
      .theme_micci(12)
    save_fig(p9, "fig_09_stratified_pattern", 14, 11)
  }

  message(sprintf("[7/7] Done. Wrote %d artefacts to %s",
                  length(written), output_dir))
  invisible(list(
    global_ci    = global_ci,
    descriptive  = desc_tab,
    qa_coverage  = qa_cov,
    qa_distribution = qa_dist,
    qa_posterior = qa_post,
    stratified   = strat,
    files        = written
  ))
}
