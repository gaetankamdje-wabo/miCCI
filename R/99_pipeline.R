# =============================================================================
# miCCI — 99_pipeline.R
# Public per-strategy wrappers, evaluation helper, full pipeline, demo runner
# ggplot2 bug fix: use asNamespace("ggplot2") instead of bare `ggplot2`
# =============================================================================

# ── Internal helpers ──────────────────────────────────────────────────────────

.load_mannheim_data <- function(path,
                                date_from = "2010-01-01",
                                date_to   = "2024-09-30") {
  if (!requireNamespace("arrow", quietly = TRUE))
    stop("Package 'arrow' is required to read Parquet files. ",
         "Install it with: install.packages('arrow')")
  dt <- as.data.table(arrow::read_parquet(path))
  need <- c("falnr","age","date_admission","date_discharge",
            "stay_in_days","diagnosen")
  miss <- setdiff(need, names(dt))
  if (length(miss)) stop("Missing columns in Parquet file: ",
                          paste(miss, collapse = ", "))
  dt <- dt[, ..need]
  dt[, date_admission := as.Date(as.character(date_admission))]
  dt[, year           := as.integer(format(date_admission, "%Y"))]
  dt[, diagnosen      := as.character(diagnosen)]
  dt <- dt[!is.na(date_admission) & !is.na(year)]
  dt <- dt[!is.na(diagnosen) & diagnosen != "" & diagnosen != "NA"]
  dt[, stay_in_days := as.numeric(stay_in_days)]
  dt <- dt[!is.na(stay_in_days) & stay_in_days >= 0]
  dt[, age := as.numeric(age)]
  dt <- dt[!is.na(age) & age >= 0]
  dt <- dt[date_admission >= as.Date(date_from) &
           date_admission <= as.Date(date_to)]
  dt[, n_diagnoses := lengths(strsplit(diagnosen, "\\|+"))]
  los_r <- frank(dt$stay_in_days, ties.method = "dense")
  dt[, los_decile := as.integer(ceiling(10 * los_r / max(los_r)))]
  dt[los_decile > 10L, los_decile := 10L]
  dt <- unique(dt, by = "falnr")
  message(sprintf("\u2705 Loaded %d encounters (%d\u2013%d)",
                  nrow(dt), min(dt$year), max(dt$year)))
  dt
}

.temporal_split <- function(dt,
                             train_years = 2010:2018,
                             cal_year    = 2019L,
                             val_years   = 2020:2024) {
  list(
    train     = dt[year %in% train_years],
    calibrate = dt[year == cal_year],
    validate  = dt[year %in% val_years]
  )
}

.timed <- function(label, expr) {
  cat(sprintf("    %s ... ", label))
  t0 <- proc.time()
  r  <- eval(expr)
  cat(sprintf("%.1fs\n", (proc.time() - t0)[3]))
  r
}

# ── Public helper: accuracy metrics ──────────────────────────────────────────

#' Compute standard accuracy metrics between estimated and gold CCI
#'
#' @description
#' Returns five metrics used throughout the paper to compare every strategy
#' against the gold standard: mean absolute error (MAE), root mean squared
#' error (RMSE), coefficient of determination (R-squared), mean bias, and
#' sample size.
#'
#' @param predicted Numeric vector of estimated CCI values from any strategy.
#' @param actual    Numeric vector of gold-standard CCI values (same length).
#' @return Named numeric vector: mae, rmse, r_squared, bias, n.
#' @examples
#' eval_metrics(c(1.2, 3.1, 0.8), c(1L, 3L, 1L))
#' @export
eval_metrics <- function(predicted, actual) {
  d   <- predicted - actual
  ss_res <- sum(d^2,                             na.rm = TRUE)
  ss_tot <- sum((actual - mean(actual, na.rm = TRUE))^2, na.rm = TRUE)
  c(
    mae       = round(mean(abs(d),           na.rm = TRUE), 4),
    rmse      = round(sqrt(mean(d^2,         na.rm = TRUE)), 4),
    r_squared = round(1 - ss_res / ss_tot,                   4),
    bias      = round(mean(d,                na.rm = TRUE), 4),
    n         = length(actual)
  )
}

# ── Public per-strategy wrappers (single patient + batch) ────────────────────
# Each wrapper prepares lookups when called standalone, so every function
# works independently without a prior setup step.

#' S1: Interval CCI for a single patient
#'
#' @description
#' Classifies every 3-character ICD prefix as *certain* (all 4-digit subcodes
#' map to the same Charlson group) or *possible* (only a subset do). Returns a
#' lower bound (CCI_min), an upper bound (CCI_max), a midpoint estimate
#' (CCI_mid), and the interval width.  No training data or external frequencies
#' are needed: this strategy works from ICD hierarchy logic alone.
#'
#' @param icd_anon  Character vector of 3-character ICD-10-GM prefixes.
#' @param quan_map  Quan mapping list, from \code{load_quan_map()}.
#' @param cache     Precomputed lookup environment, from \code{precompute_lookups()}.
#' @return Named list: cci_min (integer), cci_max (integer), cci_mid (numeric),
#'   interval_width (integer).
#' @examples
#' \dontrun{
#' qm  <- load_quan_map()
#' dst <- load_destatis()
#' cch <- precompute_lookups(dst, qm)
#' cci_interval(c("E11","N18","I50"), qm, cch)
#' }
#' @export
cci_interval <- function(icd_anon, quan_map, cache) {
  codes  <- unique(icd3(icd_anon))
  gn     <- get(".group_names", envir = cache)
  dm     <- get(".dep_map",     envir = cache)
  wm     <- get(".wt_map",      envir = cache)
  status <- setNames(rep("none", length(gn)), gn)
  for (code in codes) {
    pc <- get_prefix_cache(code, cache)
    if (is.null(pc)) next
    for (gk in gn) {
      s <- pc$group_status[gk]
      if (s == "certain")  { status[gk] <- "certain" }
      else if (s == "possible" && status[gk] == "none") { status[gk] <- "possible" }
    }
  }
  certain <- names(status)[status == "certain"]
  for (gk in certain) {
    deps <- dm[[gk]]
    if (length(deps) > 0 && any(deps %in% certain))
      certain <- setdiff(certain, gk)
  }
  any_hit  <- names(status)[status %in% c("possible", "certain")]
  cci_min  <- sum(wm[certain])
  cci_max  <- sum(wm[any_hit])
  list(
    cci_min        = as.integer(cci_min),
    cci_max        = as.integer(cci_max),
    cci_mid        = (cci_min + cci_max) / 2,
    interval_width = as.integer(cci_max - cci_min)
  )
}

#' S1: Interval CCI for a batch of patients
#'
#' @param dt       data.table with a column named \code{diagnosen} containing
#'   pipe-separated 3-character ICD-10-GM codes.
#' @param quan_map Quan mapping list.
#' @param cache    Precomputed lookup environment.
#' @return data.table with columns cci_min, cci_max, cci_mid, interval_width
#'   (one row per input row, same order).
#' @export
cci_interval_batch <- function(dt, quan_map, cache) {
  gn <- get(".group_names",    envir = cache)
  wm <- get(".wt_map",         envir = cache)
  dl <- get(".dep_lookup_dt",  envir = cache)

  all_pref <- ls(cache)
  all_pref <- all_pref[!startsWith(all_pref, ".")]
  pg_rows  <- list()
  for (pref in all_pref) {
    pc <- get(pref, envir = cache)
    for (gk in gn) {
      if (pc$group_status[gk] != "none")
        pg_rows[[length(pg_rows)+1L]] <-
          list(pref = pref, gk = gk, st = pc$group_status[gk])
    }
  }
  pg <- rbindlist(pg_rows)

  dx_list <- strsplit(as.character(dt$diagnosen), "\\|+")
  long <- data.table(
    idx  = rep(seq_len(nrow(dt)), lengths(dx_list)),
    code = toupper(gsub("[^A-Z0-9]", "", unlist(dx_list, use.names = FALSE)))
  )
  long <- long[nchar(code) >= 3]
  long[, pref := substr(code, 1, 3)]
  long <- unique(long[, .(idx, pref)])

  hits    <- merge(long, pg, by = "pref", allow.cartesian = TRUE)
  enc_grp <- hits[, .(st = if (any(st == "certain")) "certain" else "possible"),
                  by = .(idx, gk)]

  wt_dt   <- data.table(gk = names(wm), w = unname(wm))
  max_dt  <- merge(enc_grp, wt_dt, by = "gk")
  cmax_dt <- max_dt[, .(cci_max = sum(w)), by = idx]

  cert_dt <- enc_grp[st == "certain"]
  if (nrow(cert_dt) > 0 && nrow(dl) > 0) {
    sup <- merge(cert_dt[, .(idx, gk)], dl,
                 by.x = "gk", by.y = "child", allow.cartesian = TRUE)
    sup <- merge(sup, cert_dt[, .(idx, gk)],
                 by.x = c("idx","parent"), by.y = c("idx","gk"))
    if (nrow(sup) > 0)
      cert_dt <- fsetdiff(cert_dt[, .(idx,gk)], unique(sup[, .(idx,gk)]))
    else
      cert_dt <- cert_dt[, .(idx, gk)]
  } else {
    cert_dt <- cert_dt[, .(idx, gk)]
  }
  min_wt  <- merge(cert_dt, wt_dt, by = "gk")
  cmin_dt <- min_wt[, .(cci_min = sum(w)), by = idx]

  out <- data.table(idx = seq_len(nrow(dt)))
  out <- merge(out, cmin_dt, by = "idx", all.x = TRUE)
  out <- merge(out, cmax_dt, by = "idx", all.x = TRUE)
  out[is.na(cci_min), cci_min := 0L]
  out[is.na(cci_max), cci_max := 0L]
  out[, cci_mid        := (cci_min + cci_max) / 2]
  out[, interval_width := cci_max - cci_min]
  setorder(out, idx)
  out[, .(cci_min        = as.integer(cci_min),
          cci_max        = as.integer(cci_max),
          cci_mid,
          interval_width = as.integer(interval_width))]
}

#' S2: Probabilistic (E-CCI) for a single patient
#'
#' @description
#' For every 3-character prefix, computes the expected probability that the
#' true 4-digit subcode maps to each Charlson group, using national Destatis
#' frequencies as the prior.  The expected CCI contribution of each group is
#' P(group | prefix) × weight(group).  Fully deterministic; no draws needed.
#'
#' @param icd_anon  Character vector of 3-character ICD-10-GM prefixes.
#' @param quan_map  Quan mapping list.
#' @param cache     Precomputed lookup environment.
#' @param age       Optional numeric age (currently reserved for future use).
#' @return Named list with element \code{e_cci} (numeric, rounded to 2 dp).
#' @export
cci_probabilistic <- function(icd_anon, quan_map, cache, age = NULL) {
  codes <- unique(icd3(icd_anon))
  gn    <- get(".group_names", envir = cache)
  dm    <- get(".dep_map",     envir = cache)
  wm    <- get(".wt_map",      envir = cache)
  gp    <- setNames(numeric(length(gn)), gn)
  for (code in codes) {
    pc <- get_prefix_cache(code, cache)
    if (is.null(pc)) next
    probs <- pc$child_freqs
    for (gk in gn) {
      matching <- vapply(probs$code_nodot,
        function(cn) gk %in% (pc$child_groups[[cn]] %||% character(0)),
        logical(1))
      gp[gk] <- max(gp[gk], sum(probs$prob[matching], na.rm = TRUE))
    }
  }
  cont <- numeric(length(gn)); names(cont) <- gn
  for (gk in gn) {
    if (gp[gk] <= 0) next
    deps <- dm[[gk]]
    if (length(deps) > 0 && any(gp[deps] > 0, na.rm = TRUE)) next
    cont[gk] <- gp[gk] * wm[gk]
  }
  list(e_cci = round(sum(cont), 2))
}

#' S2: Probabilistic (E-CCI) for a batch of patients
#' @param dt       data.table with column \code{diagnosen}.
#' @param quan_map Quan mapping list.
#' @param cache    Precomputed lookup environment.
#' @return Numeric vector of expected CCI, one value per row.
#' @export
cci_probabilistic_batch <- function(dt, quan_map, cache) {
  gn <- get(".group_names",   envir = cache)
  wm <- get(".wt_map",        envir = cache)
  dl <- get(".dep_lookup_dt", envir = cache)

  all_pref <- ls(cache)
  all_pref <- all_pref[!startsWith(all_pref, ".")]
  pg_rows  <- list()
  for (pref in all_pref) {
    pc  <- get(pref, envir = cache)
    cdt <- pc$child_freqs
    cg  <- pc$child_groups
    for (gk in gn) {
      p <- sum(cdt$prob[vapply(cdt$code_nodot,
        function(cn) gk %in% (cg[[cn]] %||% character(0)), logical(1))],
        na.rm = TRUE)
      if (p > 0)
        pg_rows[[length(pg_rows)+1L]] <- list(pref = pref, gk = gk, gp = p)
    }
  }
  pgp <- rbindlist(pg_rows)

  dx_list <- strsplit(as.character(dt$diagnosen), "\\|+")
  long <- data.table(
    idx  = rep(seq_len(nrow(dt)), lengths(dx_list)),
    code = toupper(gsub("[^A-Z0-9]", "", unlist(dx_list, use.names = FALSE)))
  )
  long <- long[nchar(code) >= 3]
  long[, pref := substr(code, 1, 3)]
  long <- unique(long[, .(idx, pref)])

  hits   <- merge(long, pgp, by = "pref", allow.cartesian = TRUE)
  enc_gp <- hits[, .(gp = max(gp)), by = .(idx, gk)]

  if (nrow(dl) > 0 && nrow(enc_gp) > 0) {
    chk <- merge(enc_gp[, .(idx, child_gk = gk)], dl,
                 by.x = "child_gk", by.y = "child", allow.cartesian = TRUE)
    chk <- chk[enc_gp, on = .(idx = idx, parent = gk), nomatch = 0L]
    if (nrow(chk) > 0) {
      sup <- unique(chk[, .(idx, gk = child_gk)])
      enc_gp <- enc_gp[!sup, on = c("idx","gk")]
    }
  }

  wt_dt <- data.table(gk = names(wm), w = unname(wm))
  sc    <- merge(enc_gp, wt_dt, by = "gk")
  sc[, contribution := gp * w]
  ecci  <- sc[, .(e_cci = round(sum(contribution), 2)), by = idx]

  out        <- rep(0, nrow(dt))
  out[ecci$idx] <- ecci$e_cci
  out
}

#' S3: Multiple Imputation CCI for a single patient
#'
#' @description
#' Treats code truncation as a Missing At Random problem.  For each 3-character
#' prefix, one 4-digit subcode is drawn from the Destatis frequency distribution
#' m times; the gold-standard CCI is computed on each imputed code set; Rubin's
#' Rules pooling (arithmetic mean) gives the final estimate.
#'
#' @param icd_anon  Character vector of 3-character ICD-10-GM prefixes.
#' @param quan_map  Quan mapping list.
#' @param cache     Precomputed lookup environment.
#' @param age       Optional numeric age (reserved).
#' @param m         Number of imputations (default 20 for production).
#' @param seed      Random seed for reproducibility.
#' @return Named list with element \code{mi_cci} (numeric).
#' @export
cci_mi <- function(icd_anon, quan_map, cache,
                   age = NULL, m = 20L, seed = 42L) {
  set.seed(seed)
  codes <- unique(icd3(icd_anon))
  if (length(codes) == 0) return(list(mi_cci = 0))
  pools <- lapply(codes, function(code) {
    pc <- get_prefix_cache(code, cache)
    if (is.null(pc)) return(data.table(code_nodot = drop_dot(code), prob = 1))
    p  <- pc$child_freqs[prob > 0]
    if (nrow(p) == 0) return(data.table(code_nodot = drop_dot(code), prob = 1))
    p
  })
  vals <- integer(m)
  for (i in seq_len(m)) {
    imp    <- vapply(pools, function(pool)
      add_dot4(pool$code_nodot[sample.int(nrow(pool), 1L, prob = pool$prob)]),
      character(1))
    vals[i] <- cci_gold(imp, quan_map)$cci
  }
  list(mi_cci = round(mean(vals), 2))
}

#' S3: Multiple Imputation CCI for a batch of patients
#' @param dt     data.table with column \code{diagnosen}.
#' @param quan_map Quan mapping list.
#' @param cache    Precomputed lookup environment.
#' @param m        Number of imputations (default 20).
#' @param seed     Random seed.
#' @return Numeric vector of MI-CCI estimates, one per row.
#' @export
cci_mi_batch <- function(dt, quan_map, cache, m = 20L, seed = 42L) {
  set.seed(seed)
  pl <- build_pattern_lookup(quan_map)
  dl <- build_dep_lookup(quan_map)
  n  <- nrow(dt)

  dx_list   <- strsplit(as.character(dt$diagnosen), "\\|+")
  code_sets <- lapply(dx_list, function(x)
    unique(substr(toupper(gsub("[^A-Z0-9]","",x)), 1, 3)))

  all_prefs  <- unique(unlist(code_sets))
  donor_pools <- list()
  for (pref in all_prefs) {
    pc <- get_prefix_cache(pref, cache)
    if (is.null(pc)) {
      donor_pools[[pref]] <- data.table(code_nodot = pref, prob = 1)
    } else {
      p <- pc$child_freqs[prob > 0]
      donor_pools[[pref]] <- if (nrow(p) > 0) p else
                             data.table(code_nodot = pref, prob = 1)
    }
  }
  message(sprintf("  S3: %d patients, %d prefixes, m=%d",
                  n, length(all_prefs), m))

  cci_sum <- numeric(n)
  for (imp in seq_len(m)) {
    imp_dx <- vapply(seq_len(n), function(i) {
      prefs <- code_sets[[i]]
      drawn <- vapply(prefs, function(pref) {
        pool <- donor_pools[[pref]]
        add_dot4(pool$code_nodot[sample.int(nrow(pool), 1L, prob = pool$prob)])
      }, character(1))
      paste(drawn, collapse = "|")
    }, character(1))
    imp_dt  <- data.table(diagnosen = imp_dx)
    cci_sum <- cci_sum + cci_gold_batch(imp_dt, quan_map, pl, dl)
  }
  round(cci_sum / m, 2)
}

#' S4: Bayesian CCI for a single patient
#'
#' @description
#' Extends multiple imputation with a Dirichlet prior.  Concentration
#' parameters are set to alpha_k = alpha_0 * P(subcode | prefix) + 0.5.
#' A Gamma representation draws posterior samples; the point estimate is the
#' posterior median.  Compared to S3, S4 offers slightly more regularization
#' for rare subcodes at the cost of a marginally wider error distribution.
#'
#' @param icd_anon  Character vector of 3-character ICD-10-GM prefixes.
#' @param quan_map  Quan mapping list.
#' @param cache     Precomputed lookup environment.
#' @param age       Optional numeric age (reserved).
#' @param n_draws   Number of posterior draws (default 25).
#' @param alpha_0   Dirichlet concentration multiplier (default 10).
#' @param seed      Random seed.
#' @return Named list with element \code{posterior_median} (numeric).
#' @export
cci_bayesian <- function(icd_anon, quan_map, cache,
                         age = NULL, n_draws = 25L,
                         alpha_0 = 10, seed = 42L) {
  set.seed(seed)
  codes <- unique(icd3(icd_anon))
  if (length(codes) == 0) return(list(posterior_median = 0))
  posts <- lapply(codes, function(code) {
    pc <- get_prefix_cache(code, cache)
    if (is.null(pc)) return(NULL)
    p  <- pc$child_freqs[prob > 0]
    if (nrow(p) == 0) return(NULL)
    list(alpha = alpha_0 * p$prob + 0.5, names = p$code_nodot)
  })
  posts <- Filter(Negate(is.null), posts)
  if (length(posts) == 0) return(list(posterior_median = 0))
  draws <- integer(n_draws)
  for (d in seq_len(n_draws)) {
    imp <- character(0)
    for (post in posts) {
      theta <- rgamma(length(post$alpha), shape = post$alpha, rate = 1)
      theta <- theta / sum(theta)
      imp   <- c(imp, add_dot4(post$names[sample.int(length(theta),1L,prob=theta)]))
    }
    draws[d] <- cci_gold(imp, quan_map)$cci
  }
  list(posterior_median = round(median(draws), 2))
}

#' S4: Bayesian CCI for a batch of patients
#' @param dt       data.table with column \code{diagnosen}.
#' @param quan_map Quan mapping list.
#' @param cache    Precomputed lookup environment.
#' @param n_draws  Number of Dirichlet draws per patient (default 25).
#' @param alpha_0  Prior concentration (default 10).
#' @param seed     Random seed.
#' @return Numeric vector of posterior medians, one per row.
#' @export
cci_bayesian_batch <- function(dt, quan_map, cache,
                               n_draws = 25L, alpha_0 = 10, seed = 42L) {
  set.seed(seed)
  pl <- build_pattern_lookup(quan_map)
  dl <- build_dep_lookup(quan_map)
  n  <- nrow(dt)

  dx_list   <- strsplit(as.character(dt$diagnosen), "\\|+")
  code_sets <- lapply(dx_list, function(x)
    unique(substr(toupper(gsub("[^A-Z0-9]","",x)), 1, 3)))

  all_prefs  <- unique(unlist(code_sets))
  dir_params <- list()
  for (pref in all_prefs) {
    pc <- get_prefix_cache(pref, cache)
    if (is.null(pc)) {
      dir_params[[pref]] <- list(alpha = 1, names = pref)
    } else {
      p <- pc$child_freqs[prob > 0]
      if (nrow(p) == 0) {
        dir_params[[pref]] <- list(alpha = 1, names = pref)
      } else {
        dir_params[[pref]] <- list(alpha = alpha_0 * p$prob + 0.5,
                                   names = p$code_nodot)
      }
    }
  }
  message(sprintf("  S4: %d patients, %d prefixes, n_draws=%d",
                  n, length(all_prefs), n_draws))

  cci_mat <- matrix(0L, nrow = n, ncol = n_draws)
  for (d in seq_len(n_draws)) {
    imp_dx <- vapply(seq_len(n), function(i) {
      prefs <- code_sets[[i]]
      drawn <- vapply(prefs, function(pref) {
        par   <- dir_params[[pref]]
        theta <- rgamma(length(par$alpha), shape = par$alpha, rate = 1)
        theta <- theta / sum(theta)
        add_dot4(par$names[sample.int(length(theta),1L,prob=theta)])
      }, character(1))
      paste(drawn, collapse = "|")
    }, character(1))
    cci_mat[, d] <- cci_gold_batch(data.table(diagnosen = imp_dx),
                                   quan_map, pl, dl)
  }
  apply(cci_mat, 1, function(x) round(median(x), 2))
}

#' S5: Train the XGBoost ML-CCI model
#'
#' @description
#' Trains an XGBoost gradient boosting regressor on gold-standard CCI labels.
#' Features are binary 3-character ICD indicator columns, patient age,
#' number of diagnoses, and length-of-stay decile.  An optional calibration
#' set triggers isotonic regression post-processing to correct systematic bias.
#'
#' @param dt_train      data.table with columns \code{diagnosen}, \code{age},
#'   \code{stay_in_days}.
#' @param quan_map      Quan mapping list.
#' @param nrounds       XGBoost iterations (default 500).
#' @param max_depth     Tree depth (default 6).
#' @param eta           Learning rate (default 0.05).
#' @param subsample     Row subsampling ratio (default 0.8).
#' @param colsample_bytree Column subsampling ratio (default 0.8).
#' @param dt_calibrate  Optional data.table for isotonic calibration.
#' @param seed          Random seed.
#' @return A list with elements \code{model}, \code{vocab},
#'   \code{calibration_fn}, \code{pattern_lookup}, \code{dep_lookup}.
#' @export
cci_ml_train <- function(dt_train, quan_map,
                         nrounds = 500L, max_depth = 6L, eta = 0.05,
                         subsample = 0.8, colsample_bytree = 0.8,
                         dt_calibrate = NULL, seed = 42L) {
  if (!requireNamespace("xgboost", quietly = TRUE))
    stop("Package 'xgboost' is required. Install with: install.packages('xgboost')")
  set.seed(seed)
  pl   <- build_pattern_lookup(quan_map)
  dl   <- build_dep_lookup(quan_map)
  feat <- .prepare_ml_features(dt_train, quan_map,
                                pattern_lookup = pl, dep_lookup = dl)
  dtr  <- xgboost::xgb.DMatrix(data = feat$X, label = feat$y)
  params <- list(
    objective         = "reg:squarederror",
    eval_metric       = "mae",
    max_depth         = max_depth,
    eta               = eta,
    subsample         = subsample,
    colsample_bytree  = colsample_bytree,
    nthread           = max(1L, parallel::detectCores() - 1L)
  )
  model  <- xgboost::xgb.train(params = params, data = dtr,
                                nrounds = nrounds, verbose = 0)
  cal_fn <- identity
  if (!is.null(dt_calibrate) && nrow(dt_calibrate) > 0) {
    fc   <- .prepare_ml_features(dt_calibrate, quan_map, vocab = feat$vocab,
                                  pattern_lookup = pl, dep_lookup = dl)
    raw  <- predict(model, xgboost::xgb.DMatrix(data = fc$X))
    iso  <- isoreg(raw, fc$y)
    cal_fn <- approxfun(iso$x[iso$ord], iso$yf, rule = 2, ties = mean)
  }
  list(model = model, vocab = feat$vocab, calibration_fn = cal_fn,
       pattern_lookup = pl, dep_lookup = dl)
}

#' S5: Predict ML-CCI for new patients
#' @param ml_obj  Object returned by \code{cci_ml_train()}.
#' @param dt_new  data.table with columns \code{diagnosen}, \code{age},
#'   \code{stay_in_days}.
#' @param quan_map Quan mapping list.
#' @return Numeric vector of calibrated CCI estimates (clamped to [0, 29]).
#' @export
cci_ml_predict <- function(ml_obj, dt_new, quan_map) {
  if (!requireNamespace("xgboost", quietly = TRUE))
    stop("Package 'xgboost' is required.")
  feat <- .prepare_ml_features(dt_new, quan_map,
                                vocab         = ml_obj$vocab,
                                pattern_lookup = ml_obj$pattern_lookup,
                                dep_lookup    = ml_obj$dep_lookup)
  raw  <- predict(ml_obj$model,
                  xgboost::xgb.DMatrix(data = feat$X))
  pred <- ml_obj$calibration_fn(raw)
  pmin(pmax(round(pred, 2), 0), 29)
}

.prepare_ml_features <- function(dt, quan_map, vocab = NULL,
                                  pattern_lookup = NULL, dep_lookup = NULL) {
  if (is.null(pattern_lookup)) pattern_lookup <- build_pattern_lookup(quan_map)
  if (is.null(dep_lookup))     dep_lookup     <- build_dep_lookup(quan_map)
  y <- cci_gold_batch(dt, quan_map, pattern_lookup, dep_lookup)

  code_sets <- strsplit(as.character(dt$diagnosen), "\\|+")
  code_sets <- lapply(code_sets, function(x)
    unique(substr(toupper(gsub("[^A-Z0-9]","",x)), 1, 3)))

  if (is.null(vocab)) {
    freq  <- sort(table(unlist(code_sets)), decreasing = TRUE)
    mfreq <- max(2L, ceiling(nrow(dt) * 0.01))
    vocab <- sort(names(freq[freq >= mfreq]))
    if (length(vocab) == 0) vocab <- sort(names(freq[freq >= 1]))
  }

  n        <- nrow(dt)
  code_mat <- matrix(0L, nrow = n, ncol = length(vocab))
  colnames(code_mat) <- paste0("icd_", vocab)
  for (i in seq_len(n)) {
    m2 <- code_sets[[i]] %in% vocab
    if (any(m2)) code_mat[i, match(code_sets[[i]][m2], vocab)] <- 1L
  }

  los_r  <- frank(dt$stay_in_days, ties.method = "dense")
  los_d  <- as.integer(ceiling(10 * los_r / max(los_r, na.rm = TRUE)))
  los_d[los_d > 10L] <- 10L

  X <- cbind(code_mat,
             age         = as.numeric(dt$age),
             n_diagnoses = lengths(code_sets),
             los_decile  = los_d)
  list(X = X, y = y, vocab = vocab)
}

#' S6: Train NNLS ensemble weights on a calibration set
#'
#' @description
#' IMPORTANT: S6 is a meta-learner.  Its inputs are the predictions of
#' S1 through S5 — not raw ICD codes.  You must therefore run S1, S2, S3,
#' S4, and S5 on a calibration set first, collect their outputs into a
#' five-column matrix, and then call this function.
#'
#' Non-Negative Least Squares ensures every weight is >= 0.  Weights are
#' normalized so they sum to 1.  The ensemble is then applied to new
#' patients by calling \code{cci_ensemble()}.
#'
#' @param estimates_matrix Numeric matrix with exactly five columns:
#'   S1_mid, S2_ecci, S3_mi, S4_bayes, S5_ml (one row per patient).
#' @param cci_gold_vec Numeric or integer vector of gold-standard CCI
#'   values for the calibration patients.
#' @return A list with elements \code{weights} (named numeric vector of
#'   length 5) and \code{residual_norm} (NNLS residual).
#' @export
cci_ensemble_train <- function(estimates_matrix, cci_gold_vec) {
  if (!requireNamespace("nnls", quietly = TRUE))
    stop("Package 'nnls' is required. Install with: install.packages('nnls')")
  stopifnot(ncol(estimates_matrix) == 5,
            length(cci_gold_vec) == nrow(estimates_matrix))
  colnames(estimates_matrix) <- c("S1_mid","S2_ecci","S3_mi","S4_bayes","S5_ml")
  fit  <- nnls::nnls(estimates_matrix, cci_gold_vec)
  w    <- fit$x; names(w) <- colnames(estimates_matrix)
  wtot <- sum(w)
  if (wtot > 0) w <- w / wtot
  list(weights = round(w, 4), residual_norm = fit$deviance)
}

#' S6: Apply trained ensemble weights to new predictions
#'
#' @description
#' After training with \code{cci_ensemble_train()}, use this function to
#' produce final ensemble estimates.  The input is the same five-column
#' matrix of S1-S5 outputs for the new patients.
#'
#' @param estimates  Either a five-element numeric vector (single patient)
#'   or a matrix with five columns (batch).  Columns must be in order:
#'   S1_mid, S2_ecci, S3_mi, S4_bayes, S5_ml.
#' @param weights    Named numeric vector of length 5, from
#'   \code{cci_ensemble_train()}.
#' @return Numeric vector of ensemble CCI estimates.
#' @export
cci_ensemble <- function(estimates, weights) {
  if (is.vector(estimates) && length(estimates) == 5)
    estimates <- matrix(estimates, nrow = 1)
  as.numeric(estimates %*% weights)
}

# ── Plot generator: ggplot2 namespace bug fixed ───────────────────────────────
# Root cause: `gg <- ggplot2` fails when ggplot2 is only in Imports (not
# attached). Fix: asNamespace("ggplot2") always returns the namespace
# regardless of whether the package has been explicitly attached.

.generate_plots <- function(val, strategies, sn, output_dir) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 is required for plots. Run: install.packages('ggplot2')")

  # FIXED: was `gg <- ggplot2` (only works if ggplot2 is attached via library())
  gg <- asNamespace("ggplot2")

  pal <- c(
    "Gold Standard"    = "#1B1B1B",
    "S1 Interval"      = "#E69F00",
    "S2 Probabilistic" = "#56B4E9",
    "S3 MI-CCI"        = "#009E73",
    "S4 Bayesian"      = "#CC79A7",
    "S5 ML-CCI"        = "#0072B2",
    "S6 Ensemble"      = "#D55E00"
  )

  # Figure 1: MAE bar chart
  mae_v <- vapply(strategies, function(s)
    mean(abs(val[[s]] - val$cci_gold), na.rm = TRUE), numeric(1))
  p1d <- data.table(Strategy = factor(sn, levels = sn[order(mae_v)]),
                    MAE = mae_v)
  p1 <- gg$ggplot(p1d, gg$aes(x = Strategy, y = MAE)) +
    gg$geom_col(fill = "#2C3E50", width = 0.6) +
    gg$geom_text(gg$aes(label = sprintf("%.3f", MAE)),
                 vjust = -0.3, size = 3.5) +
    gg$labs(title = "Mean Absolute Error by Strategy (Validation Cohort)",
            y = "MAE", x = NULL) +
    gg$theme_minimal(base_size = 12) +
    gg$theme(axis.text.x = gg$element_text(angle = 30, hjust = 1))
  gg$ggsave(file.path(output_dir, "fig_01_mae_comparison.pdf"),
            p1, width = 8, height = 5)

  # Figure 2: Predicted vs Gold (hexbin faceted)
  ld <- melt(val[, c("cci_gold", strategies), with = FALSE],
             id.vars = "cci_gold",
             variable.name = "strategy", value.name = "predicted")
  ld[, strategy_label := sn[match(strategy, strategies)]]
  p2 <- gg$ggplot(ld, gg$aes(x = cci_gold, y = predicted)) +
    gg$geom_abline(slope = 1, intercept = 0,
                   linetype = "dashed", color = "red", linewidth = 0.5) +
    gg$geom_hex(bins = 30) +
    gg$scale_fill_viridis_c(trans = "log10") +
    gg$facet_wrap(~strategy_label, ncol = 3) +
    gg$labs(title = "Predicted vs. Gold-Standard CCI",
            x = "Gold CCI", y = "Predicted CCI") +
    gg$theme_minimal(base_size = 11) +
    gg$coord_equal(xlim = c(0, 15), ylim = c(0, 15))
  gg$ggsave(file.path(output_dir, "fig_02_predicted_vs_gold.pdf"),
            p2, width = 12, height = 8)

  # Figure 3: Temporal stability
  ym <- val[, lapply(.SD, function(x) mean(abs(x - cci_gold), na.rm = TRUE)),
            by = year, .SDcols = strategies]
  yl <- melt(ym, id.vars = "year",
             variable.name = "strategy", value.name = "MAE")
  yl[, strategy_label := sn[match(strategy, strategies)]]
  p3 <- gg$ggplot(yl, gg$aes(x = year, y = MAE, color = strategy_label)) +
    gg$geom_line(linewidth = 0.8) + gg$geom_point(size = 2) +
    gg$scale_color_manual(values = pal[-1]) +
    gg$labs(title = "Temporal Stability: MAE by Year",
            x = "Year", y = "MAE", color = "Strategy") +
    gg$theme_minimal(base_size = 12) +
    gg$theme(legend.position = "bottom")
  gg$ggsave(file.path(output_dir, "fig_03_temporal_stability.pdf"),
            p3, width = 10, height = 6)

  # Figure 4: CCI density overlay
  dens_dt <- data.table(value = val$cci_gold, Source = "Gold Standard")
  for (j in seq_along(strategies))
    dens_dt <- rbind(dens_dt,
      data.table(value = val[[strategies[j]]], Source = sn[j]))
  dens_dt[, Source := factor(Source, levels = c("Gold Standard", sn))]
  p4 <- gg$ggplot(dens_dt, gg$aes(x = value, fill = Source, color = Source)) +
    gg$geom_density(alpha = 0.15, linewidth = 0.7, adjust = 1.5) +
    gg$scale_fill_manual(values  = pal) +
    gg$scale_color_manual(values = pal) +
    gg$xlim(-0.5, 15) +
    gg$labs(title = "CCI Distribution: Gold Standard vs. All Strategies",
            x = "CCI", y = "Density") +
    gg$theme_minimal(base_size = 12) +
    gg$theme(legend.position = "bottom",
             legend.title    = gg$element_blank())
  gg$ggsave(file.path(output_dir, "fig_04_density_overlay.pdf"),
            p4, width = 10, height = 6)

  # Figure 5: Bland-Altman
  ba_rows <- lapply(seq_along(strategies), function(j) {
    data.table(
      mean_val = (val$cci_gold + val[[strategies[j]]]) / 2,
      diff_val = val[[strategies[j]]] - val$cci_gold,
      Strategy = sn[j]
    )
  })
  ba_dt <- rbindlist(ba_rows)
  ba_dt[, Strategy := factor(Strategy, levels = sn)]
  ba_stats <- ba_dt[, .(bias    = mean(diff_val, na.rm = TRUE),
                         sd_diff = sd(diff_val,   na.rm = TRUE)), by = Strategy]
  ba_stats[, loa_lo := bias - 1.96 * sd_diff]
  ba_stats[, loa_hi := bias + 1.96 * sd_diff]
  p5 <- gg$ggplot(ba_dt, gg$aes(x = mean_val, y = diff_val)) +
    gg$geom_hline(yintercept = 0, linetype = "solid", color = "grey50") +
    gg$geom_point(alpha = 0.3, size = 0.8, color = "#2C3E50") +
    gg$geom_hline(data = ba_stats, gg$aes(yintercept = bias),
                  linetype = "dashed", color = "#D55E00", linewidth = 0.6) +
    gg$geom_hline(data = ba_stats, gg$aes(yintercept = loa_lo),
                  linetype = "dotted", color = "#0072B2", linewidth = 0.5) +
    gg$geom_hline(data = ba_stats, gg$aes(yintercept = loa_hi),
                  linetype = "dotted", color = "#0072B2", linewidth = 0.5) +
    gg$facet_wrap(~Strategy, ncol = 3, scales = "free_y") +
    gg$labs(title    = "Bland-Altman Agreement: Predicted vs. Gold CCI",
            subtitle = "Dashed = mean bias     Dotted = 95% Limits of Agreement",
            x = "Mean of Gold and Predicted CCI",
            y = "Difference (Predicted minus Gold)") +
    gg$theme_minimal(base_size = 11)
  gg$ggsave(file.path(output_dir, "fig_05_bland_altman.pdf"),
            p5, width = 12, height = 8)

  # Table 1: Descriptive statistics
  desc_fn <- function(x, nm) {
    data.table(
      Source = nm,
      N      = sum(!is.na(x)),
      Min    = round(min(x,            na.rm = TRUE), 2),
      Q1     = round(quantile(x, 0.25, na.rm = TRUE), 2),
      Median = round(median(x,          na.rm = TRUE), 2),
      Mean   = round(mean(x,            na.rm = TRUE), 2),
      Q3     = round(quantile(x, 0.75, na.rm = TRUE), 2),
      Max    = round(max(x,            na.rm = TRUE), 2),
      SD     = round(sd(x,             na.rm = TRUE), 2)
    )
  }
  tab1 <- rbind(
    desc_fn(val$cci_gold, "Gold Standard"),
    rbindlist(lapply(seq_along(strategies), function(j)
      desc_fn(val[[strategies[j]]], sn[j])))
  )
  fwrite(tab1, file.path(output_dir, "tab_01_descriptive_statistics.csv"))

  # Table 2: Agreement metrics
  tab2 <- rbindlist(lapply(seq_along(strategies), function(j) {
    d  <- val[[strategies[j]]] - val$cci_gold
    data.table(
      Strategy     = sn[j],
      MAE          = round(mean(abs(d), na.rm = TRUE), 4),
      RMSE         = round(sqrt(mean(d^2, na.rm = TRUE)), 4),
      R_squared    = round(1 - sum(d^2, na.rm=TRUE) /
                      sum((val$cci_gold - mean(val$cci_gold,na.rm=TRUE))^2,
                           na.rm=TRUE), 4),
      Bias         = round(mean(d, na.rm = TRUE), 4),
      SD_diff      = round(sd(d, na.rm = TRUE), 4),
      LoA_lower    = round(mean(d,na.rm=TRUE) - 1.96*sd(d,na.rm=TRUE), 4),
      LoA_upper    = round(mean(d,na.rm=TRUE) + 1.96*sd(d,na.rm=TRUE), 4),
      Pct_exact    = round(100 * mean(abs(d) < 0.5, na.rm=TRUE), 1),
      Pct_within_1 = round(100 * mean(abs(d) <= 1,  na.rm=TRUE), 1)
    )
  }))
  fwrite(tab2, file.path(output_dir, "tab_02_agreement_metrics.csv"))
  cat("  5 figures + 2 tables saved.\n")
}

# ── Full pipeline (requires Parquet data from UMM §21 extract) ────────────────

#' Run the complete miCCI validation pipeline
#'
#' @description
#' Loads raw Parquet data, applies temporal splits (2010-2018 train /
#' 2019 calibrate / 2020-2024 validate), runs all six strategies, prints
#' metrics, saves CSVs, and generates five publication-quality PDF figures.
#'
#' @param data_path   Path to the Parquet file.
#' @param output_dir  Destination folder for all outputs.
#' @param sample_size Optional integer for quick testing (subsamples each split).
#' @return Invisibly, a list with \code{splits}, \code{ml_obj}, \code{ensemble}.
#' @export
run_pipeline <- function(data_path, output_dir, sample_size = NULL) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  cat(paste0(
    strrep("=", 60), "\n",
    " miCCI Validation Pipeline  v0.4.3\n",
    " Data:   ", data_path, "\n",
    " Output: ", output_dir, "\n",
    strrep("=", 60), "\n"
  ))

  cat("[1/8] Loading Parquet data\n")
  dt     <- .load_mannheim_data(data_path)
  splits <- .temporal_split(dt)
  cat(sprintf("  Train=%d  Cal=%d  Val=%d\n",
              nrow(splits$train), nrow(splits$calibrate),
              nrow(splits$validate)))
  if (!is.null(sample_size)) {
    set.seed(42)
    for (nm in names(splits))
      if (nrow(splits[[nm]]) > sample_size)
        splits[[nm]] <- splits[[nm]][sample(.N, sample_size)]
    cat(sprintf("  Subsampled to %d per split\n", sample_size))
  }

  cat("[2/8] Precomputing lookups\n")
  quan_map    <- load_quan_map()
  dt_dest     <- .timed("Destatis", load_destatis())
  cache       <- .timed("Cache",    precompute_lookups(dt_dest, quan_map))
  pl          <- build_pattern_lookup(quan_map)
  dl          <- build_dep_lookup(quan_map)

  cat("[3/8] Gold CCI\n")
  for (nm in names(splits))
    splits[[nm]][, cci_gold :=
      .timed(nm, cci_gold_batch(.SD, quan_map, pl, dl))]
  val <- splits$validate
  cat(sprintf("  Val: mean=%.1f  median=%d  max=%d\n",
              mean(val$cci_gold), as.integer(median(val$cci_gold)),
              max(val$cci_gold)))

  cat("[4/8] S1: Interval\n")
  for (nm in names(splits)) {
    res <- .timed(nm, cci_interval_batch(splits[[nm]], quan_map, cache))
    splits[[nm]][, c("s1_min","s1_max","s1_mid","s1_width") :=
                   .(res$cci_min, res$cci_max, res$cci_mid, res$interval_width)]
  }

  cat("[5/8] S2: Probabilistic\n")
  for (nm in names(splits))
    splits[[nm]][, s2_ecci :=
      .timed(nm, cci_probabilistic_batch(.SD, quan_map, cache))]

  cat("[6/8] S3: MI-CCI\n")
  for (nm in names(splits))
    splits[[nm]][, s3_mi :=
      .timed(nm, cci_mi_batch(.SD, quan_map, cache, m = 20L))]

  cat("[6.5/8] S4: Bayesian\n")
  for (nm in names(splits))
    splits[[nm]][, s4_bayes :=
      .timed(nm, cci_bayesian_batch(.SD, quan_map, cache, n_draws = 25L))]

  cat("[7/8] S5: ML\n")
  ml_obj <- .timed("train",
    cci_ml_train(splits$train, quan_map, dt_calibrate = splits$calibrate))
  for (nm in names(splits))
    splits[[nm]][, s5_ml :=
      .timed(nm, cci_ml_predict(ml_obj, .SD, quan_map))]

  cat("[8/8] S6: Ensemble\n")
  cal     <- splits$calibrate
  est_mat <- as.matrix(cal[, .(s1_mid, s2_ecci, s3_mi, s4_bayes, s5_ml)])
  ens     <- cci_ensemble_train(est_mat, cal$cci_gold)
  cat(sprintf("  Weights: %s\n",
    paste(sprintf("%s=%.3f", names(ens$weights), ens$weights),
          collapse = " | ")))
  for (nm in names(splits)) {
    mat <- as.matrix(
      splits[[nm]][, .(s1_mid, s2_ecci, s3_mi, s4_bayes, s5_ml)])
    splits[[nm]][, s6_ensemble := cci_ensemble(mat, ens$weights)]
  }

  strategies <- c("s1_mid","s2_ecci","s3_mi","s4_bayes","s5_ml","s6_ensemble")
  sn <- c("S1 Interval","S2 Probabilistic","S3 MI-CCI",
          "S4 Bayesian","S5 ML-CCI","S6 Ensemble")

  cat("\n=== RESULTS ===\n")
  rl <- list()
  for (nm in names(splits)) {
    cat(sprintf("\n--- %s (n=%d) ---\n", toupper(nm), nrow(splits[[nm]])))
    for (j in seq_along(strategies)) {
      m <- eval_metrics(splits[[nm]][[strategies[j]]], splits[[nm]]$cci_gold)
      cat(sprintf("  %-20s MAE=%.3f  RMSE=%.3f  R\u00b2=%.3f\n",
                  sn[j], m["mae"], m["rmse"], m["r_squared"]))
      rl[[paste(nm, strategies[j], sep = "_")]] <- m
    }
    cov <- mean(splits[[nm]]$cci_gold >= splits[[nm]]$s1_min &
                splits[[nm]]$cci_gold <= splits[[nm]]$s1_max, na.rm = TRUE)
    cat(sprintf("  S1 Coverage: %.1f%%\n", 100 * cov))
  }

  cat("\n=== SAVING ===\n")
  out_cols <- c("falnr","year","age","cci_gold",
                "s1_min","s1_max","s1_mid",
                "s2_ecci","s3_mi","s4_bayes","s5_ml","s6_ensemble")
  fwrite(splits$validate[, out_cols, with = FALSE],
         file.path(output_dir, "cci_predictions_validation.csv"))
  mdt <- rbindlist(lapply(names(rl), function(nm) {
    p <- strsplit(nm, "_", fixed = FALSE)[[1]]
    data.table(split = p[1], strategy = paste(p[-1], collapse="_"), t(rl[[nm]]))
  }))
  fwrite(mdt,   file.path(output_dir, "cci_metrics_summary.csv"))
  fwrite(data.table(strategy = names(ens$weights), weight = ens$weights),
         file.path(output_dir, "ensemble_weights.csv"))
  .generate_plots(splits$validate, strategies, sn, output_dir)
  cat(sprintf("  All outputs written to: %s\n", output_dir))
  invisible(list(splits = splits, ml_obj = ml_obj, ensemble = ens))
}

# ── run_demo(): works entirely on built-in data — no external files needed ────

#' Quick demonstration using the built-in dummy dataset
#'
#' @description
#' Loads the \code{micci_patients} dataset shipped with the package, simulates
#' §21 KHEntgG anonymization (truncation to 3-character codes), runs all six
#' strategies, and prints metrics.  No Parquet file or Destatis download is
#' required.  Use this function first after installation to verify everything
#' works.
#'
#' @param m       MI imputations for S3 (default 5 for speed; use 20 for
#'   production).
#' @param n_draws Dirichlet draws for S4 (default 5; use 25 for production).
#' @return Invisibly, the annotated data.table with all strategy scores.
#' @examples
#' result <- run_demo()
#' print(result[, .(p_id, cci_gold, s1_mid, s2_ecci, s3_mi)])
#' @export
run_demo <- function(m = 5L, n_draws = 5L) {
  cat(paste0(
    strrep("=", 55), "\n",
    " miCCI v0.4.3  --  Quick Demo\n",
    " Built-in dataset: micci_patients (", nrow(micci_patients),
    " patients)\n",
    strrep("=", 55), "\n\n"
  ))

  dt <- as.data.table(micci_patients)

  # Prepare batch-compatible columns
  dt[, diagnosen    := list_full_icd_code]
  dt[, stay_in_days := los_in_days]
  dt[, age          := 60]   # placeholder when age is not in the dataset
  dt[, n_diagnoses  := lengths(strsplit(diagnosen, "\\|+"))]
  los_r <- frank(dt$stay_in_days, ties.method = "dense")
  dt[, los_decile := as.integer(ceiling(10 * los_r / max(los_r, 1)))]

  # Precompute lookups once
  cat("Step 1/8  Loading Quan mapping and Destatis frequencies...\n")
  quan_map <- load_quan_map()
  dt_dest  <- load_destatis()
  cache    <- precompute_lookups(dt_dest, quan_map)
  pl       <- build_pattern_lookup(quan_map)
  dl       <- build_dep_lookup(quan_map)

  # Gold CCI from full 4-digit codes
  cat("Step 2/8  Computing gold-standard CCI from full codes...\n")
  dt[, cci_gold := cci_gold_batch(dt, quan_map, pl, dl)]
  cat(sprintf("  mean=%.2f  median=%d  max=%d\n",
              mean(dt$cci_gold), as.integer(median(dt$cci_gold)),
              max(dt$cci_gold)))

  # Simulate §21 KHEntgG: truncate to 3-character prefixes
  dt[, diagnosen_anon := sapply(list_full_icd_code, function(x) {
    codes <- split_icd(x)
    paste(unique(substr(drop_dot(normalize_icd(codes)), 1, 3)), collapse = "|")
  })]
  dt_anon <- copy(dt)
  dt_anon[, diagnosen := diagnosen_anon]

  # S1
  cat("Step 3/8  S1: Interval...\n")
  res_s1 <- cci_interval_batch(dt_anon, quan_map, cache)
  dt[, c("s1_min","s1_max","s1_mid","s1_width") :=
       .(res_s1$cci_min, res_s1$cci_max, res_s1$cci_mid, res_s1$interval_width)]

  # S2
  cat("Step 4/8  S2: Probabilistic...\n")
  dt[, s2_ecci := cci_probabilistic_batch(dt_anon, quan_map, cache)]

  # S3
  cat(sprintf("Step 5/8  S3: Multiple Imputation (m=%d)...\n", m))
  dt[, s3_mi := cci_mi_batch(dt_anon, quan_map, cache, m = m)]

  # S4
  cat(sprintf("Step 6/8  S4: Bayesian (n_draws=%d)...\n", n_draws))
  dt[, s4_bayes := cci_bayesian_batch(dt_anon, quan_map, cache,
                                      n_draws = n_draws)]

  # S5 (trains on full demo set — illustrative only, not a proper train/test split)
  cat("Step 7/8  S5: ML-CCI (XGBoost, training on demo set)...\n")
  ml <- cci_ml_train(dt, quan_map)
  dt[, s5_ml := cci_ml_predict(ml, dt_anon, quan_map)]

  # S6
  cat("Step 8/8  S6: NNLS Ensemble...\n")
  mat <- as.matrix(dt[, .(s1_mid, s2_ecci, s3_mi, s4_bayes, s5_ml)])
  ens <- cci_ensemble_train(mat, dt$cci_gold)
  dt[, s6_ensemble := cci_ensemble(mat, ens$weights)]

  # Summary
  strategies <- c("s1_mid","s2_ecci","s3_mi","s4_bayes","s5_ml","s6_ensemble")
  sn <- c("S1 Interval","S2 Probabilistic","S3 MI-CCI",
          "S4 Bayesian","S5 ML-CCI","S6 Ensemble")
  cat("\n--- Demo Results ---\n")
  for (j in seq_along(strategies)) {
    mv <- eval_metrics(dt[[strategies[j]]], dt$cci_gold)
    cat(sprintf("  %-20s MAE=%.3f  RMSE=%.3f  R\u00b2=%.3f\n",
                sn[j], mv["mae"], mv["rmse"], mv["r_squared"]))
  }
  cov <- mean(dt$cci_gold >= dt$s1_min & dt$cci_gold <= dt$s1_max)
  cat(sprintf("  S1 interval coverage: %.1f%%\n\n", 100 * cov))
  cat("Done. Use run_pipeline() for the full Parquet-based validation.\n")
  invisible(dt)
}
