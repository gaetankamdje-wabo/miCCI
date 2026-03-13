# =============================================================================
# miCCI — 50_cci_ml.R
# S5: ML-CCI — XGBoost (already vectorised natively)
# =============================================================================

.prepare_ml_features <- function(dt, quan_map, vocab = NULL,
                                  pattern_lookup = NULL, dep_lookup = NULL) {
  if (is.null(pattern_lookup)) pattern_lookup <- build_pattern_lookup(quan_map)
  if (is.null(dep_lookup)) dep_lookup <- build_dep_lookup(quan_map)
  y <- cci_gold_batch(dt, quan_map, pattern_lookup, dep_lookup)

  code_sets <- strsplit(as.character(dt$diagnosen), "\\|+")
  code_sets <- lapply(code_sets, function(x) {
    unique(substr(toupper(gsub("[^A-Z0-9]", "", x)), 1, 3))
  })

  if (is.null(vocab)) {
    freq <- sort(table(unlist(code_sets)), decreasing = TRUE)
    min_freq <- max(2L, ceiling(nrow(dt) * 0.01))
    vocab <- sort(names(freq[freq >= min_freq]))
    if (length(vocab) == 0) vocab <- sort(names(freq[freq >= 1]))
  }

  n <- nrow(dt)
  code_mat <- matrix(0L, nrow = n, ncol = length(vocab))
  colnames(code_mat) <- paste0("icd_", vocab)
  for (i in seq_len(n)) {
    m <- code_sets[[i]] %in% vocab
    if (any(m)) code_mat[i, match(code_sets[[i]][m], vocab)] <- 1L
  }

  n_diag <- lengths(code_sets)
  los_r <- frank(dt$stay_in_days, ties.method = "dense")
  los_dec <- as.integer(ceiling(10 * los_r / max(los_r, na.rm = TRUE)))
  los_dec[los_dec > 10L] <- 10L

  X <- cbind(code_mat, age = as.numeric(dt$age), n_diagnoses = n_diag, los_decile = los_dec)
  list(X = X, y = y, vocab = vocab)
}

#' @export
cci_ml_train <- function(dt_train, quan_map, nrounds = 500L, max_depth = 6L, eta = 0.05,
                         subsample = 0.8, colsample_bytree = 0.8,
                         dt_calibrate = NULL, seed = 42L) {
  requireNamespace("xgboost", quietly = TRUE)
  set.seed(seed)
  pl <- build_pattern_lookup(quan_map); dl <- build_dep_lookup(quan_map)
  feat <- .prepare_ml_features(dt_train, quan_map, pattern_lookup = pl, dep_lookup = dl)
  dtrain <- xgboost::xgb.DMatrix(data = feat$X, label = feat$y)
  params <- list(objective = "reg:squarederror", eval_metric = "mae",
                 max_depth = max_depth, eta = eta, subsample = subsample,
                 colsample_bytree = colsample_bytree,
                 nthread = max(1L, parallel::detectCores() - 1L))
  model <- xgboost::xgb.train(params = params, data = dtrain, nrounds = nrounds, verbose = 0)
  cal_fn <- identity
  if (!is.null(dt_calibrate) && nrow(dt_calibrate) > 0) {
    feat_cal <- .prepare_ml_features(dt_calibrate, quan_map, vocab = feat$vocab,
                                      pattern_lookup = pl, dep_lookup = dl)
    pred_raw <- predict(model, xgboost::xgb.DMatrix(data = feat_cal$X))
    iso <- isoreg(pred_raw, feat_cal$y)
    cal_fn <- approxfun(iso$x[iso$ord], iso$yf, rule = 2, ties = mean)
  }
  list(model = model, vocab = feat$vocab, calibration_fn = cal_fn,
       pattern_lookup = pl, dep_lookup = dl)
}

#' @export
cci_ml_predict <- function(ml_obj, dt_new, quan_map) {
  requireNamespace("xgboost", quietly = TRUE)
  feat <- .prepare_ml_features(dt_new, quan_map, vocab = ml_obj$vocab,
                                pattern_lookup = ml_obj$pattern_lookup,
                                dep_lookup = ml_obj$dep_lookup)
  pred <- ml_obj$calibration_fn(predict(ml_obj$model, xgboost::xgb.DMatrix(data = feat$X)))
  pmin(pmax(round(pred, 2), 0), 29)
}
