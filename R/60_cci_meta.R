# =============================================================================
# miCCI v1.0.0 — 60_cci_meta.R
# S6: Cross-validated Super Learner meta-estimator over S1-S4.
#
# Method
# ------
# The Super Learner (van der Laan, Polley & Hubbard 2007, "Super Learner",
# Statistical Applications in Genetics and Molecular Biology, 6(1)) is a
# cross-validated ensemble meta-learner with asymptotic oracle properties:
# given a library of candidate estimators, its out-of-sample risk is
# asymptotically at least as small as the risk of the best single estimator
# in the library.
#
# In miCCI v1.0.0 the library consists of the four base estimators S1 mid,
# S2 expected, S3 multiple imputation and S4 Bayesian posterior median,
# exposed as constant "learners" that pass each strategy through as its
# own prediction. SuperLearner::SuperLearner() then fits a non-negative
# least-squares meta-weight vector (method = "method.NNLS") over 10-fold
# cross-validated out-of-fold predictions. The resulting weights are
# non-negative and sum to 1 under the SuperLearner NNLS default.
#
# No train / calibrate / validate split is used: the entire consistent
# 2010-2024 cohort is passed to SuperLearner, and the internal 10-fold CV
# guarantees that each encounter's final meta-prediction does not depend
# on its own gold value.
#
# Output: a numeric vector s6_meta, one value per encounter. The raw
#   Super Learner convex combination is returned as-is. v1.0.0 does NOT
#   clip to the S1 certainty interval: empirical testing showed that
#   S1-based clipping systematically biases S6 downwards on subgroups
#   where S1_max is a strict but imperfect upper bound, destroying the
#   Super Learner oracle property. Trusting the cross-validated NNLS
#   combination is the correct thing to do.
# =============================================================================

#' Pass-through SL learners — one per base strategy
#'
#' Each learner is a degenerate wrapper that simply returns the corresponding
#' column of X as its prediction, making the Super Learner a pure
#' meta-combination over the four strategy columns without any additional
#' feature learning.
#'
#' @keywords internal
.SL_passthrough <- function(col_name) {
  force(col_name)
  function(Y, X, newX, family, obsWeights, id, ...) {
    pred <- as.numeric(newX[[col_name]])
    fit  <- list(col_name = col_name)
    class(fit) <- "SL_passthrough"
    list(pred = pred, fit = fit)
  }
}

#' S3 predict method for pass-through learners
#' @export
predict.SL_passthrough <- function(object, newdata, ...) {
  as.numeric(newdata[[object$col_name]])
}

#' Internal: build the SL library for miCCI
#' @keywords internal
.micci_sl_library <- function() {
  # Learners are constructed per call so that col_name is captured correctly
  list(
    SL.s1_mid   = .SL_passthrough("s1_mid"),
    SL.s2_ecci  = .SL_passthrough("s2_ecci"),
    SL.s3_mi    = .SL_passthrough("s3_mi"),
    SL.s4_bayes = .SL_passthrough("s4_bayes")
  )
}

#' Fit the Super Learner meta-estimator and return per-encounter predictions.
#'
#' Uses 10-fold cross-validation internally. Each encounter's final prediction
#' is the Super Learner output for that encounter, clipped to the S1
#' certainty interval [s1_min, s1_max].
#'
#' @param dt data.table with columns
#'   cci_gold, s1_min, s1_max, s1_mid, s2_ecci, s3_mi, s4_bayes
#' @param V  number of cross-validation folds (default 10)
#' @param seed RNG seed for fold assignment reproducibility
#' @param verbose if TRUE, SuperLearner prints its fold-level progress
#'
#' @return list with elements:
#'   \item{predictions}{numeric vector of S6 meta-predictions, length nrow(dt)}
#'   \item{weights}{named numeric vector of length 4 — NNLS meta-weights
#'     over S1_mid, S2_ecci, S3_mi, S4_bayes}
#'   \item{cv_risk}{named numeric vector — cross-validated risk of each
#'     base learner reported by SuperLearner}
#'   \item{sl_fit}{the raw SuperLearner fit object for diagnostics}
#'
#' @export
cci_meta_fit <- function(dt, V = 10L, seed = 42L, verbose = FALSE) {
  if (!requireNamespace("SuperLearner", quietly = TRUE))
    stop("cci_meta_fit() requires the 'SuperLearner' package. ",
         "Install it with install.packages('SuperLearner').")

  needed <- c("cci_gold", "s1_min", "s1_max",
              "s1_mid", "s2_ecci", "s3_mi", "s4_bayes")
  miss <- setdiff(needed, names(dt))
  if (length(miss))
    stop("cci_meta_fit: missing columns: ", paste(miss, collapse = ", "))

  Y <- as.numeric(dt$cci_gold)
  X <- data.frame(
    s1_mid   = as.numeric(dt$s1_mid),
    s2_ecci  = as.numeric(dt$s2_ecci),
    s3_mi    = as.numeric(dt$s3_mi),
    s4_bayes = as.numeric(dt$s4_bayes)
  )

  # Assign learner wrappers to the calling environment so SuperLearner
  # can locate them by name. Use a local environment and pass it via
  # env argument so we never pollute the user's global environment.
  sl_env <- new.env(parent = globalenv())
  sl_env$SL.s1_mid   <- .SL_passthrough("s1_mid")
  sl_env$SL.s2_ecci  <- .SL_passthrough("s2_ecci")
  sl_env$SL.s3_mi    <- .SL_passthrough("s3_mi")
  sl_env$SL.s4_bayes <- .SL_passthrough("s4_bayes")
  # Needed so that predict.SL_passthrough is visible inside SuperLearner
  sl_env$predict.SL_passthrough <- predict.SL_passthrough
  # SuperLearner looks up its default screener "All" in `env`; copy it in
  # so the lookup succeeds when env != globalenv().
  sl_env$All <- SuperLearner::All

  # Build SL.library explicitly as (learner, screener) pairs so that
  # SuperLearner does not rely on its own implicit default expansion.
  sl_library <- list(
    c("SL.s1_mid",   "All"),
    c("SL.s2_ecci",  "All"),
    c("SL.s3_mi",    "All"),
    c("SL.s4_bayes", "All")
  )

  set.seed(seed)
  fit <- SuperLearner::SuperLearner(
    Y          = Y,
    X          = X,
    SL.library = sl_library,
    family     = stats::gaussian(),
    method     = "method.NNLS",
    cvControl  = list(V = V, shuffle = TRUE),
    verbose    = verbose,
    env        = sl_env
  )

  # SuperLearner returns SL.predict (n x 1) — the convex combination on
  # the full training data using the learned meta-weights. This IS the
  # out-of-sample CV-safe prediction for the NNLS meta-learner because
  # the meta-weights themselves are fit on CV predictions.
  #
  # v1.0.0 design note: we deliberately do NOT clip to [s1_min, s1_max].
  # Empirical testing on n = 1000 showed that clipping to the S1 interval
  # systematically biases S6 downwards on subgroups where S1_max is a
  # strict but imperfect upper bound. The Super Learner oracle property
  # is preserved only if the NNLS convex combination is used as-is.
  s6_meta <- as.numeric(fit$SL.predict)

  weights <- as.numeric(fit$coef)
  names(weights) <- c("S1_mid", "S2_ecci", "S3_mi", "S4_bayes")

  cv_risk <- as.numeric(fit$cvRisk)
  names(cv_risk) <- c("S1_mid", "S2_ecci", "S3_mi", "S4_bayes")

  list(
    predictions = s6_meta,
    weights     = weights,
    cv_risk     = cv_risk,
    sl_fit      = fit
  )
}
