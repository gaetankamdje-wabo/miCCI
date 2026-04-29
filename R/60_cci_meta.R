# =============================================================================
# miCCI / 60_cci_meta.R
# Meta-learner: cross-validated Super Learner over S1..S4.
#
# The Super Learner (van der Laan, Polley & Hubbard 2007) is a
# cross-validated convex ensemble with an asymptotic oracle property:
# given a library of base estimators, its out-of-sample risk is at least
# as small as the risk of the best single estimator in the library.
#
# In miCCI, the library consists of pass-through wrappers around S1 mid,
# S2 expected, S3 multiple-imputation mean and S4 Bayesian posterior
# median. SuperLearner::SuperLearner() fits a non-negative least-squares
# meta-weight vector over 10-fold cross-validated out-of-fold predictions.
# Weights are non-negative and sum to one.
#
# No train/calibrate/validate split is used: the entire cohort is passed
# in, and the internal 10-fold CV guarantees that each encounter's final
# meta-prediction does not depend on its own gold value.
#
# Output column name across the package: `meta`.
# =============================================================================

#' Pass-through SuperLearner learner: returns a single column of newX
#' as its prediction, used to make the SuperLearner a pure meta-combination
#' over the four base-strategy columns without any additional feature
#' learning.
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

#' S3 predict method for pass-through SuperLearner learners.
#' @param object  fitted pass-through learner.
#' @param newdata data.frame with the column referenced by `object$col_name`.
#' @param ...     unused.
#' @export
predict.SL_passthrough <- function(object, newdata, ...) {
  as.numeric(newdata[[object$col_name]])
}

#' Fit the cross-validated meta-learner.
#'
#' @param dt   data.table with columns `cci_gold`, `s1_min`, `s1_max`,
#'   `s1_mid`, `s2_ecci`, `s3_mi`, `s4_bayes`.
#' @param V    number of cross-validation folds (default 10).
#' @param seed RNG seed for fold assignment reproducibility.
#' @param verbose if TRUE, SuperLearner prints fold-level progress.
#'
#' @return list with elements
#' \describe{
#'   \item{predictions}{numeric vector of meta predictions, length nrow(dt).}
#'   \item{weights}{named numeric vector of length 4 - NNLS weights over
#'     S1 mid, S2 ecci, S3 mi, S4 bayes.}
#'   \item{cv_risk}{named numeric vector - cross-validated risk of each
#'     base learner.}
#'   \item{sl_fit}{the raw SuperLearner fit object.}
#' }
#' @export
cci_meta_fit <- function(dt, V = 10L, seed = 42L, verbose = FALSE) {
  if (!isTRUE(requireNamespace("SuperLearner", quietly = TRUE)))
    stop("cci_meta_fit() requires the 'SuperLearner' package.")

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

  # Local environment for SuperLearner's by-name lookups; never touches
  # the user's global environment.
  sl_env <- new.env(parent = globalenv())
  sl_env$SL.s1_mid   <- .SL_passthrough("s1_mid")
  sl_env$SL.s2_ecci  <- .SL_passthrough("s2_ecci")
  sl_env$SL.s3_mi    <- .SL_passthrough("s3_mi")
  sl_env$SL.s4_bayes <- .SL_passthrough("s4_bayes")
  sl_env$predict.SL_passthrough <- predict.SL_passthrough
  sl_env$All <- SuperLearner::All

  sl_library <- list(
    c("SL.s1_mid",   "All"),
    c("SL.s2_ecci",  "All"),
    c("SL.s3_mi",    "All"),
    c("SL.s4_bayes", "All")
  )

  fit <- .with_local_seed(seed, {
    SuperLearner::SuperLearner(
      Y          = Y,
      X          = X,
      SL.library = sl_library,
      family     = stats::gaussian(),
      method     = "method.NNLS",
      cvControl  = list(V = V, shuffle = TRUE),
      verbose    = verbose,
      env        = sl_env
    )
  })

  weights <- as.numeric(fit$coef)
  names(weights) <- c("S1_mid", "S2_ecci", "S3_mi", "S4_bayes")
  cv_risk <- as.numeric(fit$cvRisk)
  names(cv_risk) <- c("S1_mid", "S2_ecci", "S3_mi", "S4_bayes")

  list(
    predictions = as.numeric(fit$SL.predict),
    weights     = weights,
    cv_risk     = cv_risk,
    sl_fit      = fit
  )
}
