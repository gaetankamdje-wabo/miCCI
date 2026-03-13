# =============================================================================
# miCCI — 60_cci_ensemble.R
# Strategy 6: NNLS Ensemble Meta-Learner
# =============================================================================

#' Train NNLS ensemble weights on calibration stratum
cci_ensemble_train <- function(estimates_matrix, cci_gold) {
  requireNamespace("nnls", quietly = TRUE)
  stopifnot(ncol(estimates_matrix) == 5, length(cci_gold) == nrow(estimates_matrix))
  colnames(estimates_matrix) <- c("S1_mid", "S2_ecci", "S3_mi", "S4_bayes", "S5_ml")
  fit <- nnls::nnls(estimates_matrix, cci_gold)
  w <- fit$x; names(w) <- colnames(estimates_matrix)
  w_sum <- sum(w)
  if (w_sum > 0) w <- w / w_sum
  list(weights = round(w, 4), residual_norm = fit$deviance)
}

#' Predict ensemble CCI
#' @export
cci_ensemble <- function(estimates, weights) {
  if (is.vector(estimates) && length(estimates) == 5) estimates <- matrix(estimates, nrow = 1)
  as.numeric(estimates %*% weights)
}
