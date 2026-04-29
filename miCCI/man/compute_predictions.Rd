\name{compute_predictions}
\alias{compute_predictions}
\title{Stage 1 of the pipeline: run S1..S4 + Meta and persist the artefact.}
\description{Stage 1 of the pipeline: run S1..S4 + Meta and persist the artefact.}
\usage{compute_predictions(data_path, output_dir, freq_table, sample_size, mi_m = 20L, bayes_draws = 25L, alpha_0 = 10, sl_cv_folds = 10L, seed = 42L)}
\arguments{
  \item{data_path}{Argument data_path.}
  \item{output_dir}{Argument output_dir.}
  \item{freq_table}{Argument freq_table.}
  \item{sample_size}{Argument sample_size.}
  \item{mi_m}{Argument mi_m.}
  \item{bayes_draws}{Argument bayes_draws.}
  \item{alpha_0}{Argument alpha_0.}
  \item{sl_cv_folds}{Argument sl_cv_folds.}
  \item{seed}{Argument seed.}
}
\value{See package vignette and source for details.}
\examples{\dontrun{
# See README for usage examples.
}}

