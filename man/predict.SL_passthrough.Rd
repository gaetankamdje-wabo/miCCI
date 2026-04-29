\name{predict.SL_passthrough}
\alias{predict.SL_passthrough}
\title{Predict method for pass-through SuperLearner learners}
\description{S3 predict method used internally by cci_meta_fit().}
\usage{\method{predict}{SL_passthrough}(object, newdata, ...)}
\arguments{
  \item{object}{Fitted pass-through learner.}
  \item{newdata}{data.frame with the column referenced by object$col_name.}
  \item{...}{Unused.}
}
\value{Numeric vector of predictions.}

