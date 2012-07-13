\name{PredictionMatrix.Model}
\alias{PredictionMatrix}
\alias{PredictionMatrix.Model}
\title{Matrix connecting noisy data to true function}
\arguments{
  \item{this}{The Model object.}

  \item{d}{The Dataset to train the Model on.}

  \item{X.out}{matrix; the X-points where we want to
  predict the function.}

  \item{contributions}{character vector; id's of the
  contributing Covariances we want to predict (other
  contributions are considered noise).}

  \item{...}{Not used.}
}
\value{
  The matrix \code{M} which yields the predictions at
  \code{X.out} when multiplied by the datapoint vector.
}
\description{
  A matrix relating function values at the output points to
  function values at the input points.
}
\seealso{
  \code{\link{Model}}
}

