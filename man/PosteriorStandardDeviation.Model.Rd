\name{PosteriorStandardDeviation.Model}
\alias{PosteriorStandardDeviation.Model}
\title{Pointwise uncertainty}
\arguments{
  \item{this}{The Model object.}

  \item{d}{The Dataset to train the Model on.}

  \item{X.out}{(matrix) the X-points where we want to
  predict the function.}

  \item{contributions}{(character vector) Id's of the
  contributing Covariances we want to predict (other
  contributions are considered noise); default is every
  contribution not named 'noise'.}

  \item{...}{Not used.}
}
\value{
  A numeric vector with the posterior standard deviation at
  every point in X.out.
}
\description{
  Computes the posterior \dQuote{sigma} at a given point.
}
\seealso{
  \code{\link{PredictionMatrix.Model}}

  \code{\link{PosteriorMean.Model}}

  \code{\link{Model}}
}

