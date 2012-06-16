\name{PosteriorInterval.Model}
\alias{PosteriorInterval.Model}
\title{Best estimate, including uncertainty}
\arguments{
  \item{this}{The Model object.}

  \item{d}{The Dataset to train the Model on.}

  \item{X.out}{(matrix) The X-points where we want to
  predict the function.}

  \item{num.sd}{The number of standard deviations from the
  mean our interval should include (defaults to 1).}

  \item{contributions}{(character vector) Id's of the
  contributing Covariances we want to predict (other
  contributions are considered noise); default is every
  contribution not named 'noise'.}

  \item{...}{Not used.}
}
\value{
  A numeric vector with optimal predictions at every point
  in X.out.
}
\description{
  Computes bounds on the uncertainty (in terms of the
  standard deviation) in the prediction at a given point,
  along with the prediction.
}
\seealso{
  \code{\link{PredictionMatrix.Model}}

  \code{\link{PosteriorMean.Model}}

  \code{\link{Model}}
}

