\name{PosteriorMean.Model}
\alias{PosteriorMean.Model}
\title{Best estimate of the true function}
\arguments{
  \item{this}{The Model object.}

  \item{d}{The Dataset to train the Model on.}

  \item{X.out}{(matrix) The X-points where we want to
  predict the function.}

  \item{contributions}{(character vector) Id's of the
  contributing Covariances we want to predict (other
  contributions are considered noise); default is every
  contribution not named 'noise'.}

  \item{untransform.result}{logical; if TRUE, we transform
  back to the space of \code{\link{dpts}} (as opposed to
  the space of \code{\link{xformedDpts}} where training
  takes place).}

  \item{...}{Not used.}
}
\value{
  A numeric vector with optimal predictions at every point
  in X.out.
}
\description{
  Computes this Model's optimal prediction of the
  underlying function's value at every point in 'X.out'.
}
\seealso{
  \code{\link{Model}}
}

