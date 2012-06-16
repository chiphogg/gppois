\name{L.Model}
\alias{L.Model}
\title{Lower Cholesky root of covariance matrix}
\arguments{
  \item{this}{The Model object.}

  \item{d}{The Dataset we're training on.}

  \item{X.out}{A matrix (with d$d columns) of X-locations
  where we want predictions.}

  \item{contributions}{The names of the parameters which
  are considered \dQuote{signal}.}

  \item{...}{Not used.}
}
\value{
  The lower-triangular Cholesky decomposition of the
  covariance matrix K.  (i.e., the matrix L such that L %*%
  t(L) = K.)
}
\description{
  Compute the lower-triangular Cholesky decomposition of
  the covariance matrix. This is useful for taking random
  draws from the posterior.
}
\seealso{
  \code{\link{Model}}
}

