\name{K.Covariance}
\alias{K.Covariance}
\alias{KInIn}
\alias{KInOut}
\alias{KOutIn}
\alias{KOutOut}
\title{Covariance matrix}
\arguments{
  \item{d}{A Dataset object encapsulating the data to train
  on.}

  \item{X.out}{A matrix (with d$d columns) of X-locations
  where we want predictions.}

  \item{...}{Not used.}
}
\value{
  A \code{ncol(X2) x ncol(X1)} matrix: the covariance from
  each input point to each output point.
}
\description{
  Gives the covariance matrix K between an input X1 and
  output X2, where K[i, j] is the covariance function
  evaluated at datapoints X1[i] and X2[j]:
}
\details{
  \describe{ \item{\code{KInIn}:}{X1 = d$X, X2 = d$X}
  \item{\code{KInOut}:}{X1 = X.out, X2 = d$X}
  \item{\code{KOutIn}:}{X1 = d$X, X2 = X.out}
  \item{\code{KOutOut}:}{X1 = X.out, X2 = X.out} }
}
\seealso{
  \code{\link{KInInDeriv}}

  \code{\link{Covariance}}
}

