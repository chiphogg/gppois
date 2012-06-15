\name{KInInDeriv.Covariance}
\alias{KInInDeriv}
\alias{KInInDeriv.Covariance}
\title{Element-wise derivatives of covariance matrix}
\arguments{
  \item{d}{The Dataset whose X-values we are training on.}

  \item{param}{The name of the parameter with respect to
  which we're differentiating.}

  \item{...}{Not used.}
}
\value{
  The element-by-element derivative of the total K-matrix
  for Dataset \code{d}, with respect to the parameter named
  \code{param}.
}
\description{
  Computes element-wise derivatives of the "In-In"
  covariance matrix (i.e., the covariance matrix
  \emph{from} the observed datapoints \code{d$X}, \emph{to}
  those same datapoints \code{d$X}).  These derivatives are
  the only ones important for computing derivatives of the
  posterior probability function.
}
\note{
  As this is the superclass, all we do here is to check
  whether this is our problem (i.e., whether this object
  actually has a parameter named \code{param}).  If not,
  just return a zero matrix.  Otherwise, we let the
  subclasses handle the actual computation.
}
\seealso{
  \code{\link{K.Covariance}}

  \code{\link{Covariance}}
}

