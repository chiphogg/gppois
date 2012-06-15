\name{KDerivImplementation.CovarianceNoise}
\alias{KDerivImplementation.CovarianceNoise}
\title{Element-wise derivatives of Covariance matrix}
\arguments{
  \item{d}{The Dataset whose X-values determine KInIn.}

  \item{param}{The (plain) name of the parameter with
  respect to which we're differentiating.}

  \item{...}{Not used.}
}
\value{
  A matrix whose elements are the derivatives of the
  corresponding elements in KInIn, with respect to the
  parameter \code{param}.
}
\description{
  Calculate the element-wise derivative of \code{KInIn},
  with respect to the parameter whose (plain) name is
  \code{param}.
}
\seealso{
  \code{\link{CovarianceNoise}}
}

