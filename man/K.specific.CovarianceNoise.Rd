\name{K.specific.CovarianceNoise}
\alias{K.specific.CovarianceNoise}
\title{Noise Covariance matrix}
\arguments{
  \item{X}{X-values for the input points (i.e., where we
  have data)}

  \item{X.out}{X-values for the points desired to predict}

  \item{...}{Not used.}
}
\value{
  The covariance matrix taking \code{X} into \code{X.out},
  based on the parameter values in \code{this}.
}
\description{
  Calculates a covariance matrix for the noise covariance
  specifically.
}
\note{
  This is very similar to a \code{\link{CovarianceSE}} in
  the limit as \code{ell} goes to zero.  However, there is
  a very important difference!  Consider the case of
  multiple datapoints at the same X-value.  CovarianceSE
  assigns a correlation of 1; CovarianceNoise assigns a
  correlation of 0.  Caution is required!
}
\seealso{
  \code{\link{CovarianceNoise}}
}

