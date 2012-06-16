\name{K.specific.CovarianceSE}
\alias{K.specific.CovarianceSE}
\title{Squared-Exponential Covariance matrix}
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
  Calculates a covariance matrix for the
  squared-exponential covariance function.
}
\seealso{
  \code{\link{CovarianceSE}}
}

