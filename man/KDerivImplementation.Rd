\name{KDerivImplementation}
\alias{KDerivImplementation}
\title{Element-wise derivatives of Covariance matrix}
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
  "Pure virtual" method to calculate the element-wise
  derivative of \code{\link{KInIn}}, with respect to the
  parameter whose (plain) name is \code{param}.
}

