\name{K.specific}
\alias{K.specific}
\title{Covariance matrix implementation}
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
  "Pure virtual" method to calculate a covariance matrix
  for a specific subclass.
}

