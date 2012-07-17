\name{print.Covariance}
\alias{print.Covariance}
\title{Pretty-printing for Covariance objects}
\usage{
  \method{print}{Covariance} (x, indent = 0, ...)
}
\arguments{
  \item{x}{The Covariance object to print.}

  \item{indent}{Aids in formatting: the number of spaces to
  print before every line.}

  \item{\dots}{Not used.}
}
\description{
  Prints out the id of the Covariance, followed by a list
  of its parameters. For each parameter, we print the name,
  lower bound, current value, and upper bound.
}
\seealso{
  \code{\link{Covariance}}
}

