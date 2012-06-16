\name{print.Model}
\alias{print.Model}
\title{Pretty-printing for Model objects}
\usage{
  \method{print}{Model} (this, indent = 0, ...)
}
\arguments{
  \item{this}{The Model object to print.}

  \item{indent}{Aids in formatting: the number of spaces to
  print before every line.}

  \item{...}{Not used.}
}
\description{
  Prints out the id of the Model, followed by a list of its
  contributing Covariances and their IDs.  Then, prints
  each parameter, including the name, lower bound, current
  value, and upper bound.
}
\seealso{
  \code{\link{Model}}
}

