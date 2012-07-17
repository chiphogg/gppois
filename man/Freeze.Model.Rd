\name{Freeze.Model}
\alias{Freeze}
\alias{Freeze.Model}
\title{Make some parameters constant}
\arguments{
  \item{this}{The Model object.}

  \item{p.names}{The names of the parameters to freeze.}

  \item{...}{Not used.}
}
\description{
  Makes a subset of parameters constant, by setting the
  upper and lower bounds equal to the current parameter
  value.  Defaults to freezing all parameters.
}
\seealso{
  \code{\link{Model}}
}

