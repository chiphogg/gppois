\name{getUpper.Model}
\alias{getUpper.Model}
\alias{Model$upper}
\alias{setUpper.Model}
\title{Upper bounds for parameters}
\arguments{
  \item{this}{The Model object.}

  \item{L}{A (named) vector of new parameter values (we
  ONLY use ones which are named, and whose names match up
  with names of parameters.)}

  \item{for.training}{If TRUE, we ignore "constant"
  parameters (i.e., where lower=upper) and return the *log*
  of any "scale" parameters (such as ell or sigma.f for the
  SE model).}

  \item{...}{Not used.}
}
\value{
  The upper bounds for the parameters for this model.
}
\description{
  Upper bounds for the parameter values.
}
\seealso{
  \code{\link{getLower.Model}}

  \code{\link{getParams.Model}}

  \code{\link{Model}}
}

