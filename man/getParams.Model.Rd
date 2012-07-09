\name{getParams.Model}
\alias{Model$params}
\alias{getParams.Model}
\alias{setParams.Model}
\title{Parameters for the Model}
\arguments{
  \item{this}{The Model object.}

  \item{p}{A (named) vector of new parameter values (we
  ONLY use ones which are named, and whose names match up
  with names of parameters.)}

  \item{for.training}{If TRUE, we ignore "constant"
  parameters (i.e., where lower=upper) and return the *log*
  of any "scale" parameters (such as ell or sigma.f for the
  SE model).}

  \item{...}{Not used.}
}
\value{
  A vector with the current values for each parameter.
}
\description{
  A named vector of parameters governing this Model.  Note
  that the names are \emph{decorated} by prepending the
  Covariance id; this prevents namespace collisions.
}
\seealso{
  \code{\link{getParams.Covariance}} for more about
  \code{for.training}

  \code{\link{getUpper.Model}}

  \code{\link{getLower.Model}}

  \code{\link{Model}}
}

