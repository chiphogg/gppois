\name{getLowerPlain.CovarianceSEVaryingEll}
\alias{CovarianceSEVaryingEll$lowerPlain}
\alias{getLowerPlain.CovarianceSEVaryingEll}
\alias{setLowerPlain.CovarianceSEVaryingEll}
\title{Lower bounds for params, with plain names}
\arguments{
  \item{L}{A (named) vector of new lower bounds (we ONLY
  use ones which are named, and whose names match up with
  names of parameters.)}

  \item{...}{Not used.}
}
\value{
  The lower bounds for the parameters for this covariance
  function, but with names undecorated by its id.
}
\description{
  Gives a vector of lower bounds for the parameter values,
  whose names are NOT decorated by the id of this
  Covariance object.
}
\seealso{
  \code{\link{getUpperPlain.CovarianceSEVaryingEll}}

  \code{\link{CovarianceSEVaryingEll}}
}

