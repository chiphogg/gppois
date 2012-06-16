\name{getUpperPlain.CovarianceSEAniso2D}
\alias{CovarianceSEAniso2D$upperPlain}
\alias{getUpperPlain.CovarianceSEAniso2D}
\alias{setUpperPlain.CovarianceSEAniso2D}
\title{Upper bounds for params, with plain names}
\arguments{
  \item{U}{A (named) vector of new upper bounds (we ONLY
  use ones which are named, and whose names match up with
  names of parameters.)}

  \item{...}{Not used.}
}
\value{
  The upper bounds for the parameters for this covariance
  function, but with names undecorated by its id.
}
\description{
  Gives a vector of upper bounds for the parameter values,
  whose names are NOT decorated by the id of this
  Covariance object.
}
\seealso{
  \code{\link{getLowerPlain.CovarianceSEAniso2D}}

  \code{\link{CovarianceSEAniso2D}}
}

