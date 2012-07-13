\name{params}
\alias{getParams}
\alias{params}
\alias{setParams}
\title{Parameter values}
\description{
  These R.oo functions retrieve and/or set the values of
  parameters which govern an object (i.e., typically a
  \code{\link{Covariance}}, or a \code{\link{Model}} which
  is a collection of Covariances).
}
\details{
  When parameters are moved, they cannot push the
  boundaries aside.  For example, suppose a parameter and
  its bounds are set as (lower=3, param=4, upper=5).
  Setting param to 6 will result in values (lower=3,
  param=5, upper=5): it is clamped to lie within the
  bounds.
}
\seealso{
  \code{\link{lower}}

  \code{\link{upper}}

  \code{\link{Covariance$params}}

  \code{\link{Model$params}}
}

