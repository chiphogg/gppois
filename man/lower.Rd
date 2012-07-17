\name{lower}
\alias{getLower}
\alias{lower}
\alias{setLower}
\title{Lower bounds for parameters}
\description{
  Objects governed by parameters also have boundaries on
  those parameters. These R.oo functions retrieve and/or
  set these lower bounds.
}
\details{
  Boundaries push obstacles aside when they are changed.
  For example, suppose a parameter and its bounds are set
  as (lower=3, param=4, upper=5).  Setting lower to 6 will
  result in values (lower=6, param=6, upper=6): the other
  values are "dragged along".
}
\seealso{
  \code{\link{upper}}

  \code{\link{params}}

  \code{\link{Covariance$upper}}

  \code{\link{Model$upper}}
}

