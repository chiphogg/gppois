\name{upperPlain}
\alias{getUpperPlain}
\alias{setUpperPlain}
\alias{upperPlain}
\title{Upper bounds for plain-named parameters}
\description{
  Objects governed by parameters also have boundaries on
  those parameters. These R.oo functions retrieve and/or
  set these upper bounds.
}
\details{
  Covariance objects usually decorate their parameter names
  by prepending the id.  This lets us combine several
  Covariance objects into a single model, without worrying
  about parameter names colliding.  Often, though, we need
  to refer to the \dQuote{plain} names -- something like
  \code{ell} instead of \code{SE_1.ell}.  That is what
  these \dQuote{plain-named} versions are for.

  Boundaries push obstacles aside when they are changed.
  For example, suppose a parameter and its bounds are set
  as (lower=3, param=4, upper=5).  Setting lower to 6 will
  result in values (lower=6, param=6, upper=6): the other
  values are "dragged along".
}
\seealso{
  \code{\link{lowerPlain}}

  \code{\link{paramsPlain}}

  \code{\link{upper}}
}

