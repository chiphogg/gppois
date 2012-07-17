\name{paramNamesPlain}
\alias{getParamNamesPlain}
\alias{paramNamesPlain}
\title{Undecorated parameter names}
\description{
  Covariance objects usually decorate their parameter names
  by prepending the id.  This lets us combine several
  Covariance objects into a single model, without worrying
  about parameter names colliding.  Often, though, we need
  to refer to the \dQuote{plain} names -- something like
  \code{ell} instead of \code{SE_1.ell}.  That is what
  these \dQuote{plain-named} versions are for.
}
\seealso{
  \code{\link{lowerPlain}}

  \code{\link{paramsPlain}}

  \code{\link{upperPlain}}
}

