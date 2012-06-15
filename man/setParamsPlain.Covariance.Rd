\name{setParamsPlain.Covariance}
\alias{Covariance$paramsPlain}
\alias{getParamsPlain}
\alias{getParamsPlain.Covariance}
\alias{paramsPlain}
\alias{setParamsPlain.Covariance}
\title{Set parameters using undecorated names}
\arguments{
  \item{p}{A (named) vector of new parameter values (we
  \emph{only} use ones which are named, and whose names
  match up with names of parameters.)}

  \item{...}{Not used.}
}
\description{
  Set parameters of this covariance using "undecorated" or
  "plain" names.
}
\details{
  \section{Undecorated Names}{ Covariance objects decorate
  their parameter names with their ID: e.g., a Covariance
  with id "SE" will have a parameter named "SE.ell".  This
  is important when you have multiple Covariance objects;
  decorating with the ID keeps the names unique.  However,
  it does make it less convenient to refer to them.  This
  function allows to access them with the more intuitive,
  \emph{undecorated} names.

  I do not provide anything like the \code{for.training}
  option from \code{\link{getParams}}, because it doesn't
  make any sense.  The training functions always have to
  use the \emph{decorated} names, or else risk name
  collisions. }

  \section{Note for developers}{ For developers who plan to
  write new Covariance classes: the corresponding
  \emph{accessor} method is the responsibility of
  subclasses.  i.e., you must define a method
  "getParamsPlain" for your class.  (I suggest looking at
  the corresponding methods for the
  \code{\link{CovarianceSE}} class to get started.) }
}
\seealso{
  \code{\link{getParamsPlain}}

  \code{\link{getLowerPlain}}

  \code{\link{setLowerPlain}}

  \code{\link{Covariance}}
}

