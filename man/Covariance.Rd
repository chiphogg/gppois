\name{Covariance}
\alias{Covariance}
\title{Covariance: superclass for covariance functions}
\arguments{
  \item{id}{(character) A string to identify this
  covariance object.}

  \item{...}{Not used.}
}
\description{
  \code{Covariance} is the superclass for more specific
  types (SE, Matern, etc.).  Nobody will actually make an
  object of type \code{Covariance}; instead, they will use
  one of the derived classes (such as
  \code{\link{CovarianceSE}}).
}
\details{
  { Regarding parameter names: the \emph{subclass} has the
  responsibility to provide \sQuote{plain-named} versions
  of all the virtual fields (\code{paramNamesPlain},
  \code{paramsPlain}, \code{lowerPlain}, etc.).  The
  \emph{superclass} will automatically handle the
  "decorated" versions (\code{paramNames}, \code{params},
  \code{lower}, etc.).

  All \code{Covariance} subclasses should remember the
  param-values and \code{\link{Dataset}} they last used to
  compute their K-matrix.  If they get asked to compute it
  again, they will simply return the previously-computed
  result if these values have not changed. }
}

