\name{CovarianceSE}
\alias{CovarianceSE}
\title{CovarianceSE: (S)quared-(E)xponential covariance}
\arguments{
  \item{id}{(character) A string to identify this
  covariance object.}

  \item{ell}{(numeric) A characteristic horizontal scale
  for features in functions being modeled.}

  \item{sigma.f}{(numeric) A characteristic vertical scale
  for features in functions being modeled.}

  \item{ell.bounds}{(numeric) The range of values which
  \code{ell} might assume.}

  \item{sigma.f.bounds}{(numeric) The range of values which
  \code{sigma.f} might assume.}

  \item{...}{Not used.}
}
\description{
  The standard squared-exponential covariance.  Governed by
  two parameters: a horizontal and a vertical lengthscale.
}
\section{Covariance Parameters}{
  This section lists the fit parameters corresponding to
  this type of Covariance.  Any parameters marked as
  \dQuote{(Scale parameter)} will be optimized in
  log-space, consistent with the Jeffreys prior.

  \describe{ \item{ell}{(Scale parameter) The horizontal
  feature lengthscale.} \item{sigma.f}{(Scale parameter)
  The vertical feature lengthscale.} }
}

