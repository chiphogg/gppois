\name{CovarianceSELocalized}
\alias{CovarianceSELocalized}
\title{Localized Squared-exponential Covariance}
\arguments{
  \item{id}{(character) A string to identify this
  covariance object.}

  \item{ell}{(numeric) A characteristic horizontal scale
  for features in functions being modeled.}

  \item{sigma.f}{(numeric) A characteristic vertical scale
  for features in functions being modeled.}

  \item{X.L}{(numeric) The left boundary of the localized
  region.}

  \item{X.R}{(numeric) The right boundary of the localized
  region.}

  \item{ell.bounds}{(numeric) The range of values which
  \code{ell} might assume.}

  \item{sigma.f.bounds}{(numeric) The range of values which
  \code{sigma.f} might assume.}

  \item{X.L.bounds}{(numeric) The range of values which
  \code{X.L} might assume.}

  \item{X.R.bounds}{(numeric) The range of values which
  \code{X.R} might assume.}

  \item{...}{Not used.}
}
\description{
  A squared-exponential covariance which is limited in
  spatial extent.  In addition to the two usual parameters
  (i.e., the horizontal and vertical lengthscales), there
  are two \dQuote{boundary} parameters, \code{X.L} and
  \code{X.R}.  \code{sigma.f} transitions smoothly to zero
  outside the region between \code{X.L} and \code{X.R}.
  The transition lengthscale is \code{ell} (if it were any
  smaller, it could introduce sub-\code{ell} features).
}
\section{Covariance Parameters}{
  This section lists the fit parameters corresponding to
  this type of Covariance.  Any parameters marked as
  \dQuote{(Scale parameter)} will be optimized in
  log-space, consistent with the Jeffreys prior.

  \describe{ \item{ell}{(Scale parameter) The horizontal
  feature lengthscale.} \item{sigma.f}{(Scale parameter)
  The vertical feature lengthscale.} \item{X.L}{The left
  boundary of the localized region.} \item{X.R}{The right
  boundary of the localized region.} }
}

