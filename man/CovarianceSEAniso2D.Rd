\name{CovarianceSEAniso2D}
\alias{CovarianceSEAniso2D}
\title{2D anisotropic SE covariance}
\arguments{
  \item{id}{(character) A string to identify this
  covariance object.}

  \item{ell.1}{(numeric) One characteristic horizontal
  scale for features in the functions.}

  \item{ell.2}{(numeric) Another characteristic horizontal
  scale for features in the functions.}

  \item{theta.1}{(numeric) The angle of the ell.1
  direction.}

  \item{sigma.f}{(numeric) A characteristic vertical scale
  for features in functions being modeled.}

  \item{ell.1.bounds}{(numeric) The range of values which
  \code{ell.1} might assume.}

  \item{ell.2.bounds}{(numeric) The range of values which
  \code{ell.2} might assume.}

  \item{theta.1.bounds}{(numeric) The range of values which
  \code{theta.1} might assume.}

  \item{sigma.f.bounds}{(numeric) The range of values which
  \code{sigma.f} might assume.}

  \item{...}{Not used.}
}
\description{
  A 2D squared-exponential covariance whose eigenvalues of
  the covariance matrix are not assumed identical.  In
  other words, features might vary more rapidly in one
  direction than in the orthogonal direction.
}
\section{Covariance Parameters}{
  This section lists the fit parameters corresponding to
  this type of Covariance.  Any parameters marked as
  \dQuote{(Scale parameter)} will be optimized in
  log-space, consistent with the Jeffreys prior.

  \describe{ \item{ell.1}{(Scale parameter) One horizontal
  feature lengthscale.} \item{ell.2}{(Scale parameter)
  Another horizontal feature lengthscale.}
  \item{theta.1}{The angle of the ell.1 axis}
  \item{sigma.f}{(Scale parameter) The vertical feature
  lengthscale.} }
}

