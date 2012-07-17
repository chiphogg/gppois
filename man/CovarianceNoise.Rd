\name{CovarianceNoise}
\alias{CovarianceNoise}
\title{CovarianceNoise: i.i.d. Gaussian noise}
\arguments{
  \item{id}{(character) A string to identify this
  covariance object.}

  \item{sigma}{(numeric) The default value of the noise.}

  \item{sigma.bounds}{(numeric) The range of values which
  \code{sigma} might assume.}

  \item{...}{Not used.}
}
\description{
  This subclass lets us treat noise in a unified way within
  our Model.
}
\section{Covariance Parameters}{
  This section lists the fit parameters corresponding to
  this type of Covariance.  Any parameters marked as
  \dQuote{(Scale parameter)} will be optimized in
  log-space, consistent with the Jeffreys prior.

  \describe{ \item{sigma}{(Scale parameter) The magnitude
  (standard deviation) of the noise.} }
}

