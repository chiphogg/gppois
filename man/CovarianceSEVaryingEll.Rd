\name{CovarianceSEVaryingEll}
\alias{CovarianceSEVaryingEll}
\title{Nonstationary Squared-Exponential Covariance}
\arguments{
  \item{id}{(character) A string to identify this
  covariance object.}

  \item{X.ell}{(numeric) The X-values where \code{ell} is
  specified.}

  \item{ell}{(numeric) A series of local values for the
  characteristic horizontal scale for features in functions
  being modeled.}

  \item{sigma.f}{(numeric) A characteristic vertical scale
  for features in functions being modeled.}

  \item{...}{Not used.}
}
\description{
  A \emph{nonstationary} squared-exponential covariance
  (meaning that the parameters are \emph{no longer}
  independent of the covariates).  It is governed by the
  two usual parameters (i.e., the horizontal and vertical
  lengthscales \code{ell} and \code{sigma.f}), but
  \code{ell} now depends on \code{X}.
}
\details{
  Note that it cannot easily be optimized within the
  current paradigm!  Usually in \pkg{gppois}, a Model is
  trained on a Dataset by varying the hyperparameters,
  which are assumed to be small and finite in number.  But
  now \code{ell} is a continuous function -- there are
  infinitely many hyperparameters!  Even if we only care
  about the denoised values of each datapoint, we still
  have one hyperparameter per datapoint, which is far too
  many for most optimizers.

  The solution is to optimize it with a different paradigm:
  the \emph{focus regions} approach (Hogg et al. 2012).  If
  \code{ell(X)} varies slowly, we can break our function
  into pieces, perform local fits, and interpolate the
  resulting values of \code{ell}.  This functionality is
  not yet integrated into the main package.
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
\references{
  Hogg, C., K. Mullen, and I. Levin (2012). A Bayesian
  approach for denoising one-dimensional data. Journal of
  Applied Crystallography, 45(3), pp. 471-481.
}

