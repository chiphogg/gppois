%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  R/CovarianceSEVaryingEll.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{CovarianceSEVaryingEll}
\docType{class}
\alias{CovarianceSEVaryingEll}


\title{Nonstationary Squared-Exponential Covariance}

\description{
  A \emph{nonstationary} squared-exponential covariance (meaning that the
  parameters are \emph{no longer} independent of the covariates).  It is
  governed by the two usual parameters (i.e., the horizontal and vertical
  lengthscales \code{ell} and \code{sigma.f}), but \code{ell} now
  depends on \code{X}.

  Note that it cannot easily be optimized within the current paradigm!
  Usually in \pkg{gppois}, a Model is trained on a Dataset by varying the
  hyperparameters, which are assumed to be small and finite in number.  But
  now \code{ell} is a continuous function -- there are infinitely many
  hyperparameters!  Even if we only care about the denoised values of each
  datapoint, we still have one hyperparameter per datapoint, which is far too
  many for most optimizers.

  The solution is to optimize it with a different paradigm: the \emph{focus
  regions} approach (Hogg et al. 2012).  If \code{ell(X)} varies slowly, we
  can break our function into pieces, perform local fits, and interpolate the
  resulting values of \code{ell}.  This functionality is not yet integrated
  into the main package.

  Package:   \cr
\bold{Class CovarianceSEVaryingEll}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link{Covariance}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{CovarianceSEVaryingEll}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{CovarianceSEVaryingEll}\cr
extends \link{Covariance}\cr


}

\usage{CovarianceSEVaryingEll(..., id="SEVaryingEll", X.ell=NA, ell=NA, sigma.f=NA)}

\arguments{
  \item{id}{(character) A string to identify this covariance object.}
  \item{ell}{(numeric) A characteristic horizontal scale for features in
     functions being modeled.}
  \item{sigma.f}{(numeric) A characteristic vertical scale for features in
     functions being modeled.}
  \item{ell.bounds}{(numeric) The range of values which \code{ell} might
     assume.}
  \item{sigma.f.bounds}{(numeric) The range of values which \code{sigma.f}
     might assume.}
  \item{...}{Not used.}
}

\section{Covariance Parameters}{
  This section lists the fit parameters corresponding to this type of
  Covariance.  Any parameters marked as \dQuote{(Scale parameter)} will be
  optimized in log-space, consistent with the Jeffreys prior.

  \describe{
    \item{ell}{(Scale parameter) The horizontal feature lengthscale.}
    \item{sigma.f}{(Scale parameter) The vertical feature lengthscale.}
  }
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[R:ell.CovarianceSEVaryingEll]{ell}} \tab ell(X).\cr
 \tab \code{\link[R:getLogspaceNames.CovarianceSEVaryingEll]{getLogspaceNames}} \tab Names of "scale"-type parameters.\cr
 \tab \code{\link[R:getLowerPlain.CovarianceSEVaryingEll]{getLowerPlain}} \tab Lower bounds for params, with plain names.\cr
 \tab \code{\link[R:getParamNamesPlain.CovarianceSEVaryingEll]{getParamNamesPlain}} \tab Basenames of parameters.\cr
 \tab \code{\link[R:getParamsPlain.CovarianceSEVaryingEll]{getParamsPlain}} \tab Parameter values with plain names.\cr
 \tab \code{\link[R:getUpperPlain.CovarianceSEVaryingEll]{getUpperPlain}} \tab Upper bounds for params, with plain names.\cr
 \tab \code{\link[R:KDerivImplementation.CovarianceSEVaryingEll]{KDerivImplementation}} \tab Element-wise derivatives of Covariance matrix.\cr
 \tab \code{\link[R:K.specific.CovarianceSEVaryingEll]{K.specific}} \tab Nonstationary Squared-exponential Covariance matrix.\cr
 \tab \code{paramsPlainImplementation} \tab  -\cr
 \tab \code{setLowerPlain} \tab  -\cr
 \tab \code{setUpperPlain} \tab  -\cr
 \tab \code{\link[R:sigma.f.CovarianceSEVaryingEll]{sigma.f}} \tab sigma.f(X).\cr
 \tab \code{\link[R:Variance.CovarianceSEVaryingEll]{Variance}} \tab Nonstationary SE variance at each point.\cr
 \tab \code{\link[R:ell.CovarianceSEVaryingEll]{ell}} \tab ell(X).\cr
 \tab \code{\link[R:getLogspaceNames.CovarianceSEVaryingEll]{getLogspaceNames}} \tab Names of "scale"-type parameters.\cr
 \tab \code{\link[R:getLowerPlain.CovarianceSEVaryingEll]{getLowerPlain}} \tab Lower bounds for params, with plain names.\cr
 \tab \code{\link[R:getParamNamesPlain.CovarianceSEVaryingEll]{getParamNamesPlain}} \tab Basenames of parameters.\cr
 \tab \code{\link[R:getParamsPlain.CovarianceSEVaryingEll]{getParamsPlain}} \tab Parameter values with plain names.\cr
 \tab \code{\link[R:getUpperPlain.CovarianceSEVaryingEll]{getUpperPlain}} \tab Upper bounds for params, with plain names.\cr
 \tab \code{\link[R:KDerivImplementation.CovarianceSEVaryingEll]{KDerivImplementation}} \tab Element-wise derivatives of Covariance matrix.\cr
 \tab \code{\link[R:K.specific.CovarianceSEVaryingEll]{K.specific}} \tab Nonstationary Squared-exponential Covariance matrix.\cr
 \tab \code{paramsPlainImplementation} \tab  -\cr
 \tab \code{setLowerPlain} \tab  -\cr
 \tab \code{setUpperPlain} \tab  -\cr
 \tab \code{\link[R:sigma.f.CovarianceSEVaryingEll]{sigma.f}} \tab sigma.f(X).\cr
 \tab \code{\link[R:Variance.CovarianceSEVaryingEll]{Variance}} \tab Nonstationary SE variance at each point.\cr
}


 \bold{Methods inherited from Covariance}:\cr
clone, EncodeForTraining, FixConstParam, getId, getLower, getParamNames, getParams, getUpper, KInIn, KInInDeriv, KInOut, KOutIn, KOutOut, PrependId, print, setId, setLower, setParams, setParamsPlain, setUpper, clone, EncodeForTraining, FixConstParam, getId, getLower, getParamNames, getParams, getUpper, KInIn, KInInDeriv, KInOut, KOutIn, KOutOut, PrependId, print, setId, setLower, setParams, setParamsPlain, setUpper

 \bold{Methods inherited from Object}:\cr
as.character, attachLocally, attach, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, [[<-, [[, $<-, $, objectSize, print, registerFinalizer, save


}

\references{
   Hogg, C., K. Mullen, and I. Levin (2012). A Bayesian approach for
   denoising one-dimensional data. Journal of Applied Crystallography, 45(3),
   pp. 471-481.
}

\author{Charles R. Hogg III}
\keyword{classes}
