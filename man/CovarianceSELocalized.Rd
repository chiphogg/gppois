%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  R/CovarianceSELocalized.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{CovarianceSELocalized}
\docType{class}
\alias{CovarianceSELocalized}


\title{Localized Squared-exponential Covariance}

\description{
  A squared-exponential covariance which is limited in spatial extent.
  In addition to the two usual parameters (i.e., the horizontal and vertical
  lengthscales), there are two \dQuote{boundary} parameters, \code{X.L} and
  \code{X.R}.  \code{sigma.f} transitions smoothly to zero outside the region
  between \code{X.L} and \code{X.R}.  The transition lengthscale is
  \code{ell} (if it were any smaller, it could introduce sub-\code{ell}
  features).

  Package:   \cr
\bold{Class CovarianceSELocalized}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link{Covariance}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{CovarianceSELocalized}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{CovarianceSELocalized}\cr
extends \link{Covariance}\cr


}

\usage{CovarianceSELocalized(..., id="SE", ell=NA, sigma.f=NA, X.L=NA, X.R=NA, ell.bounds=NA, sigma.f.bounds=NA, X.L.bounds=NA, X.R.bounds=NA)}

\arguments{
  \item{id}{(character) A string to identify this covariance object.}
  \item{ell}{(numeric) A characteristic horizontal scale for features in
     functions being modeled.}
  \item{sigma.f}{(numeric) A characteristic vertical scale for features in
     functions being modeled.}
  \item{X.L}{(numeric) The left boundary of the localized region.}
  \item{X.R}{(numeric) The right boundary of the localized region.}
  \item{ell.bounds}{(numeric) The range of values which \code{ell} might
     assume.}
  \item{sigma.f.bounds}{(numeric) The range of values which \code{sigma.f}
     might assume.}
  \item{X.L.bounds}{(numeric) The range of values which \code{X.L} might
     assume.}
  \item{X.R.bounds}{(numeric) The range of values which \code{X.R} might
     assume.}
  \item{...}{Not used.}
}

\section{Covariance Parameters}{
  This section lists the fit parameters corresponding to this type of
  Covariance.  Any parameters marked as \dQuote{(Scale parameter)} will be
  optimized in log-space, consistent with the Jeffreys prior.

  \describe{
    \item{ell}{(Scale parameter) The horizontal feature lengthscale.}
    \item{sigma.f}{(Scale parameter) The vertical feature lengthscale.}
    \item{X.L}{The left boundary of the localized region.}
    \item{X.R}{The right boundary of the localized region.}
  }
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[R:getLogspaceNames.CovarianceSELocalized]{getLogspaceNames}} \tab Names of "scale"-type parameters.\cr
 \tab \code{\link[R:getLowerPlain.CovarianceSELocalized]{getLowerPlain}} \tab Lower bounds for params, with plain names.\cr
 \tab \code{\link[R:getParamNamesPlain.CovarianceSELocalized]{getParamNamesPlain}} \tab Basenames of parameters.\cr
 \tab \code{\link[R:getParamsPlain.CovarianceSELocalized]{getParamsPlain}} \tab Parameter values with plain names.\cr
 \tab \code{\link[R:getUpperPlain.CovarianceSELocalized]{getUpperPlain}} \tab Upper bounds for params, with plain names.\cr
 \tab \code{\link[R:KDerivImplementation.CovarianceSELocalized]{KDerivImplementation}} \tab Element-wise derivatives of Covariance matrix.\cr
 \tab \code{\link[R:K.specific.CovarianceSELocalized]{K.specific}} \tab Localized Squared-Exponential Covariance matrix.\cr
 \tab \code{paramsPlainImplementation} \tab  -\cr
 \tab \code{setLowerPlain} \tab  -\cr
 \tab \code{setUpperPlain} \tab  -\cr
 \tab \code{\link[R:Variance.CovarianceSELocalized]{Variance}} \tab Localized SE variance at each point.\cr
}


 \bold{Methods inherited from Covariance}:\cr
clone, EncodeForTraining, FixConstParam, getId, getLower, getParamNames, getParams, getUpper, KInIn, KInInDeriv, KInOut, KOutIn, KOutOut, PrependId, print, setId, setLower, setParams, setParamsPlain, setUpper, clone, EncodeForTraining, FixConstParam, getId, getLower, getParamNames, getParams, getUpper, KInIn, KInInDeriv, KInOut, KOutIn, KOutOut, PrependId, print, setId, setLower, setParams, setParamsPlain, setUpper

 \bold{Methods inherited from Object}:\cr
as.character, attachLocally, attach, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, [[<-, [[, $<-, $, objectSize, print, registerFinalizer, save


}

\author{Charles R. Hogg III}
\keyword{classes}
