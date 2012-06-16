%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  R/CovarianceSE.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{CovarianceSE}
\docType{class}
\alias{CovarianceSE}


\title{(S)quared-(E)xponential Covariance}

\description{
  The standard squared-exponential covariance.  Governed by two parameters: a
  horizontal and a vertical lengthscale.

  Package:   \cr
\bold{Class CovarianceSE}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link{Covariance}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{CovarianceSE}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{CovarianceSE}\cr
extends \link{Covariance}\cr


}

\usage{CovarianceSE(..., id="SE", ell=NA, sigma.f=NA, ell.bounds=NA, sigma.f.bounds=NA)}

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
 \tab \code{\link[R:getLogspaceNames.CovarianceSE]{getLogspaceNames}} \tab Names of "scale"-type parameters.\cr
 \tab \code{\link[R:getLowerPlain.CovarianceSE]{getLowerPlain}} \tab Lower bounds for params, with plain names.\cr
 \tab \code{\link[R:getParamNamesPlain.CovarianceSE]{getParamNamesPlain}} \tab Basenames of parameters.\cr
 \tab \code{\link[R:getParamsPlain.CovarianceSE]{getParamsPlain}} \tab Parameter values with plain names.\cr
 \tab \code{\link[R:getUpperPlain.CovarianceSE]{getUpperPlain}} \tab Upper bounds for params, with plain names.\cr
 \tab \code{\link[R:KDerivImplementation.CovarianceSE]{KDerivImplementation}} \tab Element-wise derivatives of Covariance matrix.\cr
 \tab \code{\link[R:K.specific.CovarianceSE]{K.specific}} \tab Squared-Exponential Covariance matrix.\cr
 \tab \code{paramsPlainImplementation} \tab  -\cr
 \tab \code{setLowerPlain} \tab  -\cr
 \tab \code{setUpperPlain} \tab  -\cr
 \tab \code{\link[R:Variance.CovarianceSE]{Variance}} \tab SE variance at each point.\cr
}


 \bold{Methods inherited from Covariance}:\cr
clone, EncodeForTraining, FixConstParam, getId, getLower, getParamNames, getParams, getUpper, KInIn, KInInDeriv, KInOut, KOutIn, KOutOut, PrependId, print, setId, setLower, setParams, setParamsPlain, setUpper, clone, EncodeForTraining, FixConstParam, getId, getLower, getParamNames, getParams, getUpper, KInIn, KInInDeriv, KInOut, KOutIn, KOutOut, PrependId, print, setId, setLower, setParams, setParamsPlain, setUpper

 \bold{Methods inherited from Object}:\cr
as.character, attachLocally, attach, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, [[<-, [[, $<-, $, objectSize, print, registerFinalizer, save


}

\author{Charles R. Hogg III}
\keyword{classes}
