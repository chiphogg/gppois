%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  R/CovarianceNoise.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{CovarianceNoise}
\docType{class}
\alias{CovarianceNoise}


\title{Covariance describing random noise}

\description{
  This subclass lets us treat noise in a unified way within our Model.

  Package:   \cr
\bold{Class CovarianceNoise}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link{Covariance}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{CovarianceNoise}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{CovarianceNoise}\cr
extends \link{Covariance}\cr


}

\usage{CovarianceNoise(id="noise", sigma=NA, sigma.bounds=NA, ...)}

\arguments{
  \item{id}{(character) A string to identify this covariance object.}
  \item{sigma}{(numeric) The default value of the noise.}
  \item{sigma.bounds}{(numeric) The range of values which \code{sigma} might assume.}
  \item{...}{Not used.}
}

\section{Covariance Parameters}{
  This section lists the fit parameters corresponding to this type of
  Covariance.  Any parameters marked as \dQuote{(Scale parameter)} will be
  optimized in log-space, consistent with the Jeffreys prior.

  \describe{
    \item{sigma}{(Scale parameter) The magnitude (standard deviation) of the
       noise.}
  }
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[R:getLogspaceNames.CovarianceNoise]{getLogspaceNames}} \tab Names of "scale"-type parameters.\cr
 \tab \code{\link[R:getLowerPlain.CovarianceNoise]{getLowerPlain}} \tab Lower bounds for params, with plain names.\cr
 \tab \code{\link[R:getParamNamesPlain.CovarianceNoise]{getParamNamesPlain}} \tab Basenames of parameters.\cr
 \tab \code{\link[R:getParamsPlain.CovarianceNoise]{getParamsPlain}} \tab Parameter values with plain names.\cr
 \tab \code{\link[R:getUpperPlain.CovarianceNoise]{getUpperPlain}} \tab Upper bounds for params, with plain names.\cr
 \tab \code{\link[R:KDerivImplementation.CovarianceNoise]{KDerivImplementation}} \tab Element-wise derivatives of Covariance matrix.\cr
 \tab \code{\link[R:K.specific.CovarianceNoise]{K.specific}} \tab Noise Covariance matrix.\cr
 \tab \code{paramsPlainImplementation} \tab  -\cr
 \tab \code{setLowerPlain} \tab  -\cr
 \tab \code{setUpperPlain} \tab  -\cr
 \tab \code{\link[R:Variance.CovarianceNoise]{Variance}} \tab Noise variance at each point.\cr
 \tab \code{\link[R:getLogspaceNames.CovarianceNoise]{getLogspaceNames}} \tab Names of "scale"-type parameters.\cr
 \tab \code{\link[R:getLowerPlain.CovarianceNoise]{getLowerPlain}} \tab Lower bounds for params, with plain names.\cr
 \tab \code{\link[R:getParamNamesPlain.CovarianceNoise]{getParamNamesPlain}} \tab Basenames of parameters.\cr
 \tab \code{\link[R:getParamsPlain.CovarianceNoise]{getParamsPlain}} \tab Parameter values with plain names.\cr
 \tab \code{\link[R:getUpperPlain.CovarianceNoise]{getUpperPlain}} \tab Upper bounds for params, with plain names.\cr
 \tab \code{\link[R:KDerivImplementation.CovarianceNoise]{KDerivImplementation}} \tab Element-wise derivatives of Covariance matrix.\cr
 \tab \code{\link[R:K.specific.CovarianceNoise]{K.specific}} \tab Noise Covariance matrix.\cr
 \tab \code{paramsPlainImplementation} \tab  -\cr
 \tab \code{setLowerPlain} \tab  -\cr
 \tab \code{setUpperPlain} \tab  -\cr
 \tab \code{\link[R:Variance.CovarianceNoise]{Variance}} \tab Noise variance at each point.\cr
}


 \bold{Methods inherited from Covariance}:\cr
clone, EncodeForTraining, FixConstParam, getId, getLower, getParamNames, getParams, getUpper, KInIn, KInInDeriv, KInOut, KOutIn, KOutOut, PrependId, print, setId, setLower, setParams, setParamsPlain, setUpper, clone, EncodeForTraining, FixConstParam, getId, getLower, getParamNames, getParams, getUpper, KInIn, KInInDeriv, KInOut, KOutIn, KOutOut, PrependId, print, setId, setLower, setParams, setParamsPlain, setUpper

 \bold{Methods inherited from Object}:\cr
as.character, attachLocally, attach, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, [[<-, [[, $<-, $, objectSize, print, registerFinalizer, save


}

\author{Charles R. Hogg III}
\keyword{classes}
