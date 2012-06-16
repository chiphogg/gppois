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
 \tab \code{getLogspaceNames} \tab  -\cr
 \tab \code{getLowerPlain} \tab  -\cr
 \tab \code{getParamNamesPlain} \tab  -\cr
 \tab \code{getParamsPlain} \tab  -\cr
 \tab \code{getUpperPlain} \tab  -\cr
 \tab \code{KDerivImplementation} \tab  -\cr
 \tab \code{K.specific} \tab  -\cr
 \tab \code{paramsPlainImplementation} \tab  -\cr
 \tab \code{setLowerPlain} \tab  -\cr
 \tab \code{setUpperPlain} \tab  -\cr
 \tab \code{Variance} \tab  -\cr
}


 \bold{Methods inherited from Covariance}:\cr
clone, EncodeForTraining, FixConstParam, getId, getLower, getParamNames, getParams, getUpper, KInIn, KInInDeriv, KInOut, KOutIn, KOutOut, PrependId, print, setId, setLower, setParams, setParamsPlain, setUpper

 \bold{Methods inherited from Object}:\cr
as.character, attachLocally, attach, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, [[<-, [[, $<-, $, objectSize, print, registerFinalizer, save


}

\author{Charles R. Hogg III}
\keyword{classes}
