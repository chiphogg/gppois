%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  R/CovarianceSEAniso2D.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{CovarianceSEAniso2D}
\docType{class}
\alias{CovarianceSEAniso2D}


\title{2D Anisotropic SE Covariance}

\description{
  A 2D squared-exponential covariance whose eigenvalues of the covariance
  matrix are not assumed identical.  In other words, features might vary more
  rapidly in one direction than in the orthogonal direction.

  Package:   \cr
\bold{Class CovarianceSEAniso2D}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link{Covariance}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{CovarianceSEAniso2D}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{CovarianceSEAniso2D}\cr
extends \link{Covariance}\cr


}

\usage{CovarianceSEAniso2D(..., id="Aniso2D", ell.1=NA, ell.2=NA, theta.1=NA, sigma.f=NA, ell.1.bounds=NA, ell.2.bounds=NA, theta.1.bounds=NA, sigma.f.bounds=NA)}

\arguments{
  \item{id}{(character) A string to identify this covariance object.}
  \item{ell.1}{(numeric) One characteristic horizontal scale for features in
     the functions.}
  \item{ell.2}{(numeric) Another characteristic horizontal scale for features
     in the functions.}
  \item{theta.1}{(numeric) The angle of the ell.1 direction.}
  \item{sigma.f}{(numeric) A characteristic vertical scale for features in
     functions being modeled.}
  \item{ell.1.bounds}{(numeric) The range of values which \code{ell.1} might
     assume.}
  \item{ell.2.bounds}{(numeric) The range of values which \code{ell.2} might
     assume.}
  \item{theta.1.bounds}{(numeric) The range of values which \code{theta.1}
     might assume.}
  \item{sigma.f.bounds}{(numeric) The range of values which \code{sigma.f}
     might assume.}
  \item{...}{Not used.}
}

\section{Covariance Parameters}{
  This section lists the fit parameters corresponding to this type of
  Covariance.  Any parameters marked as \dQuote{(Scale parameter)} will be
  optimized in log-space, consistent with the Jeffreys prior.

  \describe{
    \item{ell.1}{(Scale parameter) One horizontal feature lengthscale.}
    \item{ell.2}{(Scale parameter) Another horizontal feature lengthscale.}
    \item{theta.1}{The angle of the ell.1 axis}
    \item{sigma.f}{(Scale parameter) The vertical feature lengthscale.}
  }
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[R:getLogspaceNames.CovarianceSEAniso2D]{getLogspaceNames}} \tab Names of "scale"-type parameters.\cr
 \tab \code{\link[R:getLowerPlain.CovarianceSEAniso2D]{getLowerPlain}} \tab Lower bounds for params, with plain names.\cr
 \tab \code{\link[R:getParamNamesPlain.CovarianceSEAniso2D]{getParamNamesPlain}} \tab Basenames of parameters.\cr
 \tab \code{\link[R:getParamsPlain.CovarianceSEAniso2D]{getParamsPlain}} \tab Parameter values with plain names.\cr
 \tab \code{\link[R:getUpperPlain.CovarianceSEAniso2D]{getUpperPlain}} \tab Upper bounds for params, with plain names.\cr
 \tab \code{\link[R:KDerivImplementation.CovarianceSEAniso2D]{KDerivImplementation}} \tab Element-wise derivatives of Covariance matrix.\cr
 \tab \code{\link[R:K.specific.CovarianceSEAniso2D]{K.specific}} \tab Anisotropic 2D SE Covariance matrix.\cr
 \tab \code{paramsPlainImplementation} \tab  -\cr
 \tab \code{setLowerPlain} \tab  -\cr
 \tab \code{setUpperPlain} \tab  -\cr
 \tab \code{\link[R:Variance.CovarianceSEAniso2D]{Variance}} \tab Anisotropic 2D SE variance at each point.\cr
}


 \bold{Methods inherited from Covariance}:\cr
clone, EncodeForTraining, FixConstParam, getId, getLower, getParamNames, getParams, getUpper, KInIn, KInInDeriv, KInOut, KOutIn, KOutOut, PrependId, print, setId, setLower, setParams, setParamsPlain, setUpper, clone, EncodeForTraining, FixConstParam, getId, getLower, getParamNames, getParams, getUpper, KInIn, KInInDeriv, KInOut, KOutIn, KOutOut, PrependId, print, setId, setLower, setParams, setParamsPlain, setUpper

 \bold{Methods inherited from Object}:\cr
as.character, attachLocally, attach, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, [[<-, [[, $<-, $, objectSize, print, registerFinalizer, save


}

\author{Charles R. Hogg III}
\keyword{classes}
