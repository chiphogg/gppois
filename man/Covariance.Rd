%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  R/Covariance.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{Covariance}
\docType{class}
\alias{Covariance}


\title{Superclass for all Covariance function objects}

\description{
  \code{Covariance} is the superclass for more specific types (SE, Matern,
  etc.).  Nobody will actually make an object of type \code{Covariance};
  instead, they will use one of the derived classes (such as
  \code{link{CovarianceSE}}).

  Superclass for all Covariance function objects.

  Package:   \cr
\bold{Class Covariance}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{Covariance}\cr

\bold{Directly known subclasses:}\cr
\link{CovarianceNoise}, \link{CovarianceSE}, \link{CovarianceSEAniso2D}, \link{CovarianceSELocalized}, \link{CovarianceSEVaryingEll}\cr

public static class \bold{Covariance}\cr
extends \link[R.oo]{Object}\cr


}

\usage{Covariance(id="", ...)}

\arguments{
  \item{id}{(character) A string to identify this covariance object.}
  \item{...}{Not used.}
}

\details{
   Regarding parameter names: the SUBCLASS has the responsibility to provide
   "plain-named" versions of all the virtual fields (paramNamesPlain,
   paramsPlain, lowerPlain, etc.).  The SUPERCLASS will automatically handle
   the "decorated" versions (paramNames, params, lower, etc.).

   All Covariance subclasses should remember the param-values and Dataset
   they last used to compute their K-matrix.  If they get asked to compute it
   again, they will simply return the previously-computed result if these
   values have not changed.
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{ClampedParamVals} \tab  -\cr
 \tab \code{ClampParams} \tab  -\cr
 \tab \code{\link[R:clone.Covariance]{clone}} \tab Deep-clone a Dataset.\cr
 \tab \code{EncodeForTraining} \tab  -\cr
 \tab \code{\link[R:FixConstParam.Covariance]{FixConstParam}} \tab Set a parameter to a constant value..\cr
 \tab \code{\link[R:getId.Covariance]{getId}} \tab The ID for this Covariance.\cr
 \tab \code{getLower} \tab  -\cr
 \tab \code{\link[R:getParamNames.Covariance]{getParamNames}} \tab ID-decorated names of Covariance parameters.\cr
 \tab \code{\link[R:getParams.Covariance]{getParams}} \tab Parameter values for this Covariance.\cr
 \tab \code{getUpper} \tab  -\cr
 \tab \code{KInIn} \tab  -\cr
 \tab \code{\link[R:KInInDeriv.Covariance]{KInInDeriv}} \tab Element-wise derivatives of covariance matrix.\cr
 \tab \code{KInOut} \tab  -\cr
 \tab \code{KOutIn} \tab  -\cr
 \tab \code{KOutOut} \tab  -\cr
 \tab \code{PrependId} \tab  -\cr
 \tab \code{\link[R:print.Covariance]{print}} \tab Pretty-printing for Covariance objects.\cr
 \tab \code{PushLowerBounds} \tab  -\cr
 \tab \code{PushUpperBounds} \tab  -\cr
 \tab \code{setId} \tab  -\cr
 \tab \code{\link[R:setLower.Covariance]{setLower}} \tab Lower bounds for parameters using full names.\cr
 \tab \code{setParams} \tab  -\cr
 \tab \code{\link[R:setParamsPlain.Covariance]{setParamsPlain}} \tab Set parameters using undecorated names.\cr
 \tab \code{\link[R:setUpper.Covariance]{setUpper}} \tab Upper bounds for parameters using full names.\cr
 \tab \code{UndecorateNames} \tab  -\cr
}


 \bold{Methods inherited from Object}:\cr
as.character, attachLocally, attach, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, [[<-, [[, $<-, $, objectSize, print, registerFinalizer, save


}

\author{Charles R. Hogg III}
\keyword{classes}
