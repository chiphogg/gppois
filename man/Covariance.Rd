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
\cr

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
 \tab \code{clone} \tab  -\cr
 \tab \code{EncodeForTraining} \tab  -\cr
 \tab \code{FixConstParam} \tab  -\cr
 \tab \code{getId} \tab  -\cr
 \tab \code{getLower} \tab  -\cr
 \tab \code{getParamNames} \tab  -\cr
 \tab \code{getParams} \tab  -\cr
 \tab \code{getUpper} \tab  -\cr
 \tab \code{KInIn} \tab  -\cr
 \tab \code{KInInDeriv} \tab  -\cr
 \tab \code{KInOut} \tab  -\cr
 \tab \code{KOutIn} \tab  -\cr
 \tab \code{KOutOut} \tab  -\cr
 \tab \code{PrependId} \tab  -\cr
 \tab \code{print} \tab  -\cr
 \tab \code{PushLowerBounds} \tab  -\cr
 \tab \code{PushUpperBounds} \tab  -\cr
 \tab \code{setId} \tab  -\cr
 \tab \code{setLower} \tab  -\cr
 \tab \code{setParams} \tab  -\cr
 \tab \code{setParamsPlain} \tab  -\cr
 \tab \code{setUpper} \tab  -\cr
 \tab \code{UndecorateNames} \tab  -\cr
}


 \bold{Methods inherited from Object}:\cr
as.character, attachLocally, attach, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, [[<-, [[, $<-, $, objectSize, print, registerFinalizer, save


}

\author{Charles R. Hogg III}
\keyword{classes}
