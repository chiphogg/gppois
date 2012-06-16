%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  R/LazyMatrix.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{LazyMatrix}
\docType{class}
\alias{LazyMatrix}


\title{Wrapper to avoid recomputing matrices}

\description{
  This is the \emph{constructor} for a \code{LazyMatrix} object:
  Wrapper to avoid recomputing matrices.

  Note that this class doesn't do any actual computation!  It just stores the
  results, and also tells other code whether or not it needs to.

  Here is the class hierarchy:
  Package:   \cr
\bold{Class LazyMatrix}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{LazyMatrix}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{LazyMatrix}\cr
extends \link[R.oo]{Object}\cr



}

\usage{LazyMatrix(...)}

\arguments{
  \item{...}{Not used.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{getM} \tab  -\cr
 \tab \code{NeedToRecalculate} \tab  -\cr
 \tab \code{StoreMatrix} \tab  -\cr
}


 \bold{Methods inherited from Object}:\cr
as.character, attachLocally, attach, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, [[<-, [[, $<-, $, objectSize, print, registerFinalizer, save


}

\author{Charles R. Hogg III}
\keyword{classes}
