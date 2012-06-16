%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  R/Model.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{Model}
\docType{class}
\alias{Model}


\title{A trainable collection of Covariances}

\description{
  This is the \emph{constructor} for a \code{Model} object:
  A trainable collection of Covariances.

  The \code{Model} expresses our beliefs or knowledge about a Dataset.  We
  assume a Dataset can be modeled as the sum of one or more Gaussian Process
  Covariance functions, each of which is governed by parameters.  The Model
  can then be trained on a Dataset, a process which selects the parameter
  values that best describe the data.  The model can then be used to make
  predictions about the true function -- either at noisy datapoints, or
  interpolating into data-free regions, or both.

  This class should work just fine, as long as
    a) we are training on all the datapoints together (i.e., not breaking
       them up into subregions to divide-and-conquer), and
    b) this$params returns a vector which is amenable to simple optimization
       routines (i.e., none of the parameters require special treatment).
  If either of these conditions fail, a new approach is needed: either a
  specialized subclass should be created, or the problem should be broken
  into smaller pieces where these assumptions are good.

  Ironically, \emph{both} these conditions fail for \emph{both} scenarios
  considered in our Journal of Applied Crystallography paper, despite the
  fact that I wrote this software to perform the analysis for that paper.
  I hope to remedy this in a future version.  However, even in the meantime,
  having these classes still makes it very much easier to build the
  specialized functions I need.  Moreover, experience has shown that the
  plain-vanilla class structure is already good enough for other applications
  unrelated to denoising of scattering curves.

  Here is the class hierarchy:\cr
  Package:   \cr
\bold{Class Model}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{Model}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{Model}\cr
extends \link[R.oo]{Object}\cr



}

\usage{Model(id="", ...)}

\arguments{
  \item{id}{(character) An id which identifies this Model.}
  \item{...}{Not used.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[R:AddCovariance.Model]{AddCovariance}} \tab Add a new Covariance to this Model.\cr
 \tab \code{CheckContributionsAndWarn} \tab  -\cr
 \tab \code{\link[R:clone.Model]{clone}} \tab Deep-clone a Model.\cr
 \tab \code{ComputeKChol} \tab  -\cr
 \tab \code{ComputeL} \tab  -\cr
 \tab \code{\link[R:Forget.Model]{Forget}} \tab Clear precomputed matrices from memory.\cr
 \tab \code{\link[R:Freeze.Model]{Freeze}} \tab Make some parameters constant.\cr
 \tab \code{\link[R:getContributionIds.Model]{getContributionIds}} \tab Id's of each contributing Covariance.\cr
 \tab \code{\link[R:getId.Model]{getId}} \tab ID string for this Model.\cr
 \tab \code{\link[R:getLower.Model]{getLower}} \tab Lower bounds for parameters.\cr
 \tab \code{\link[R:getParams.Model]{getParams}} \tab Parameters for the Model.\cr
 \tab \code{\link[R:getSignalIds.Model]{getSignalIds}} \tab ID's of non-noise Covariances.\cr
 \tab \code{\link[R:getUpper.Model]{getUpper}} \tab Upper bounds for parameters.\cr
 \tab \code{getVaryingParamNames} \tab  -\cr
 \tab \code{KDeriv} \tab  -\cr
 \tab \code{KInv} \tab  -\cr
 \tab \code{K} \tab  -\cr
 \tab \code{KTotal} \tab  -\cr
 \tab \code{\link[R:L.Model]{L}} \tab Lower Cholesky root of covariance matrix.\cr
 \tab \code{LogDetK} \tab  -\cr
 \tab \code{\link[R:NamedCovariance.Model]{NamedCovariance}} \tab Retrieve one contributing Covariance.\cr
 \tab \code{\link[R:PlotBubblingSurfaces2D.Model]{PlotBubblingSurfaces2D}} \tab Animated uncertainty in a surface.\cr
 \tab \code{\link[R:PosteriorInterval.Model]{PosteriorInterval}} \tab Best estimate, including uncertainty.\cr
 \tab \code{\link[R:PosteriorMean.Model]{PosteriorMean}} \tab Best estimate of the true function.\cr
 \tab \code{\link[R:PosteriorStandardDeviation.Model]{PosteriorStandardDeviation}} \tab Pointwise uncertainty.\cr
 \tab \code{\link[R:PredictionMatrix.Model]{PredictionMatrix}} \tab Matrix connecting noisy data to true function.\cr
 \tab \code{\link[R:print.Model]{print}} \tab Pretty-printing for Model objects.\cr
 \tab \code{setId} \tab  -\cr
 \tab \code{setLower} \tab  -\cr
 \tab \code{\link[R:SetNoiseBounds.Model]{SetNoiseBounds}} \tab Uncertainty about the noise level.\cr
 \tab \code{setParams} \tab  -\cr
 \tab \code{setUpper} \tab  -\cr
 \tab \code{\link[R:Train.Model]{Train}} \tab Train a Model on a Dataset.\cr
 \tab \code{\link[R:AddCovariance.Model]{AddCovariance}} \tab Add a new Covariance to this Model.\cr
 \tab \code{CheckContributionsAndWarn} \tab  -\cr
 \tab \code{\link[R:clone.Model]{clone}} \tab Deep-clone a Model.\cr
 \tab \code{ComputeKChol} \tab  -\cr
 \tab \code{ComputeL} \tab  -\cr
 \tab \code{\link[R:Forget.Model]{Forget}} \tab Clear precomputed matrices from memory.\cr
 \tab \code{\link[R:Freeze.Model]{Freeze}} \tab Make some parameters constant.\cr
 \tab \code{\link[R:getContributionIds.Model]{getContributionIds}} \tab Id's of each contributing Covariance.\cr
 \tab \code{\link[R:getId.Model]{getId}} \tab ID string for this Model.\cr
 \tab \code{\link[R:getLower.Model]{getLower}} \tab Lower bounds for parameters.\cr
 \tab \code{\link[R:getParams.Model]{getParams}} \tab Parameters for the Model.\cr
 \tab \code{\link[R:getSignalIds.Model]{getSignalIds}} \tab ID's of non-noise Covariances.\cr
 \tab \code{\link[R:getUpper.Model]{getUpper}} \tab Upper bounds for parameters.\cr
 \tab \code{getVaryingParamNames} \tab  -\cr
 \tab \code{KDeriv} \tab  -\cr
 \tab \code{KInv} \tab  -\cr
 \tab \code{K} \tab  -\cr
 \tab \code{KTotal} \tab  -\cr
 \tab \code{\link[R:L.Model]{L}} \tab Lower Cholesky root of covariance matrix.\cr
 \tab \code{LogDetK} \tab  -\cr
 \tab \code{\link[R:NamedCovariance.Model]{NamedCovariance}} \tab Retrieve one contributing Covariance.\cr
 \tab \code{\link[R:PlotBubblingSurfaces2D.Model]{PlotBubblingSurfaces2D}} \tab Animated uncertainty in a surface.\cr
 \tab \code{\link[R:PosteriorInterval.Model]{PosteriorInterval}} \tab Best estimate, including uncertainty.\cr
 \tab \code{\link[R:PosteriorMean.Model]{PosteriorMean}} \tab Best estimate of the true function.\cr
 \tab \code{\link[R:PosteriorStandardDeviation.Model]{PosteriorStandardDeviation}} \tab Pointwise uncertainty.\cr
 \tab \code{\link[R:PredictionMatrix.Model]{PredictionMatrix}} \tab Matrix connecting noisy data to true function.\cr
 \tab \code{\link[R:print.Model]{print}} \tab Pretty-printing for Model objects.\cr
 \tab \code{setId} \tab  -\cr
 \tab \code{setLower} \tab  -\cr
 \tab \code{\link[R:SetNoiseBounds.Model]{SetNoiseBounds}} \tab Uncertainty about the noise level.\cr
 \tab \code{setParams} \tab  -\cr
 \tab \code{setUpper} \tab  -\cr
 \tab \code{\link[R:Train.Model]{Train}} \tab Train a Model on a Dataset.\cr
}


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
