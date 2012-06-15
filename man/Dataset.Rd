%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  R/Dataset.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{Dataset}
\docType{class}
\alias{Dataset}


\title{A wrapper class for data being analyzed}

\description{
  This is the \emph{constructor} for a \code{Dataset} object:
  A wrapper class for data being analyzed.

  Here is the class hierarchy:
  Package:   \cr
\bold{Class Dataset}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{Dataset}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{Dataset}\cr
extends \link[R.oo]{Object}\cr



}

\usage{Dataset(id="UNNAMED", data=data.frame(), X.names="X", noise.var=c(), column="", poisson.names=c(), tol.factor=1e-04, data.offset=c(), ...)}

\arguments{
  \item{id}{(character) An id which identifies this Dataset.}
  \item{data}{(\code{data.frame}) The raw data.}
  \item{X.names}{(character vector) The names of the columns which represent
     covariates.  (NB: the length of this vector implicitly defines the
     dimensionality of the data: if you give 2 column names, you have a 2D
     dataset.)}
  \item{noise.var}{(named numeric vector) Desired noise levels for each
     column (i.e., a threshhold for precision: we don't care about
     differences smaller than the corresponding scale).}
  \item{column}{(character) The name of the quantity to select for analysis.}
  \item{poisson.names}{(character vector) The column names which represent
     Poisson data.}
  \item{tol.factor}{(numeric) The relative tolerance factor: by default, we
     set the noise variance for each column to \code{tol.factor * sd(col)}.}
  \item{data.offset}{(named numeric vector) An optional offset which is
     subtracted from each column.  If an unnamed value is included, it is
     treated as a default.  The "default default" is to subtract off the mean
     for each column.}
  \item{...}{Not used.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[R:DeleteRows.Dataset]{DeleteRows}} \tab Delete datapoints.\cr
 \tab \code{\link[R:getDataOffset.Dataset]{getDataOffset}} \tab The offset for this dataset.\cr
 \tab \code{\link[R:getD.Dataset]{getD}} \tab The dimensionality of the Dataset.\cr
 \tab \code{\link[R:getDpts.Dataset]{getDpts}} \tab Datapoint values.\cr
 \tab \code{getId} \tab  -\cr
 \tab \code{\link[R:getIsPoisson.Dataset]{getIsPoisson}} \tab Is the current column Poisson?.\cr
 \tab \code{\link[R:getN.Dataset]{getN}} \tab The number of datapoints.\cr
 \tab \code{\link[R:getNoiseVar.Dataset]{getNoiseVar}} \tab Noise variance for current column.\cr
 \tab \code{getQuantity} \tab  -\cr
 \tab \code{\link[R:getX.Dataset]{getX}} \tab X-values.\cr
 \tab \code{\link[R:getXformedDpts.Dataset]{getXformedDpts}} \tab Transformed datapoint values.\cr
 \tab \code{\link[R:MSR.Dataset]{MSR}} \tab Mean square residuals.\cr
 \tab \code{\link[R:Plot2D.Dataset]{Plot2D}} \tab Scatterplot for 2D dataset.\cr
 \tab \code{\link[R:print.Dataset]{print}} \tab Pretty-printing for Dataset objects.\cr
 \tab \code{\link[R:RemoveRange.Dataset]{RemoveRange}} \tab Remove datapoints within a given range.\cr
 \tab \code{\link[R:Same.Dataset]{Same}} \tab Check whether two Datasets are identical.\cr
 \tab \code{setId} \tab  -\cr
 \tab \code{setQuantity} \tab  -\cr
 \tab \code{\link[R:Untransform.Dataset]{Untransform}} \tab Undo any transformations on the datapoints.\cr
 \tab \code{\link[R:DeleteRows.Dataset]{DeleteRows}} \tab Delete datapoints.\cr
 \tab \code{\link[R:getDataOffset.Dataset]{getDataOffset}} \tab The offset for this dataset.\cr
 \tab \code{\link[R:getD.Dataset]{getD}} \tab The dimensionality of the Dataset.\cr
 \tab \code{\link[R:getDpts.Dataset]{getDpts}} \tab Datapoint values.\cr
 \tab \code{getId} \tab  -\cr
 \tab \code{\link[R:getIsPoisson.Dataset]{getIsPoisson}} \tab Is the current column Poisson?.\cr
 \tab \code{\link[R:getN.Dataset]{getN}} \tab The number of datapoints.\cr
 \tab \code{\link[R:getNoiseVar.Dataset]{getNoiseVar}} \tab Noise variance for current column.\cr
 \tab \code{getQuantity} \tab  -\cr
 \tab \code{\link[R:getX.Dataset]{getX}} \tab X-values.\cr
 \tab \code{\link[R:getXformedDpts.Dataset]{getXformedDpts}} \tab Transformed datapoint values.\cr
 \tab \code{\link[R:MSR.Dataset]{MSR}} \tab Mean square residuals.\cr
 \tab \code{\link[R:Plot2D.Dataset]{Plot2D}} \tab Scatterplot for 2D dataset.\cr
 \tab \code{\link[R:print.Dataset]{print}} \tab Pretty-printing for Dataset objects.\cr
 \tab \code{\link[R:RemoveRange.Dataset]{RemoveRange}} \tab Remove datapoints within a given range.\cr
 \tab \code{\link[R:Same.Dataset]{Same}} \tab Check whether two Datasets are identical.\cr
 \tab \code{setId} \tab  -\cr
 \tab \code{setQuantity} \tab  -\cr
 \tab \code{\link[R:Untransform.Dataset]{Untransform}} \tab Undo any transformations on the datapoints.\cr
}


 \bold{Methods inherited from Object}:\cr
as.character, attachLocally, attach, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, [[<-, [[, $<-, $, objectSize, print, registerFinalizer, save


}

\author{Charles R. Hogg III}
\keyword{classes}