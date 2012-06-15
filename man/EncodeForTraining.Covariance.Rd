%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  R/Covariance.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{EncodeForTraining.Covariance}
\alias{EncodeForTraining.Covariance}
\alias{Covariance.EncodeForTraining}
\alias{EncodeForTraining.Covariance}
\alias{EncodeForTraining,Covariance-method}


\title{Makes parameters suitable for optimization}

\description{
  Makes parameters suitable for optimization.
  Specifically, we want "scale"-type parameters (i.e., where we're uncertain
  about the magnitude) to be uncertain in logspace rather than regular-space.
  This is equivalent to the Jeffreys prior.
}

\usage{\method{EncodeForTraining}{Covariance}(this, values, ...)}

\arguments{
  \item{this}{The Covariance object whose parameters to encode.}
  \item{values}{A named vector of parameter values (with ID already
     prepended).}
  \item{...}{Not used.}
}

\value{
  \code{values}, but with scale-type parameters changed to their log (and
  appropriately renamed)
}

\seealso{
  \code{\link{function}}() DecodeForTraining
  For more information see \code{\link{Covariance}}.
}

\author{Charles R. Hogg III}
\keyword{internal}
\keyword{methods}
