%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  R/Covariance.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{setUpper.Covariance}
\alias{setUpper.Covariance}
\alias{Covariance.setUpper}
\alias{setUpper.Covariance}
\alias{setUpper,Covariance-method}


\title{Upper bounds for the parameter values}

\description{
  Upper bounds for the parameter values.
}

\usage{\method{setUpper}{Covariance}(this, U, for.training=FALSE, ...)}

\arguments{
  \item{U}{A (named) vector of new upper bounds (we \emph{only} use ones
     which are named, and whose names match up with names of parameters.)}
  \item{for.training}{(logical) If TRUE, we return values more suitable for
     optimization: "scale"-type parameters are given as logarithms (which
     amounts to using the Jeffreys prior), and "constant" params (i.e., where
     the lower bound equals the upper bound) are omitted entirely.}
  \item{...}{Not used.}
}

\seealso{
  \code{\link[R:setLower.Covariance]{*setLower}()}
  \code{\link[R:getUpper.Covariance]{*getUpper}()}
  For more information see \code{\link{Covariance}}.
}

\author{Charles R. Hogg III}
\keyword{internal}
\keyword{methods}
