%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  R/Covariance.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{setLower.Covariance}
\alias{setLower.Covariance}
\alias{Covariance.setLower}
\alias{setLower.Covariance}
\alias{setLower,Covariance-method}


\title{Lower bounds on each parameter vector}

\description{
  Lower bounds on each parameter vector.
}

\usage{\method{setLower}{Covariance}(this, L, for.training=FALSE, ...)}

\arguments{
  \item{L}{A (named) vector of new lower bounds (we \emph{only} use ones
     which are named, and whose names match up with names of parameters.)}
  \item{for.training}{(logical) If TRUE, we return values more suitable for
     optimization: "scale"-type parameters are given as logarithms (which
     amounts to using the Jeffreys prior), and "constant" params (i.e., where
     the lower bound equals the upper bound) are omitted entirely.}
  \item{...}{Not used.}
}

\seealso{
  \code{\link[R:getLower.Covariance]{*getLower}()}
  \code{\link[R:setUpper.Covariance]{*setUpper}()}
  For more information see \code{\link{Covariance}}.
}

\author{Charles R. Hogg III}
\keyword{internal}
\keyword{methods}
