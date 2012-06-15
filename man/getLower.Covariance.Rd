%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  R/Covariance.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getLower.Covariance}
\alias{getLower.Covariance}
\alias{Covariance.getLower}
\alias{getLower.Covariance}
\alias{getLower,Covariance-method}


\title{Lower bounds on each parameter vector}

\description{
  Lower bounds on each parameter vector.
}

\usage{\method{getLower}{Covariance}(this, for.training=FALSE, ...)}

\arguments{
  \item{for.training}{(logical) If TRUE, we return values more suitable for
     optimization: "scale"-type parameters are given as logarithms (which
     amounts to using the Jeffreys prior), and "constant" params (i.e., where
     the lower bound equals the upper bound) are omitted entirely.}
  \item{...}{Not used.}
}

\value{
  (named numeric vector) The lower bounds for the parameters for this
  covariance function.
}

\seealso{
  \code{\link[R:setLower.Covariance]{*setLower}()}
  \code{\link[R:getUpper.Covariance]{*getUpper}()}
  For more information see \code{\link{Covariance}}.
}

\author{Charles R. Hogg III}
\keyword{internal}
\keyword{methods}
