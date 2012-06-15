%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  R/Covariance.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getUpper.Covariance}
\alias{getUpper.Covariance}
\alias{Covariance.getUpper}
\alias{getUpper.Covariance}
\alias{getUpper,Covariance-method}


\title{Upper bounds for the parameter values}

\description{
  Upper bounds for the parameter values.
}

\usage{\method{getUpper}{Covariance}(this, for.training=FALSE, ...)}

\arguments{
  \item{for.training}{(logical) If TRUE, we return values more suitable for
     optimization: "scale"-type parameters are given as logarithms (which
     amounts to using the Jeffreys prior), and "constant" params (i.e., where
     the lower bound equals the upper bound) are omitted entirely.}
  \item{...}{Not used.}
}

\value{
  (named numeric vector) The upper bounds for the parameters for this
  covariance function.
}

\seealso{
  \code{\link[R:getLower.Covariance]{*getLower}()}
  \code{\link[R:setUpper.Covariance]{*setUpper}()}
  For more information see \code{\link{Covariance}}.
}

\author{Charles R. Hogg III}
\keyword{internal}
\keyword{methods}
