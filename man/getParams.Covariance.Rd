%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  R/Covariance.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getParams.Covariance}
\alias{getParams.Covariance}
\alias{Covariance.getParams}
\alias{getParams.Covariance}
\alias{getParams,Covariance-method}


\title{Parameter values for this Covariance}

\description{
  Parameter values for this Covariance.
}

\usage{\method{getParams}{Covariance}(this, for.training=FALSE, ...)}

\arguments{
  \item{for.training}{(logical) If TRUE, we return values more suitable for
     optimization: "scale"-type parameters are given as logarithms (which
     amounts to using the Jeffreys prior), and "constant" params (i.e., where
     the lower bound equals the upper bound) are omitted entirely.}
  \item{...}{Not used.}
}

\value{
  (named numeric vector) Parameter values for this Covariance.
}

\seealso{
  \code{\link[R:setParams.Covariance]{*setParams}()}
  \code{\link[R:getParamsPlain.Covariance]{*getParamsPlain}()}
  For more information see \code{\link{Covariance}}.
}

\alias{params.Covariance}

\author{Charles R. Hogg III}
\keyword{internal}
\keyword{methods}
