%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  R/Covariance.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{setParams.Covariance}
\alias{setParams.Covariance}
\alias{Covariance.setParams}
\alias{setParams.Covariance}
\alias{setParams,Covariance-method}


\title{Parameter values for this Covariance}

\description{
  Parameter values for this Covariance.
}

\usage{\method{setParams}{Covariance}(this, p, for.training=FALSE, ...)}

\arguments{
  \item{p}{(named numeric vector) New parameter values. (N.B: we \emph{only}
     use ones which are named, and whose names match up with names of
     parameters.)}
  \item{for.training}{(logical) If TRUE, we return values more suitable for
     optimization: "scale"-type parameters are given as logarithms (which
     amounts to using the Jeffreys prior), and "constant" params (i.e., where
     the lower bound equals the upper bound) are omitted entirely.}
  \item{...}{Not used.}
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
