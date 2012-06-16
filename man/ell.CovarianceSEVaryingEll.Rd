\name{ell.CovarianceSEVaryingEll}
\alias{ell.CovarianceSEVaryingEll}
\title{ell(X)}
\arguments{
  \item{X}{The X-values where we want to know
  \code{ell(X)}.}

  \item{...}{Not used.}
}
\value{
  \code{ell} evaluated at each point in \code{X}.
}
\description{
  Evaluate ell(X) at the applied X-points (by
  spline-interpolating)
}
\note{
  The spline interpolation takes place in log-space
  (\code{ell} can never be negative), but a real-space
  result is returned.
}
\seealso{
  \code{\link{CovarianceSEVaryingEll}}
}

