\name{sigma.f.CovarianceSEVaryingEll}
\alias{sigma.f}
\alias{sigma.f.CovarianceSEVaryingEll}
\title{sigma.f(X)}
\arguments{
  \item{X}{The X-values where we want to know
  \code{sigma.f(X)}.}

  \item{...}{Not used.}
}
\value{
  \code{sigma.f} evaluated at each point in \code{X}.
}
\description{
  Evaluate sigma.f(X) at the applied X-points (by
  spline-interpolating)
}
\note{
  The spline interpolation takes place in log-space
  (\code{sigma.f} can never be negative), but a real-space
  result is returned.
}
\seealso{
  \code{\link{CovarianceSEVaryingEll}}
}

