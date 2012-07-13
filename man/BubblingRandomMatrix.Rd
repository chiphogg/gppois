\name{BubblingRandomMatrix}
\alias{BubblingRandomMatrix}
\title{Generate marginally normal smooth random timetraces}
\usage{
  BubblingRandomMatrix(n.pts, N = 10, n.times = 100, ...)
}
\arguments{
  \item{n.pts}{The number of independent oscillators.}

  \item{N}{One-half the number of independent draws to take
  for each oscillator.}

  \item{n.times}{The final number of interpolated points.}

  \item{...}{Not used.}
}
\value{
  A numeric matrix with n rows and n.times columns, showing
  the timetraces of the oscillators.
}
\description{
  This function is an "engine" for visualizing uncertainty
  using animations -- at least, as long as that uncertainty
  is normally distributed.  The output is a matrix whose
  rows represent independent oscillators, and whose columns
  represent timesteps.  Every column is i.i.d. normal.
  Every row has a CDF which approaches the normal
  distribution for large N, but which is nevertheless
  smooth to all orders of time derivatives.  Increasing N
  improves the statistical properties, while increasing
  n.times relative to N yields smoother animations.
}

