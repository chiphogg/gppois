\name{Anscombe}
\alias{Anscombe}
\alias{AnscombeInverse}
\title{Anscombe transform for Poisson-noised data}
\usage{
  Anscombe(Y.p)

  AnscombeInverse(Y.g)
}
\arguments{
  \item{Y.p}{Data (presumably Poisson-noised) to be
  transformed}

  \item{Y.g}{Data (presumably Anscombe-transformed) to be
  untransformed to the original Poisson space.}
}
\value{
  The transformed (\code{Anscombe}) or untransformed
  (\code{AnscombeInverse}) data
}
\description{
  The Anscombe transform takes i.i.d. Poisson-noised data
  into i.i.d. normal with constant variance (to a very good
  approximation).  Gaussian-noised data is much easier to
  fit with Gaussian Processes.  The inverse Anscombe
  transform takes functions from this transformed space
  back into the original space.
}
\references{
  Anscombe, F.J. (1948). The Transformation of Poisson,
  Binomial and Negative-Binomial Data. Biometrika 35, pp.
  246-254.
}

