\name{MSR.Dataset}
\alias{MSR.Dataset}
\title{Mean square residuals}
\arguments{
  \item{test.data}{A sequence of datapoints, which live in
  the same space as xformedDpts (i.e., same offset, and
  Anscombe-transformed if Poisson).}

  \item{reference.data}{The comparison function: defaults
  to our noisy datapoints, but could also be different
  (e.g., if we have access to the true function).}

  \item{...}{Not used.}
}
\value{
  The mean squared difference between \code{test.data} and
  \code{reference.data}.
}
\description{
  Mean squared residuals are taken with respect to the
  currently selected column by default.
}
\seealso{
  \code{\link{Dataset}}
}

