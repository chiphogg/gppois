\name{Untransform.Dataset}
\alias{Untransform.Dataset}
\title{Undo any transformations on the datapoints}
\arguments{
  \item{values}{A numeric vector whose values "live" in the
  \emph{transformed} scale.}

  \item{...}{Not used.}
}
\value{
  Numbers corresponding to \code{values}, but in a
  non-transformed scale (i.e., instead of being comparable
  to \code{xformedDpts}, the returned quantities are
  comparable to \code{dpts}).
}
\description{
  Datapoints are often transformed prior to being analyzed.
  Examples include\cr * the Anscombe transform for
  Poisson-noised data, which makes it approximately
  Gaussian with constant variance\cr * subtracting off the
  mean\cr This function provides a (statistical!) inverse
  for these transforms.  It is intended to be applied to
  the results of the analysis (which take place in the
  \emph{transformed} space).
}
\seealso{
  \code{\link{Dataset}}

  \code{\link{xformedDpts}}

  \code{\link{dpts}}
}

