\name{RemoveRange.Dataset}
\alias{RemoveRange.Dataset}
\title{Remove datapoints within a given range}
\arguments{
  \item{X.min}{The left boundary of the X-range to remove.}

  \item{X.max}{The right boundary of the X-range to
  remove.}

  \item{...}{Not used.}
}
\description{
  Remove datapoints within a given range, for 1D data only.
  This function could be improved by building an engine
  that works on arbitrary dimensions, maybe by inputting
  some kind of generalized metric, and using L2-norm for
  spheres, L-oo norm for boxes, etc.  But for now: 1D only!
}
\seealso{
  \code{\link{Dataset}}
}

