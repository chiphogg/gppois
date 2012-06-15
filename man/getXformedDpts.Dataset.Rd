\name{getXformedDpts.Dataset}
\alias{getXformedDpts}
\alias{getXformedDpts.Dataset}
\alias{xformedDpts}
\title{Transformed datapoint values}
\arguments{
  \item{...}{Not used.}
}
\value{
  The datapoint values for the currently-selected quantity,
  transformed such that they can be assumed to have a
  Gaussian distribution.
}
\description{
  Transformed values for the currently-selected quantity.
  Transformations include the Anscombe transform (for
  Poisson data), followed by subtracting off the mean (if
  desired).
}
\details{
  Usually accesed using \code{this$xformedDpts} syntax;
  read-only.
}
\seealso{
  \code{\link{Dataset}}

  \code{\link{dpts}}
}

