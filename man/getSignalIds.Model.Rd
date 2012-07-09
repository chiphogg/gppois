\name{getSignalIds.Model}
\alias{Model$signalIds}
\alias{getSignalIds.Model}
\title{ID's of non-noise Covariances}
\arguments{
  \item{this}{The Model object.}

  \item{...}{Not used.}
}
\value{
  The id's for this Model's non-noise contributing
  Covariance objects.
}
\description{
  A \code{Model} is a collection of \code{Covariances}.
  Some subset of these will correspond to the signal of
  actual interest; the rest are considered noise.  This
  function returns the list of contributions which are
  considered to be \dQuote{signal}.
}
\seealso{
  \code{\link{Model}}
}

