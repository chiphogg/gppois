\name{Same.Dataset}
\alias{Same}
\alias{Same.Dataset}
\title{Check whether two Datasets are identical}
\arguments{
  \item{d}{The Dataset for comparison.}

  \item{compare}{The attributes (e.g., X, noiseVar, etc.)
  to check for equality.}

  \item{...}{Not used.}
}
\value{
  (logical) TRUE if all the requested attributes are the
  same.
}
\description{
  Check whether two Datasets have identical values for some
  subset of their parameters.  Many different objects in
  \code{gppois} perform expensive calculations on Dataset
  objects, then cache the results.  These objects need to
  be able to check that nothing in a Dataset has changed...
  at least, nothing that affects their huge calculation!
  This lets these objects return the cached result, saving
  considerable time.
}
\seealso{
  \code{\link{Dataset}}
}

