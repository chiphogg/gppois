\name{Train.Model}
\alias{Train}
\alias{Train.Model}
\title{Train a Model on a Dataset}
\arguments{
  \item{this}{The Model object.}

  \item{d}{(Dataset) The data which our parameters should
  describe.}

  \item{...}{Not used.}
}
\description{
  Optimize this Model's parameters so they describe the
  given data.
}
\details{
  Presently, we use R's \code{\link{optim}} \dQuote{under
  the hood}.  Specifically, we use the L-BFGS-B approach
  because we found it gave good results, and because it
  naturally incorporates boundaries.

  Future versions may have an option to calculate the
  Hessian analytically, which would allow trust region
  approaches to be used.  This should speed convergence.
}
\seealso{
  \code{\link{Model}}
}

