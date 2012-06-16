\name{LogML}
\alias{LogML}
\title{Log of the Marginal Likelihood}
\usage{
  LogML(par = model$getParams(for.training = TRUE), model,
    d, update.params = TRUE)
}
\arguments{
  \item{par}{The parameter values to use for calculating.}

  \item{model}{The Model object we're optimizing.}

  \item{d}{The Dataset we're training on.}

  \item{update.params}{logical; if TRUE, we should change
  model's params to the values in par.}
}
\value{
  The logarithm of the marginal likelihood; also, has the
  side-effect that \code{moedel$params} are set to
  \code{par} unless \code{update.params=FALSE}.
}
\description{
  Sets the parameter values to \code{par} for Model
  \code{model}, and returns the log of the (M)arginal
  (L)ikelihood for describing Dataset \code{d}.
}
\references{
  Rasmussen, C.E. and C.K.I. Williams (2005.) Gaussian
  Processes for Machine Learning. The MIT Press.
  \url{http://www.gaussianprocess.org/gpml/}
}

