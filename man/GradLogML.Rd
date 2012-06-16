\name{GradLogML}
\alias{GradLogML}
\title{Gradient of log marginal likelihood}
\usage{
  GradLogML(par = model$getParams(for.training = TRUE),
    model, d, update.params = TRUE)
}
\arguments{
  \item{par}{The parameter values to use for calculating.}

  \item{model}{The Model object we're optimizing.}

  \item{d}{The Dataset we're training on.}

  \item{update.params}{logical; if TRUE, we should change
  model's params to the given value.}
}
\value{
  The gradient of the log of the marginal likelihood (also
  has a side-effect of setting \code{model$params <- par}).
}
\description{
  Sets the parameter values to \code{par} for Model
  \code{model}, and returns the gradient (w.r.t. par) of
  the log of the (M)arginal (L)ikelihood for describing
  Dataset \code{d}.  The gradient of LogML is very helpful
  in optimizing the parameter values.
}
\references{
  Rasmussen, C.E. and C.K.I. Williams (2005.) Gaussian
  Processes for Machine Learning. The MIT Press.
  \url{http://www.gaussianprocess.org/gpml/}
}

