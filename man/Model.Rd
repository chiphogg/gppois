\name{Model}
\alias{Model}
\title{Model: a trainable collection of Covariances}
\arguments{
  \item{id}{(character) An id which identifies this Model.}

  \item{...}{Not used.}
}
\description{
  The \code{Model} expresses our beliefs or knowledge about
  a \code{\link{Dataset}}.  We assume a
  \code{\link{Dataset}} can be modeled as the sum of one or
  more Gaussian Process Covariance functions, each of which
  is governed by parameters.  The \code{Model} can then be
  trained on a \code{\link{Dataset}}, a process which
  selects the parameter values that best describe the data.
  The Model can then be used to make predictions about the
  true function -- either at noisy datapoints, or
  interpolating into data-free regions, or both.
}
\details{
  This class should work just fine, as long as a) we are
  training on all the datapoints together (i.e., not
  breaking them up into subregions to divide-and-conquer),
  and b) this$params returns a vector which is amenable to
  simple optimization routines (i.e., none of the
  parameters require special treatment).  If either of
  these conditions fail, a new approach is needed: either a
  specialized subclass should be created, or the problem
  should be broken into smaller pieces where these
  assumptions are good.

  Ironically, \emph{both} these conditions fail for
  \emph{both} scenarios considered in our Journal of
  Applied Crystallography paper, despite the fact that I
  wrote this software to perform the analysis for that
  paper.  I hope to remedy this in a future version.
  However, even in the meantime, having these classes still
  makes it very much easier to build the specialized
  functions I need.  Moreover, experience has shown that
  the \sQuote{plain-vanilla} class structure is already
  good enough for other applications unrelated to denoising
  of scattering curves.
}
\references{
  Hogg, C., K. Mullen, and I. Levin (2012). A Bayesian
  approach for denoising one-dimensional data. Journal of
  Applied Crystallography, 45(3), pp. 471-481.
}

