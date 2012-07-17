\name{Forget.Model}
\alias{Forget}
\alias{Forget.Model}
\title{Clear precomputed matrices from memory}
\arguments{
  \item{this}{The Model object.}

  \item{...}{Not used.}
}
\description{
  Model objects store computed matrices for easy reuse.
  This is usually very convenient, but sometimes we want
  the memory more than we want the speed! This method tells
  the Model to forget its results and free up some memory.
}
\seealso{
  \code{\link{Model}}

  \code{\link{LazyMatrix}}
}

