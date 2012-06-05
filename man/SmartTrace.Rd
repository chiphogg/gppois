\name{SmartTrace}
\alias{SmartTrace}
\title{Fast trace for matrix product}
\usage{
  SmartTrace(m1, m2)
}
\arguments{
  \item{m1}{One matrix.}

  \item{m2}{The other matrix.}
}
\value{
  Trace(\code{m1 %*% m2})
}
\description{
  Computes the trace of the matrix product \code{m1 %*%
  m2}, without actually evaluating that product.
}

