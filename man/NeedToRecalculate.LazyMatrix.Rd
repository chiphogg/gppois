\name{NeedToRecalculate.LazyMatrix}
\alias{NeedToRecalculate}
\alias{NeedToRecalculate.LazyMatrix}
\title{Decides whether we need to recompute}
\arguments{
  \item{...}{Not used.}
}
\value{
  TRUE iff we will need to recompute this matrix because
  the ingredients are different.
}
\description{
  Decides whether or not we would need to recompute the
  matrix, by checking whether the supplied "ingredients"
  are equal to the ones previously supplied.
}
\seealso{
  \code{\link{LazyMatrix}}
}

