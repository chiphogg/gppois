\name{getIsPoisson.Dataset}
\alias{getIsPoisson}
\alias{getIsPoisson.Dataset}
\alias{isPoisson}
\title{Is the current column Poisson?}
\arguments{
  \item{...}{Not used.}
}
\value{
  (logical) TRUE if the raw data for the currently selected
  column is Poisson.
}
\description{
  Checks whether the currently-selected column consists of
  Poisson-noised data.
}
\details{
  Usually accesed using \code{this$isPoisson} syntax;
  read-only.
}
\seealso{
  \code{\link{Dataset}}
}

