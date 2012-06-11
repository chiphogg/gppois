\name{Clamp}
\alias{Clamp}
\title{Clamps \code{x} to within \code{bounds}}
\usage{
  Clamp(x, bounds)
}
\arguments{
  \item{x}{Some numeric object whose values must be
  clamped.}

  \item{bounds}{A numeric object whose range sets the range
  for clamping.}
}
\value{
  \code{x}, but with values clamped according to
  \code{bounds}.
}
\description{
  Clamps every value in \code{x} (numeric) to lie within
  the range of \code{bounds} (numeric).
}

