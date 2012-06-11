\name{Widths}
\alias{Widths}
\title{Computes the width of each point in X.}
\usage{
  Widths(X)
}
\arguments{
  \item{X}{A sorted numeric 1-D vector, representing points
  on a line.}
}
\value{
  The width of the "territory" for each point in X.
}
\description{
  For a sorted sequence of one-dimensional points X,
  calculates the width of each point's "territory" (i.e.,
  the set of points closer to that point than any other).
  The sum of all widths is the range of X.
}

