\name{GriddedConvexHull}
\alias{GriddedConvexHull}
\title{Cover a convex region with a hexagonal grid}
\usage{
  GriddedConvexHull(X, spacing)
}
\arguments{
  \item{X}{2-column matrix of points whose convex hull
  defines the region to grid.}

  \item{spacing}{The distance between points in the grid.}
}
\value{
  A 2-column matrix of points in a hexagonal grid covering
  this region.
}
\description{
  Computes a hexagonal grid of the given \code{spacing}
  inside the convex hull of the points in \code{X}.
}

