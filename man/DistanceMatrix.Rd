\name{DistanceMatrix}
\alias{DistanceMatrix}
\title{Pairwise distances between points in X and X.out}
\usage{
  DistanceMatrix(X, X.out = X)
}
\arguments{
  \item{X}{A numeric matrix of input points.}

  \item{X.out}{A matrix of output points, whose distance to
  every point in 'X' is desired.}
}
\value{
  A matrix whose [i, j] component gives the Euclidean
  distance from X.out[i, ] to X[j, ].
}
\description{
  Computes the distance from every point in X to every
  point in X.out.  Both arguments are assumed to be numeric
  matrices with as many columns as the dimensionality of
  the space.  (i.e., N 2D points would be represented by an
  (N x 2) matrix, etc.)  Vector arguments are assumed to be
  1D points, and are automatically converted to matrices.
}

