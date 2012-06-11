\name{HexagonalGrid}
\alias{HexagonalGrid}
\title{Hexagonal grid of points}
\usage{
  HexagonalGrid(n, base.point = c(0, 0), unit = 1)
}
\arguments{
  \item{n}{The number of shells to construct}

  \item{base.point}{The centre of the shell}

  \item{unit}{The spacing between neighboring points}
}
\value{
  A ((1+3n(n+1)) x 2) matrix, whose (1+3n(n+1)) rows give
  the points in this grid, in counterclockwise order and
  spiraling outwards.
}
\description{
  Constructs a hexagonal grid of points with \code{n}
  layers, spaced by \code{unit}, and centred at
  \code{base.point}.
}

