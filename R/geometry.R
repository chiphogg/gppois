
HexagonalShell <- function(n, base.point=c(0, 0), unit=1) {
  # Produces the n'th regular hexagonal shell, having 6*n points,
  # counterclockwise, in a coordinate system where rows are along X1 direction,
  # centred around the point [base.point+(n*unit, 0)].
  #
  # Args:
  #   n:  The number of the shell, i.e., the distance of this shell from the
  #      centre.
  #   base.point:  The centre of the shell.
  #   unit:  The spacing between neighboring points.
  #
  # Returns:
  #   A (6n x 2) matrix, whose 6n rows give the points in this shell, in
  #   counterclockwise order.

  # 'A' holds hexagonal 'basis' vectors a1, a2, a3.
  A <- unit * matrix(nrow=3, ncol=2,
    c(1, -0.5,        -0.5,
      0, sqrt(3) / 2, -sqrt(3) / 2))
  pts.in.shell <- matrix(NA, nrow=6 * n, ncol=2)
  for (b in 1:3) {
    for (j in 1:n) {
      idx <- 2 * (b - 1) * n + j
      neg.idx <- (idx + 3 * n - 1) %% (6 * n) + 1
      pt.offset <- n * A[b, ] + (j - 1) * A[(b %% 3) + 1, ]
      pts.in.shell[idx, ] <- base.point + pt.offset
      pts.in.shell[neg.idx, ] <- base.point - pt.offset
    }
  }
  return (pts.in.shell)
}

#' Hexagonal grid of points
#'
#' Constructs a hexagonal grid of points with \code{n} layers, spaced by
#' \code{unit}, and centred at \code{base.point}.
#'
#' @param n  The number of shells to construct
#' @param base.point  The centre of the shell
#' @param unit  The spacing between neighboring points
#'
#' @export
#' @return  A ((1+3n(n+1)) x 2) matrix, whose (1+3n(n+1)) rows give the points
#'    in this grid, in counterclockwise order and spiraling outwards.
HexagonalGrid <- function(n, base.point=c(0, 0), unit=1) {
  hex.grid.pts <- matrix(base.point, nrow=1, ncol=2)
  for (n.sh in 1:n) {
    hex.grid.pts <- rbind(deparse.level=0, hex.grid.pts,
      HexagonalShell(n=n.sh, base.point=base.point, unit=unit))
  }
  return (hex.grid.pts)
}

SortedConvexHullIndices <- function(X) {
  # Computes points comprising the convex hull for 'X', sorted so that they
  # occur in order (either clockwise or counterclockwise).
  #
  # Args:
  #   X:  The points whose convex hull we desire.
  #
  # Returns:
  #   Indices into 'X' corresponding to points comprising its convex hull ,
  #   sorted so that they occur in order (either clockwise or
  #   counterclockwise).
  if (!require("geometry")) {
    stop("Convex hull computations require the 'geometry' package.")
  }
  conv.hull <- convhulln(p=X + 0.0)
  vertices <- conv.hull[1, ]
  while(length(vertices) < nrow(conv.hull)) {
    curr.pt <- vertices[1]
    c.p.found.at <- which(conv.hull==curr.pt, arr.ind=TRUE)
    c.h.subset <- conv.hull[c.p.found.at[, "row"], ]
    next.pt.at <- which(!(c.h.subset %in% vertices))
    vertices <- c(c.h.subset[next.pt.at], vertices)
  }
  return (vertices)
}

InsideConvexHull <- function(point, sorted.hull) {
  # Checks whether a point falls within the convex hull 'sorted.hull'.
  #
  # Args:
  #   point:  A numeric vector of length 2 with the point's coordinates.
  #   sorted.hull:  A 2-column matrix with an ordered list of points in a
  #      convex hull.
  #
  # Returns:
  #   TRUE if 'point' is inside 'sorted.hull'; FALSE otherwise.
  P <- matrix(point, nrow=1)
  N <- nrow(sorted.hull)
  area <- 0
  for (idx in 1:(N + 1)) {
    i <- (idx %% N) + 1
    j <- ( i  %% N) + 1
    last.area <- area
    area <- det(rbind(deparse.level=0, sorted.hull[i, ] - P, sorted.hull[j, ] - P))
    if (last.area * area < 0) return (FALSE)
  }
  return (TRUE)
}

PartitionPointsUsingHull <- function(X, sorted.hull) {
  # Partition points in 'X' into two groups: those inside sorted.hull, and
  # those outside.
  #
  # Args:
  #   X:  2-column matrix holding the 2D points to partition.
  #   sorted.hull:  2-column matrix whose rows are points in a convex hull,
  #      traversed in order around the perimeter.
  #
  # Returns:
  #   A list with two elements, each of which are 2-column matrices whose rows
  #   denote 2D points: $inside is the points within the hull; $outside is all
  #   the rest.
  points.inside <- points.outside <- c()
  for (i in 1:NumPoints(X)) {
    if (InsideConvexHull(point=X[i, ], sorted.hull=sorted.hull)) {
      points.inside <- rbind(deparse.level=0, points.inside, X[i, ])
    } else {
      points.outside <- rbind(deparse.level=0, points.outside, X[i, ])
    }
  }
  return (list(inside=points.inside, outside=points.outside))
}

ClipPointsToHull <- function(X, sorted.hull) {
  # Remove points from 'X' which fall outside the region bounded by
  # 'sorted.hull'.
  #
  # Args:
  #   X:  2-column matrix holding the 2D points to partition.
  #   sorted.hull:  2-column matrix whose rows are points in a convex hull,
  #      traversed in order around the perimeter.
  #
  # Returns:
  #   A 2-column matrix with the points in 'X' which fall within 'sorted.hull'.
  return (PartitionPointsUsingHull(X=X, sorted.hull=sorted.hull)$inside)
}

#' Cover a convex region with a hexagonal grid
#'
#' Computes a hexagonal grid of the given \code{spacing} inside the convex hull
#' of the points in \code{X}.
#'
#' @param X  2-column matrix of points whose convex hull defines the region to
#'    grid.
#' @param spacing  The distance between points in the grid.
#'
#' @export
#' @return  A 2-column matrix of points in a hexagonal grid covering this
#'    region.
GriddedConvexHull <- function(X, spacing) {
  grid.centre <- matrix(ncol=2, c(mean(range(X[, 1])), mean(range(X[, 2]))))
  c.hull <- X[SortedConvexHullIndices(X), ]
  n.shells <- ceiling(max(DistanceMatrix(grid.centre, c.hull)) / (
      spacing * cos(pi / 6)))
  hex.grid <- HexagonalGrid(n=n.shells, base.point=grid.centre, unit=spacing)
  clipped.grid <- ClipPointsToHull(X=hex.grid, sorted.hull=c.hull)
  return (clipped.grid)
}

#' Flag points which are "in the gap"
#'
#' This function takes two sets of points: \code{X} is usually points where we
#' have data, and \code{X.out} is usually gridded points where we evaluate the
#' Gaussian Process.  It looks for points in \code{X.out} which are "gap"
#' points, i.e., far enough away from every point in \code{X}.  The threshhold
#' is \code{gap.thresh}, which defaults to the largest nearest-neighbor distance
#' in X.
#'
#' @param X  d-column matrix whose rows represent points where we have data.
#' @param X.out  d-column matrix whose rows represent points where we want to
#'      predict.
#' @param gap.thresh  The threshhold distance which determines \sQuote{gap}
#'    points, i.e., \sQuote{non-gap} points are closer than \code{gap.thresh} to
#'    any datapoint.
#' 
#' @export
#' @return  A logical matrix, of length NumPoints(X.out), TRUE if the
#'    corresponding point in X.out is a "gap point".
FindGapPoints <- function(X, X.out, gap.thresh=NA) {

  # Find the longest nearest-neighbor distance.
  if (!(is.numeric(gap.thresh) && gap.thresh > 0)) {
    DM <- as.matrix(dist(X))
    diag(DM) <- NA
    gap.thresh <- max(apply(DM, 2, min, na.rm=TRUE))
    rm(DM)
  }
  # Find the closest X-point for each X.out-point.
  X.out.closest <- apply(DistanceMatrix(X=X, X.out=X.out), 1, min)
  return (X.out.closest > gap.thresh)
}


