# GPPois.R: Gaussian Process-based inference for Poisson-noised data.
# Copyright (c) 2011 Charles R. Hogg III
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Author: Charles R. Hogg III (2011) <charles.r.hogg@gmail.com>

# FILE DESCRIPTION:
# Provides classes and functions to analyze data using Gaussian Processes.  Key
# examples include:
#  - Dataset: the data to be analyzed
#  - Covariance: (virtual) superclass for the different possible kinds of
#      covariance structures (e.g. Squared-Exponential, periodic, Matern, ...)
#  - Model: an interface to the total covariance structure (different
#      Covariance terms contribute *additively*)
#
# These objects *need* to be mutable, so I am using the ''R.oo'' library to
# enable pass-by-reference.
#
# A lot of computations result in really big matrices.  These should be easy to
# save and restore to disk.
#
# Also, I am trying to adhere to Google's R Style Guide as much as possible:
# http://google-styleguide.googlecode.com/svn/trunk/google-r-style.html

# Object-orientation (with pass-by-reference!!)
library("R.oo")

################################################################################
# SECTION:                      Helper Functions
################################################################################

DebugIfError <- function(FUN, ...) {
  # Try running FUN(...) normally, but if there's an error, turn on debugging,
  # run it again, and turn debugging back off for this function.
  #
  # FUN:  The function to call (either a character string or just the name of
  #    the function).
  # ...:  The argument list for FUN.
  #
  # Returns:
  #   The result of calling FUN(...).
  debug.exists <- !isTRUE(
    suppressMessages(suppressWarnings(require("debug"))) == FALSE)
  if (is.character(FUN)) {
    FUN <- get(FUN)
  }
  arg.list <- list(...)
  value <- tryCatch(
    expr=FUN(...),
    error=function(e) {
      save(arg.list, file="exact_arguments.rda")
      # Make sure we can load the debug package
      if (debug.exists) {
        cat(sprintf("Error encountered! '%s'\nTrying again in debug mode.\n",
            e$message))
        mtrace(FUN)
        value <- FUN(...)
        mtrace(FUN, FALSE)
      } else {
        value <- e
      }
      return (value)
    }, finally={
      if (debug.exists) {
        mtrace(FUN, FALSE)
      }
    })
  return (value)
}

SmartTrace <- function(m1, m2) {
  # Computes the trace of the matrix product (m1 %*% m2), without actually
  # evaluating that product.
  #
  # Args:
  #   m1: One matrix.
  #   m2: The other matrix.
  #
  # Returns:
  #   Tr(m1 %*% m2).
  return (sum(as.vector(m1) * as.vector(t(m2))))
}

Clamp <- function(x, bounds) {
  # Clamps every value in 'x' to lie within the range of 'bounds'.
  #
  # Args:
  #   x: Some numeric object whose values must be clamped.
  #   bounds: A numeric object whose range sets the range for clamping.
  #
  # Returns:
  #   'x', but with values clamped according to 'bounds'.
  return (pmax(pmin(x, max(bounds)), min(bounds)))
}

ClampNamed <- function(x, lower, upper) {
  # Clamps every value in 'x' to lie between the values of 'lower' and of
  # 'upper' with the same name.
  #
  # Args:
  #   x: A named numeric vector to be clamped.
  #   lower: A numeric vector holding the lower bounds.
  #   upper: A numeric vector holding the upper bounds.
  #
  # Returns:
  #   'x', but with values clamped between 'lower' and 'upper'.
  n.x <- names(x)
  return (pmax(pmin(x[n.x], upper[n.x]), lower[n.x]))
}

ItTakes <- function(my.task, how.i.do.it) {
  # Measure how long it takes to do a given task, and print it out in
  # human-readable format.
  #
  # Args:
  #   my.task:  A short label describing what I'm trying to do
  #   how.i.do.it:  A code block; the implementation of my.task.
  #
  # Returns:
  #   Used for its side-effect.
  p.env <- parent.frame()
  cat(sprintf("%60s: ", my.task))
  elapsed <- system.time(eval(substitute(how.i.do.it), env=p.env))
  cat(sprintf("%12g sec\n", elapsed["user.self"]))
}

PlotMatrixQuickAndDirty <- function(M, ...) {
  # Plots a matrix using base R graphics.
  #
  # Args:
  #   M:  The matrix to plot.
  #
  # Returns:
  #   Used for its side-effect.
  #
  # Adapted from:
  # http://www.phaget4.org/R/image_matrix.html 2011-01-06 15:35

  min <- min(M)
  max <- max(M)
  yLabels <- rownames(M)
  xLabels <- colnames(M)
  title <-c()
  # check for additional function arguments
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
      min <- Lst$zlim[1]
      max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ){
      yLabels <- c(Lst$yLabels)
    }
    if( !is.null(Lst$xLabels) ){
      xLabels <- c(Lst$xLabels)
    }
    if( !is.null(Lst$title) ){
      title <- Lst$title
    }
  }
  # check for null values
  if( is.null(xLabels) ){
    xLabels <- c(1:ncol(M))
  }
  if( is.null(yLabels) ){
    yLabels <- c(1:nrow(M))
  }

  layout(matrix(data=c(1, 2), nrow=1, ncol=2), widths=c(4, 1), heights=c(1, 1))

  # Red and green range from 0 to 1 while Blue ranges from 1 to 0
  ColorRamp <- rgb( seq(0, 1, length=256),  # Red
    seq(0, 1, length=256),  # Green
    seq(1, 0, length=256))  # Blue
  ColorLevels <- seq(min, max, length=length(ColorRamp))

  # Reverse Y axis
  reverse <- nrow(M) : 1
  yLabels <- yLabels[reverse]
  M <- M[reverse, ]

  # Data Map
  par(mar = c(3, 5, 2.5, 2))
  image(1:length(xLabels), 1:length(yLabels), t(M), col=ColorRamp, xlab="",
    ylab="", axes=FALSE, zlim=c(min, max))
  if( !is.null(title) ){
    title(main=title)
  }
  axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=0.7)
  axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,
    cex.axis=0.7)

  # Color Scale
  par(mar = c(3, 2.5, 2.5, 2))
  image(1, ColorLevels,
    matrix(data=ColorLevels, ncol=length(ColorLevels), nrow=1),
    col=ColorRamp,
    xlab="", ylab="",
    xaxt="n")

  layout(1)
}

Widths <- function(X) {
  # For a sorted sequence of one-dimensional points X, calculates the width of
  # each point's "territory" (i.e., the set of points closer to that point than
  # any other).
  #
  # Args:
  #   X:  A sorted numeric 1-D vector, representing points on a line.
  #
  # Returns:
  #   The width of the "territory" for each point in X.
  N <- length(X)
  return (diff(c(X[1], (X[-N] + X[-1]) * 0.5, X[N])))
}

DistanceMatrix <- function(X, X.out=X) {
  # Computes the distance from every point in X to every point in X.out.
  #
  # Args:
  #   X:  A matrix of input points.
  #   X.out.:  A matrix of output points, whose distance to every point in 'X'
  #      is desired.
  #
  # Returns:
  #   A matrix whose [i, j] component gives the Euclidean distance from
  #   X.out[i, ] to X[j, ].
  X <- as.matrix(X)
  X.out <- as.matrix(X.out)
  if (ncol(X) != ncol(X.out)) {
    stop("Trying to find distance between points of different dimensions!")
  }
  X.dim <- ncol(X)
  N <- NumPoints(X)
  N.out <- NumPoints(X.out)

  # Construct two matrices with (N * N.out) rows, so that X's matrix holds
  # N.out copies of each X-point, and X.out's matrix holds N copies of each
  # X.out-point, and every point in X is matched with every point in X.out
  # exactly once.
  M <- matrix(rep(X, each=N.out), ncol=X.dim)
  M.out <- t(matrix(rep(t(X.out), N), nrow=X.dim))
  dist.sq <- rowSums((M - M.out) ^ 2)
  return (matrix(sqrt(dist.sq), nrow=N.out, ncol=N))
}

Anscombe <- function(Y.p) {
  # Calculate the Anscombe transform of some (presumably Poisson-distributed)
  # data.
  #
  # Args:
  #   Y.p:  A numeric vector of (presumably Poisson-distributed) data.
  #
  # Returns:
  #   The Anscombe transform of Y.p: very close to Gaussian distribution, with
  #   variance 1/4.
  return (sqrt(Y.p + 3.0 / 8.0))
}

AnscombeInverse <- function(Y.g) {
  # Calculate the (statistical!) inverse Anscombe transform of some data.
  #
  # Args:
  #   Y.g:  A numeric vector of (presumably Anscombe-transformed) data.
  #
  # Returns:
  #   The inverse Anscombe transform of Y.
  return (Y.g ^ 2 - (1 / 8))
}

InitializeBoundedQuantity <- function(quantity, bounds,
  default=1, ok.range=c(-Inf, Inf), logspace=FALSE) {
  # Process values for a quantity and its boundaries so they are mutually
  # consistent.
  #
  # Args:
  #   quantity:  The bounded quantity, a single numeric variable.
  #   bounds:  Numeric vector of arbitrary length, whose range() gives the
  #      lower and upper bounds.
  #   default:  The default value for quantity, in case bounds and quantity are
  #      both missing or inappropriate.
  #   ok.range:  The range of values which make sense for this quantity (e.g.,
  #      if it is positive-definite, ok.range should be c(0, Inf)).
  #   logspace:  If TRUE, we set the initial value to the *geometric* mean.
  #
  # Returns:
  #   A list x, where x$quantity is the processed value of quantity, and
  #   x$bounds is the processed value of bounds, such that:
  #   (x$bounds[1] <= x$quantity <= x$bounds[2]).
  #
  # Notes:
  #   1) If both quantity and bounds are given, clamp the quantity to the
  #      bounds.
  #   2) If only the quantity is given, set the bounds to the quantity itself
  #      (effectively making it a constant).
  #   3) If only the bounds are given, set the quantity to their mean.
  #   4) If neither is given, set bounds to range, quantity to default.

  have.quantity <- is.numeric(quantity) && length(quantity) == 1
  have.bounds <- is.numeric(bounds)
  if (have.bounds) {
    bounds <- Clamp(x=range(bounds), bounds=ok.range)
    if (!have.quantity) {
      if (logspace) {
        quantity <- exp(mean(log(bounds)))
      } else {
        quantity <- mean(bounds)
      }
    }
    quantity <- Clamp(x=quantity, bounds=bounds)
  } else {
    if (have.quantity) {
      quantity <- Clamp(x=quantity, bounds=ok.range)
      bounds <- range(quantity)
    } else {
      bounds <- ok.range
      quantity <- default
    }
  }
  return (list(quantity=quantity, bounds=bounds))
}

NumPoints <- function(X) {
  # Count the number of points represented by X, hopefully with a little
  # robustness as to its type (matrix, vector, etc.).
  #
  # Args:
  #   X:  A set of points: ideally, a matrix with the number of columns equal
  #      to the dimensionality, though for 1D points it could also turn out to
  #      be a vector.
  #
  # Returns:
  #   The number of points in X.
  N <- length(X)
  if (is.matrix(X)) {
    N <- nrow(X)
  }
  return (N)
}

SetupFileInfo <- function(name, type, name.default="DEFAULT") {
  # Given a filename and file type, either of which may be missing or
  # malformed, return a proper filename and filetype for a file to write to.
  #
  # Args:
  #   name:  Character; the name of the file to save.
  #   type:  Character; the file type (pdf, png, etc.)
  #
  # Returns:
  #   list;
  #     $name = character; the filename to save to.
  #     $type = character; the type of file to save.
  have.name <- is.character(name) && length(name) == 1
  have.type <- is.character(type) && length(type) == 1
  if (!have.name) {
    name <- name.default
  }
  ext.match <- gregexpr(name, pattern="\\.[^.]+", perl=TRUE)
  z <- length(ext.match[[1]])
  name.has.ext <- (ext.match[[1]][1] > 0)
  if (!have.type) {
    if (name.has.ext) {
      ext.start <- ext.match[[1]][z]
      ext.stop <- ext.start + attr(ext.match[[1]], "match.length")[z]
      type <- substr(x, ext.start + 1, ext.stop)
    } else {
      type <- ""
    }
  }
  if (!name.has.ext && (type != "")) {
    name <- paste(sep='', name, '.', type)
  }
  return (list(name=name, type=type))
}

BubblingRandomMatrix <- function(n.pts, N=10, n.times=100, ...) {
  # Compute a matrix whose rows appear to "bubble" independently, i.e., each
  # row is a periodic timetrace which linearly interpolates between draws from
  # a normal distribution.
  #
  # Args:
  #   n.pts:  The number of independent bubblers.
  #   N:  Half the number of independent draws to take for each bubbler.
  #   n.times:  The final number of interpolated points.
  #
  # Returns:
  #   A numeric matrix with n rows and n.times columns, showing the traces
  #   of the datapoints.

  # First: "seed matrix", which contains enough independent random draws
  max.time <- n.draws <- 2 * N
  d.t <- max.time / n.times
  times <- seq(from=d.t, to=max.time, by=d.t)
  n.seed <- n.pts * n.draws
  seed <- matrix(rnorm(n=n.seed), ncol=n.draws)

  # Next: "Fourier matrix", which gives trigonometric basis functions
  F.sin <- t(sapply(X=1:N, FUN=function(i, N, times) sin(pi * i * times / N),
    N=N, times=times))
  F.cos <- t(sapply(X=1:N, FUN=function(i, N, times) cos(pi * i * times / N),
    N=N, times=times))
  F.tot <- rbind(F.sin, F.cos) / sqrt(N)

  # Put them together and we're done!
  return (seed %*% F.tot)
}

Spaces <- function(num=0) {
  # Constructs a string consisting of the requested number of spaces.  (Useful
  # for recursive tab alignment of printed output.)
  #
  # Args:
  #   num:  The number of spaces to return.
  #
  # Returns:
  #   A string with 'num' spaces.
  return (paste(collapse='', rep(' ', num)))
}

Wrap <- function(x, wrapper) {
  # Wrap each entry in 'x' with the string in 'wrapper'
  #
  # Args:
  #   x:  Character vector whose entries are to be wrapped.
  #   wrapper:  Character string to wrap it in.
  #
  # Returns:
  #   Strings in the format "<wrapper><x><wrapper>".

  # NB: I've started working on bracket matching, in the function
  # InvertPairedSymbols() in testing.R
  return (paste(sep='', wrapper, x, wrapper))
}

PrintParams <- function(lower, params, upper, indent=0, ...) {
  # Print out a list of parameters, together with their upper and lower bounds.
  #
  # Args:
  #   lower: A numeric vector holding the lower bounds.
  #   params: A numeric vector holding the parameter values.
  #   upper: A numeric vector holding the upper bounds.
  # Ruler:
  #     1...5...10...15...20...25...30...35...40...45...50...55...60...65...70
  s <- "  PARAMETER NAME:         MINIMUM:   CURRENT VALUE:         MAXIMUM:\n"
  tab <- Spaces(num=indent)
  cat(tab, s, sep='')
  p.names <- names(params)
  for (name in p.names) {
    cat(sprintf("%s  %15s   %14.8g   %14.8g   %14.8g\n", tab, name,
        lower[name], params[name], upper[name]))
  }

}

#-------------------------------------------------------------------------------
# 2D computational geometry functions

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

HexagonalGrid <- function(n, base.point=c(0, 0), unit=1) {
  # Constructs a hexagonal grid of points with 'n' layers, spaced by 'unit',
  # and centred at 'base.point'.
  #
  # Args:
  #   n:  The number of shells to construct.
  #   base.point:  The centre of the shell.
  #   unit:  The spacing between neighboring points.
  #
  # Returns:
  #   A ((1+3n(n+1)) x 2) matrix, whose (1+3n(n+1)) rows give the points in
  #   this grid, in counterclockwise order and spiraling outwards.
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

GriddedConvexHull <- function(X, spacing) {
  # Computes a hexagonal grid of the given 'spacing' inside the convex hull of
  # the points 'X'.
  #
  # Args:
  #   X:  2-column matrix of points whose convex hull defines the region to
  #      grid.
  #   spacing:  The distance between points in the grid.
  #
  # Returns:
  #   A 2-column matrix of points in a hexagonal grid covering this region.
  grid.centre <- matrix(ncol=2, c(mean(range(X[, 1])), mean(range(X[, 2]))))
  c.hull <- X[SortedConvexHullIndices(X), ]
  n.shells <- ceiling(max(DistanceMatrix(grid.centre, c.hull)) / (
      spacing * cos(pi / 6)))
  hex.grid <- HexagonalGrid(n=n.shells, base.point=grid.centre, unit=spacing)
  clipped.grid <- ClipPointsToHull(X=hex.grid, sorted.hull=c.hull)
  return (clipped.grid)
}

FindGapPoints <- function(X, X.out, gap.width=NA) {
  # Flags "gap points": points in X.out which are further from every point in
  # X, than the largest nearest-neighbor distance in X.
  #
  # Args:
  #   X:  d-column matrix whose rows represent points where we have data.
  #   X.out:  d-column matrix whose rows represent points where we want to
  #      predict.
  #
  # Returns:
  #   A logical matrix, of length NumPoints(X.out), TRUE if the corresponding
  #   point in X.out is a "gap point".

  # Find the longest nearest-neighbor distance.
  if (!(is.numeric(gap.width) && gap.width > 0)) {
    DM <- as.matrix(dist(X))
    diag(DM) <- NA
    gap.width <- max(apply(DM, 2, min, na.rm=TRUE))
    rm(DM)
  }
  # Find the closest X-point for each X.out-point.
  X.out.closest <- apply(DistanceMatrix(X=X, X.out=X.out), 1, min)
  return (X.out.closest > gap.width)
}

#-------------------------------------------------------------------------------
# rgl plotting functions

PlotSurface <- function(X, Y, dY=NA, tri=NA, new.window=FALSE,
  Y.scale=max(dist(X)) / diff(range(Y)), ...) {
  # Plot a triangulated surface using rgl.
  #
  # Args:
  #   X:  2-column matrix holding the X-points for data.
  #   Y:  numeric vector of datapoints
  #   dY:  optional numeric vector of uncertainties in datapoints
  #   tri:  3-column matrix of indices into X, whose rows are triangles.
  #   new.window:  Whether to open a new window
  #   Y.scale:  The factor for scaling Y when plotting.
  #
  # Returns:
  #   Used for its side effect.
  X <- as.matrix(X)
  if (!require("geometry") || !require("rgl")) {
    stop("Surface plotting requires packages 'rgl' and 'geometry'.")
  }
  if (is.na(tri)) {
    tri <- delaunayn(p=X)
  }
  if (new.window || rgl.cur() == 0) {
    rgl.open()
  }
  i <- as.vector(t(tri))
  rgl.triangles(X[i, 1], Y[i], X[i, 2], ...)
  if (!isTRUE(is.na(dY))) {
    rgl.triangles(X[i, 1], (Y[i] + dY[i]), X[i, 2], alpha=0.3, ...)
    rgl.triangles(X[i, 1], (Y[i] - dY[i]), X[i, 2], alpha=0.3, ...)
  }
  # Scale the Y-axis if we haven't already
  Y.scale.old <- par3d()$scale[2]
  if (Y.scale.old == 1 && !is.na(Y.scale)) {
    par3d(scale=c(1, Y.scale, 1))
  }
}

#-------------------------------------------------------------------------------
# grid/ggplot2 layout functions

LayoutNewGridPage <- function(Layout, ...) {
  # Setup a new grid page with the specified Layout.  
  # (Adapted from http://bit.ly/9m4zyD)
  #
  # Args:
  #   Layout:  Result of calling grid.layout function with the desired
  #      arrangement.
  #
  # Returns:
  #   Used for its side-effect.
  grid.newpage()
  pushViewport(viewport(layout = Layout))
}

Subplot <- function(x, y) {
  # Terse wrapper to return a viewport into the specified row (x) and 
  # column (y) of a viewport.
  # (Adapted from http://bit.ly/9m4zyD)
  #
  # Args:
  #   x:  Numeric; the row(s) of the layout.
  #   y:  Numeric; the column(s) of the layout.
  #
  # Returns:
  #   Used for its side-effect.
  viewport(layout.pos.row=x, layout.pos.col=y)
}

################################################################################
# SECTION:                     Class Definitions
################################################################################

################################################################################
# CLASS:                           LazyMatrix
#
# A wrapper for matrices to avoid unnecessary recalculation: it remembers the
# most recent ingredients used to calculate it, and skips the recalc if they're
# identical.
#
# Virtual Fields (R = read-only):
# (R) M:  The matrix most recently computed.
#
# Methods:
#   NeedToRecalculate:  TRUE iff we will need to recompute this matrix because
#      the ingredients are different.
#   StoreMatrix:  Keep a copy of the latest version of the matrix, together
#      with the ingredients used to calculate it.
#
# This class doesn't do the actual computation.  It just tells other code
# whether or not it NEEDS to.

setConstructorS3("LazyMatrix",
  function(...) {
    # Constructs a LazyMatrix object.
    #
    # Returns:
    #   The new LazyMatrix object.

    extend(Object(), "LazyMatrix",
      .last.ingredients = list(),
      .matrix           = matrix()
      )
  })

#-------------------------------------------------------------------------------
# (LazyMatrix) PUBLIC VIRTUAL FIELDS:

# Virtual Field: M (read-only)
# The previously-computed matrix.
#
# Returns:
#   The previously-computed matrix.
setMethodS3("getM", "LazyMatrix", conflict="quiet",
  function(this, ...) {
    return (this$.matrix)
  })

#-------------------------------------------------------------------------------
# (LazyMatrix) PUBLIC METHODS:

setMethodS3("NeedToRecalculate", "LazyMatrix", conflict="quiet",
  function(this, ingredients, ...) {
    # Decides whether or not we would need to recompute the matrix based on the
    # supplied ingredients.
    #
    # Args:
    #   ingredients:  A named list of quantities this matrix depends on.
    #
    # Returns:
    #   TRUE iff we will need to recompute this matrix because the ingredients
    #   are different.

    # First, make sure we have all the old ingredients, and no new ones
    same.names <- (
      all(names(ingredients) %in% names(this$.last.ingredients)) &&
      all(names(this$.last.ingredients) %in% names(ingredients))
      )
    if (!same.names) {
      return (TRUE)
    }

    # Make sure the contents of corresponding ingredients are identical.
    for (name in names(ingredients)) {
      if (!identical(ingredients[[name]], this$.last.ingredients[[name]])) {
        return (TRUE)
      }
    }

    return (FALSE)
  })

setMethodS3("StoreMatrix", "LazyMatrix", conflict="quiet",
  function(this, M, ingredients, ...) {
    # Stores the supplied matrix in memory, along with the ingredients used to
    # calculate it.
    #
    # Args:
    #   M:  The matrix to store.
    #   ingredients:  A list of quantities used to compute 'M'.
    #
    # Returns:
    #   Used for its side-effect (i.e., storing these quantities).
    this$.matrix <- M
    this$.last.ingredients <- ingredients
    return (invisible(this))
  })

################################################################################
# CLASS:                            Dataset
#
# A wrapper class to hold data being analyzed, and facilitate interacting with
# that data.
#
# Virtual Fields (R = read only):
# (R) d:  The number of dimensions for X, the independent variable.
# (R) dataOffset:  The "zero value" for this column, i.e., the mean of the
#    Gaussian Process prediction.
# (R) dpts:  The values of the quantity being predicted.
# (R) xformedDpts:  Values of the quantity being predicted, transformed (if
#    necessary) to have a Gaussian distribution.
# ( ) id:  A character string used to identify this object.
# (R) isPoisson:  logical; TRUE if the current column has Poisson data.
# (R) n:  The number of datapoints.
# (R) noiseVar:  The uncertainty per datapoint in the current quantity, assumed
#    to be Gaussian.
# ( ) quantity:  The name (column title) of the quantity being analyzed.
# (R) X:  The covariates, or independent variables, for this Dataset.
#
# Methods:
#   DeleteRows:  Remove the rows with the specified indices.
#   MSR:  The mean square residuals, in the space of the transformed data.
#   Plot2D:  Plot a 2D (i.e. function-of-two-variables) dataset using rgl.
#   print:  Pretty-print this Dataset.
#   RemoveRange:  Like DeleteRows, but specify an X-range instead of the
#      indices.
#   Same:  Check for equality on multiple fields.
#   SameX:  Check whether this dataset has the same X-values as another dataset.
#   SaveData:  Write a binary snapshot of the current state of the object.
#   Untransform:  Apply, to the given numbers, the inverse of the
#      transformations we applied to our data.

setConstructorS3("Dataset",
  function(id="UNNAMED", data=data.frame(), X.names="X", noise.var=c(), column="",
    poisson.names=c(), tol.factor=1e-4, data.offset=c(), ...) {
    # Constructs a Dataset object with the given ID holding the given data.
    #
    # Args:
    #   id:  A character string which will be used to refer to this object.
    #   data:  A data.frame holding the data to be analyzed
    #   noise.var:  Custom values for the (assumed-Gaussian) uncertainty in
    #      each data column
    #   tol.factor:  The relative amount of noise (compared to the standard
    #      deviation of a data column, which is taken as a proxy for "amount of
    #      signal") which can be tolerated for the sake of numerical stability.
    #   X.names:  A list of names indicating which columns in 'data' should be
    #      considered as the independent variables.
    #
    # Returns:
    #   Dataset object wrapping the given data with the given ID.

    # Ignore factors
    data.trim <- data
    factor.i <- which(sapply(data, class) == 'factor')
    if (length(factor.i) > 0) {
      data.trim <- data[, -factor.i]
    }

    # Only choose the X columns which actually exist, and designate the rest as
    # candidate "data" columns.
    existing.X.names <- X.names[which(X.names %in% colnames(data.trim))]
    other.names <- colnames(data.trim)[
      -which(colnames(data.trim) %in% existing.X.names)]

    # Pick a default column for analysis.
    column <- column[1]
    if (length(other.names) > 0 && !any(column %in% other.names)) {
      column <- other.names[1]
    }

    # Compute the noise variance.
    n.var <- ComputeNoiseVariance(data=data.trim, specs=noise.var, tol=tol.factor,
      poisson.names=poisson.names)

    # Compute data offsets, according to the following logic:
    # 1. Any named entries get set to the given value.
    # 2. The first non-named entry (if any) is a default value.
    # 3. Unassigned entries get set to the default; if none exists, they are
    #    set to the mean of the data.
    if (length(data.offset) > 0 && is.null(names(data.offset))) {
      names(data.offset) <- rep('', length(data.offset))
    }
    unnamed <- which(nchar(names(data.offset)) == 0)
    default.offset <- ifelse(length(unnamed) > 0, data.offset[unnamed[1]], NA)
    offsets <- rep(default.offset, length(other.names))
    names(offsets) <- other.names
    # Add specifically-named values
    for (n in names(data.offset)[-unnamed]) {
      offsets[n] <- data.offset[n]
    }
    # Change NA to mean
    for (i in which(is.na(offsets))) {
      n <- names(offsets)[i]
      is.poisson <- n %in% poisson.names
      dpts.gauss <- dpts <- data.trim[, n]
      if (is.poisson) {
        dpts.gauss <- Anscombe(dpts)
      }
      offsets[n] <- mean(dpts.gauss)
    }

    extend(Object(), "Dataset",
      .data.frame      = data.trim,
      .data.full       = data,
      .data.offset     = offsets,
      .id              = id,
      # The computed noise variance for each column, based on the above specs:
      .noise.var       = n.var,
      # The user's original specifications about the noise variance (useful to
      # have, in case we delete rows and want to reconstruct noiseVar):
      .noise.var.specs = noise.var,
      # Columns which are Poisson distributed.
      .poisson         = colnames(data)[which(colnames(data) %in% poisson.names)],
      # The column we're actually using for predictions:
      .quantity        = column,
      .tol.factor      = tol.factor,
      .X.names         = existing.X.names
      )
  })

#-------------------------------------------------------------------------------
# (Dataset) "HELPER" FUNCTIONS:

ComputeNoiseVariance <- function(data, specs, tol, poisson.names=c()) {
  # Compute the noise variance for each column in 'data', using any custom
  # values in 'specs', or the relative tolerance 'tol' for entries lacking such
  # custom values.
  #
  # Args:
  #   data:  A data.frame holding the columns to specify.
  #   specs:  A named vector of values giving custom noise variances for
  #      correspondingly-named columns in 'data'.
  #   tol:  The default amount of noise we can tolerate in a column, relative
  #      to the sample standard deviation of the data in that column.
  #   poisson.names:  A list of columns which have the Poisson distribution;
  #      *these* columns should have a default value of 0.25, based on the
  #      Anscombe transform.
  #
  # Returns:
  #   A vector with one named entry per column of 'data', giving the Gaussian
  #   variance of the noise in this column.
  #
  # RULES FOR NOISE VARIANCE:
  # 1) Any column value manually specified takes top precedence
  # 2) Otherwise, if a column is Poisson, it gets set to 0.25
  # 3) Un-named values in 'specs':
  #  a) Any unnamed values get doled out as defaults
  #  b) If there IS NO default value, set unspecified to (tol * sd), to
  #     aid numerical stability
  #
  # NOTE that it is computationally advantageous to have nonzero noise in
  # almost every case.  Otherwise, K-matrices tend to be very close to
  # singular.  That's the idea behind specifying a relative tolerance: it makes
  # a relatively painless way to say you don't care about fluctuations smaller
  # than 'tol' (in units of your data's fluctuations, which we guess to be in
  # the same ballpark as the sample standard deviation).  Really though, users
  # should look at their data and manually decide what level fluctuations they
  # could live with for each column.

  # Setup the object
  noise.var <- rep(NA, ncol(data))
  names(noise.var) <- colnames(data)

  # Rule 1 implemented:
  given.names <- names(specs)[which(names(specs) %in% names(noise.var))]
  noise.var[given.names] <- specs[given.names]

  # Rule 2 implemented:
  is.Poisson <- names(noise.var) %in% poisson.names
  noise.var[which(is.na(noise.var) & is.Poisson)] <- 0.25

  # Rule 3 implemented:
  # 'leftover': names of columns which don't have a noise variance yet
  leftover <- names(noise.var)[which(is.na(noise.var))]
  if (length(leftover) > 0) {
    # 'defaults': initialize with the "Rule 3b" values
    defaults <- (tol * sd(data[, leftover])) ^ 2
    # If specs has any named entries, it's safe to look for unnamed ones
    # (i.e., Rule 3a)
    if (!is.null(names(specs)) ) {
      supplied.defaults <- specs[which(is.na(names(specs)) | names(specs)=="")]
      if (length(supplied.defaults) > 1) default <- supplied.defaults
    }
    noise.var[leftover] <- defaults
  }
  return (noise.var)
}

DatasetFromFile <- function(file, sep='\t', id=basename(file), ...) {
  # Read a file's contents into a data.frame; then, use that data.frame to
  # construct a Dataset object.
  #
  # Args:
  #   file:  A character string holding the name of the file to read.
  #   sep:  The character string which separates entries on a line.
  #
  # Returns:
  #   A Dataset object holding the data in the given file
  data.from.file <- read.table(header=TRUE, file=file, sep=sep)
  # Default ID is the basename of the file
  return (Dataset(data=data.from.file, id=id, ...))
}

#-------------------------------------------------------------------------------
# (Dataset) PUBLIC VIRTUAL FIELDS:

# Virtual Field: d (read-only)
# The number of dimensions for X, the independent variable.
#
# Returns:
#   The id of this Dataset object.
setMethodS3("getD", "Dataset", conflict="quiet",
  function(this, ...) {
    return (length(this$.X.names))
  })

# Virtual Field: dataOffset (read-only)
# The mean of the Gaussian Process for this data.
#
# Returns:
#   The datapoint values for the currently-selected quantity.
setMethodS3("getDataOffset", "Dataset", conflict="quiet",
  function(this, ...) {
    return (this$.data.offset[this$getQuantity()])
  })

# Virtual Field: dpts (read-only)
# The datapoint values for the currently-selected quantity.
#
# Returns:
#   The datapoint values for the currently-selected quantity.
setMethodS3("getDpts", "Dataset", conflict="quiet",
  function(this, ...) {
    return (this$.data.frame[, this$getQuantity()])
  })

# Virtual Field: xformedDpts (read-only)
# The transformed datapoint values for the currently-selected quantity, where
#   the transformation gives them a Gaussian distribution.
#
# Returns:
#   A Gaussian-noised version of the current quantity datapoints.
setMethodS3("getXformedDpts", "Dataset", conflict="quiet",
  function(this, ...) {
    gaussian.dpts <- dpts <- this$getDpts()
    if (this$getIsPoisson()) {
      # Anscombe transform: turns Poisson data into Gaussian data with constant
      # variance of 1/4.
      gaussian.dpts <- Anscombe(dpts)
    }
    return (gaussian.dpts - this$getDataOffset())
  })

# Virtual Field: id
# A character string identifying this Dataset object.
#
# Args:
#   id: the string to change the id to.
#
# Returns:
#   The id of this Dataset object.
setMethodS3("getId", "Dataset", conflict="quiet",
  function(this, ...) {
    return (this$.id)
  })
setMethodS3("setId", "Dataset", conflict="quiet",
  function(this, id, ...) {
    this$.id <- id
    return (this)
  })

# Virtual Field: isPoisson (read-only)
# A character string identifying this Dataset object.
#
# Returns:
#   logical; TRUE if the current column has Poisson data.
setMethodS3("getIsPoisson", "Dataset", conflict="quiet",
  function(this, ...) {
    return (any(this$.quantity %in% this$.poisson))
  })

# Virtual Field: n (Read-only)
# The number of datapoints.
#
# Returns:
#   The number of datapoints.
setMethodS3("getN", "Dataset", conflict="quiet",
  function(this, ...) {
    return (nrow(this$.data.frame))
  })

# Virtual Field: noiseVar (Read-only)
# The noise variance for the current column.
#
# Returns:
#   The noise variance for the current column.
setMethodS3("getNoiseVar", "Dataset", conflict="quiet",
  function(this, ...) {
    return (this$.noise.var[this$getQuantity()])
  })

# Virtual Field: quantity
# The name (i.e., column title) of the quantity currently being analyzed.
#
# Args:
#   col.name:  The name of the new quantity to be analyzed.
#
# Returns:
#   The name of the quantity currently being analyzed.
setMethodS3("getQuantity", "Dataset", conflict="quiet",
  function(this, ...) {
    return (this$.quantity)
  })
setMethodS3("setQuantity", "Dataset", conflict="quiet",
  function(this, col.name, ...) {
    if (col.name %in% colnames(this$.data.frame)) {
      this$.quantity <- col.name
    } else {
      warning(sprintf("You have selected a column ('%s') which does not exist! Column NOT changed.\n", col.name))
    }
    return (invisible(this))
  })

# Virtual Field: X (read-only)
# A (n x d) matrix giving the X-values (i.e., the independent variable) for
# this Dataset (n being the number of datapoints, d being the dimensionality).
#
# Returns:
#   The X-values for this Dataset.
setMethodS3("getX", "Dataset", conflict="quiet",
  function(this, ...) {
    return (as.matrix(this$.data.frame[, this$.X.names]))
  })

#-------------------------------------------------------------------------------
# (Dataset) PUBLIC METHODS:

setMethodS3("DeleteRows", "Dataset", conflict="quiet",
  function(this, indices, ...) {
    # Remove rows (i.e., datapoints) from this dataset.
    #
    # Args:
    #   indices:  The indices of the datapoints to be removed.
    #
    # Returns:
    #   'this', but used for its side-effect of deleting rows.
    this$.data.frame <- this$.data.frame[-indices, ]
    this$.noise.var <- ComputeNoiseVariance(data=this$.data.frame,
      specs=this$.noise.var.specs, tol=this$.tol.factor,
      poisson.names=this$.poisson)
    return (invisible(this))
  })

setMethodS3("MSR", "Dataset", conflict="quiet",
  function(this, test.data, reference.data=this$xformedDpts, ...) {
    # The mean square residuals of test.data, versus reference.data (default is
    # to compare to our datapoints).
    #
    # Args:
    #   test.data:  A sequence of datapoints, which live in the same space as
    #      xformedDpts (i.e., same offset, and Anscombe-transformed if
    #      Poisson).
    #   reference.data:  The comparison function; defaults to our noisy
    #      datapoints, but could also be different (e.g., if we have access to
    #      the true function).
    #
    # Returns:
    #   The mean square difference between test.data and reference.data.
    if (length(test.data) != length(reference.data)) {
      stop("Trying to find MSR for sequences of different lengths")
    }
    return (mean((test.data - reference.data) ^ 2))
  })

setMethodS3("Plot2D", "Dataset", conflict="quiet",
  function(this, max.points=1000, dist.factor=0.2, ...) {
    # Plot a 2D (i.e. function-of-two-variables) dataset using rgl.
    #
    # max.points:  The maximum number of points to show (these are randomly
    #    sampled from the available points).
    #
    # Returns:
    #   Used for its side-effect.
    d <- clone(this)
    if (max.points < this$n) {
      i <- sample(x=1:this$n, size=max.points, replace=FALSE)
      d$DeleteRows(-i)
    }
    unit <- min(dist(d$X))
    if (require("rgl") == FALSE) {
      stop("The Plot2D() method requires the rgl library to be installed.")
    }
    rgl.clear()
    rgl.spheres(x=d$X[, 1], z=d$X[, 2], y=d$dpts, radius=dist.factor * unit)
  })

print.Dataset <- function(this, ...) {
  # Print a summary of 'this' Dataset object.
  #
  # Returns:
  #   Used for its side-effect.
  cat(sprintf("Dataset '%s': %d datapoints in %d dimensions.\n",
      this$id, this$n, this$d))
  # Ruler:
  #     1...5...10...15...20...25...30...35...40...45...50...55...60...65...70
  s <- "            NAME:            RANGE:"
  s.quantities <- paste(sep='', s,       "      NOISE VARIANCE:\n")
  s.covariates <- paste(sep='', s,       "\n")
  format.str <- "%17s%4s%14.8g%21s\n"
  raw.data <- this$.data.frame
  quantity.names <- colnames(raw.data)
  decoration <- "------------------"
  cat(sprintf("%16s%s%s%s\n", '', decoration, 'Covariates:', decoration))
  cat(s.covariates)
  X.indices <- which(quantity.names %in% this$.X.names)
  for (X.i in X.indices) {
    cat(sprintf(format.str, quantity.names[X.i], '',
        diff(range(raw.data[, X.i])), ''))
  }
  cat(sprintf("\n%16s%s%s%s\n%18s%s\n", '', decoration, 'Quantities:',
      decoration, '', "('*' = currently selected; 'P' = Poisson)"))
  cat(s.quantities)
  for (X.i in (1:ncol(raw.data))[-X.indices]) {
    sel <- ifelse(quantity.names[X.i] == this$quantity, '*', ' ')
    pois <- ifelse(quantity.names[X.i] %in% this$.poisson, 'P', ' ')
    flag <- sprintf('(%s%s)', sel, pois)
    cat(sprintf(format.str, quantity.names[X.i], flag,
        diff(range(raw.data[, X.i])),
        sprintf('%17.8g', this$.noise.var[X.i])))
  }
  return (invisible(this))
}

setMethodS3("RemoveRange", "Dataset", conflict="quiet",
  function(this, X.min=min(this$X), X.max=max(this$X), ...) {
    # Remove datapoints whose X-values fall within the specified range.
    #
    # Args:
    #   X.min:  The left boundary of the X-range to remove.
    #   X.max:  The right boundary of the X-range to remove.
    #
    # Returns:
    #   'this', but used for its side-effect of removing datapoints.
    if (this$d > 1) {
      stop("RemoveRange() only works for 1-D data.")
    }
    if (X.min < X.max) {
      bad.indices <- which(this$X <= X.max & this$X >= X.min)
      this$.data.frame <- this$.data.frame[-bad.indices, ]
      this$.noise.var <- ComputeNoiseVariance(data=this$.data.frame,
        specs=this$.noise.var.specs, tol=this$.tol.factor,
        poisson.names=this$.poisson)
    }
    return (invisible(this))
  })

setMethodS3("Same", "Dataset", conflict="quiet",
  function(this, d, compare, ...) {
    # Check whether Dataset 'this' has the same values as Dataset 'd'.
    #
    # Args:
    #   d:  The Dataset for comparison.
    #   compare:  The attributes (e.g., X, noiseVar, etc.) to check for
    #      equality.
    #
    # Returns:
    #   TRUE if all the requested attributes are the same

    # First: guard against NULLs and the like!
    if (!("Dataset" %in% class(d) && "Dataset" %in% class(this))) {
      return (FALSE)
    }

    # A convenient helper function:
    field.passes.test <- function(field) {
      # Check whether a *single* 'field' is identical for Datasets 'this' and
      # 'd'.
      #
      # Args:
      #   field:  The field (e.g., X, noiseVar, etc.) to check.
      #
      # Returns:
      #   FALSE if this field is in the list that gets checked, *and* its value
      #   differs between 'this' and 'd'; otherwise, TRUE.
      if (!(field %in% compare)) {
        return (TRUE)
      }
      # Locate accessor function; check that both objects return equal values.
      substr(field, 1, 1) <- toupper(substr(field, 1, 1))
      funct <- get(paste(sep='', 'get', field))
      return (identical(funct(this), funct(d)))
    }

    # Check all the requested fields
    ok <- TRUE
    fields.to.check <- c(
      'X', 'noiseVar', 'dpts', 'xformedDpts', 'id', 'quantity')
    for (field in fields.to.check) {
      ok <- ok && field.passes.test(field)
    }

    # Warn user if any fields were requested that we don't check!
    not.checked <- compare[-which(compare %in% fields.to.check)]
    if (length(not.checked) > 0) {
      offenders <- paste(sep='', "'", not.checked, "'", collapse=', ')
      warning(sprintf(
          "The following fields are not checked by Dataset$Same:\n%s\n",
          offenders))
    }
    return (ok)
  })

setMethodS3("SameX", "Dataset", conflict="quiet",
  function(this, d, ...) {
    # Check whether Dataset 'this' has the same X-values as Dataset 'd'.
    #
    # Args:
    #   d:  The Dataset for comparison.
    #
    # Returns:
    #   TRUE if all the X-values are the same (and in the same order)
    return (this$Same(d, "X"))
  })

setMethodS3("Untransform", "Dataset", conflict="quiet",
  function(this, values, ...) {
    # Apply, to the given numbers, the inverse of the transformations we
    # applied to our data.
    #
    # Args:
    #   values:  A numeric vector whose values "live" in the *transformed*
    #      scale.
    #
    # Returns:
    #   Numbers corresponding to 'values', but in a non-transformed scale
    #   (i.e., instead of being comparable to xformedDpts, the returned
    #   quantities are comparable to dpts).
    values <- values + this$dataOffset
    if (this$getIsPoisson()) {
      values <- AnscombeInverse(values)
    }
    return (values)
  })

################################################################################
# CLASS:                           Covariance
#
# The superclass for all the various kinds of covariance objects (Squared
# Exponential, periodic, Matern, etc.).
#
# Helper Functions:
#   DecodeForTraining:  Takes a named vector of parameters (or lower or upper
#      bounds); transforms any tagged as "logspace" back to linear (changing names appropriately).
#   LogspaceTag:  A tag we append to variable names to indicate they have been
#      transformed into logspace.
#
# Virtual Fields (V = pure virtual, R = read only):
# (  ) id:  A character string used to identify this object.
# ( V) logspaceNames:  A vector of names indicating which parameters should be
#        optimized in logspace.
# ( V) params:  A named vector of parameters governing this Covariance (names
#        prefaced by 'id').
# ( V) paramsPlain:  A named vector of parameters governing this Covariance
#        names NOT prefaced by 'id').
# ( V) lower:  A named vector of lower bounds on parameters governing this
#        Covariance (names prefaced by 'id').
# ( V) lowerPlain:  A named vector of lower bounds on parameters governing this
#        Covariance (names NOT prefaced by 'id').
# ( V) upper:  A named vector of upper bounds on parameters governing this
#        Covariance (names prefaced by 'id').
# ( V) upperPlain:  A named vector of parameters governing this Covariance
#        (names NOT prefaced by 'id').
# (RV) paramNames:  The names of the parameters, decorated by 'id'.
# (RV) paramNamesPlain:  The names of the parameters, NOT decorated by 'id'.
#
# Methods:
#   EncodeForTraining:  Takes a named vector of parameters (or lower or upper
#      bounds); transforms "scale"-type parameters into logspace and change
#      their names.
#   FixConstParam:  Set lower bounds, upper bounds, and current value (for one
#     particular parameter) to the value specified by the user.
#   KInIn, KInOut, KOutIn, KOutOut:  Covariance matrices, based on the current
#      parameter values and a supplied Dataset.  (Function name indicates which
#      set of X-values we are calculating the covariance for: i.e., KInOut is a
#      matrix whose rows correspond to X.out (the points where we want
#      predictions), and whose columns correspond to X (the points where we
#      have data).)
#   KInInDeriv:  Element-wise derivative of 'KInIn', with respect to a given
#      named parameter.
#   PrependId:  Prepend the Covariance id to a bunch of strings.
#   print:  Prettied-up summary of this Covariance object.
#   Variance:  The diagonal of the covariance matrix.
#
# Regarding parameter names: the SUBCLASS has the responsibility to provide
# "plain-named" versions of all the virtual fields
# (paramNamesPlain, paramsPlain, lowerPlain, etc.).
# The SUPERCLASS will automatically handle the "decorated" versions
# (paramNames, params, lower, etc.).
#
# All Covariance subclasses should remember the param-values and Dataset they
# last used to compute their K-matrix.  If they get asked to compute it again,
# they will simply return the previously-computed result if these values have
# not changed.

setConstructorS3("Covariance",
  function(id="", ...) {
    # Constructs a Covariance object with the given ID.
    #
    # Args:
    #   id: A character string which will be used to refer to this object.
    #
    # Returns:
    #   Covariance object with the given ID.

    extend(Object(), "Covariance",
      .id = id
      )
  })

#-------------------------------------------------------------------------------
# (Covariance) PUBLIC VIRTUAL FIELDS:

# Virtual Field: id
# A character string identifying this Covariance object.
#
# Args:
#   id: the string to change the id to.
#
# Returns:
#   The id of this Covariance object.
setMethodS3("getId", "Covariance", conflict="quiet",
  function(this, ...) {
    return (this$.id)
  })
setMethodS3("setId", "Covariance", conflict="quiet",
  function(this, id, ...) {
    this$.id <- id
    return (this)
  })

# Virtual Field: paramNames (ReadOnly)
# The names of the parameters, decorated by the id (if it exists) of this
# Covariance object.
#
# Returns:
#   Parameter names in the form "id.basename"
setMethodS3("getParamNames", "Covariance", conflict="quiet",
  function(this, ...) {
    if (!is.character(getId(this)) || nchar(getId(this)) < 1) {
      return (getParamNamesPlain(this))
    }
    return (this$PrependId(getParamNamesPlain(this)))
  })

# Virtual Field: params
# A named vector of parameters governing this Covariance.
#
# Args:
#   p: A (named) vector of new parameter values (we ONLY use ones which are
#      named, and whose names match up with names of parameters.)
#   for.training:  If TRUE, we ignore "constant" parameters (i.e., where
#      lower=upper) and return the *log* of any "scale" parameters (such as ell
#      or sigma.f for the SE model).
#
# Returns:
#   A vector with values for this Covariance's parameters.
setMethodS3("getParams", "Covariance", conflict="quiet",
  function(this, for.training=FALSE, ...) {
    p <- getParamsPlain(this)
    names(p) <- getParamNames(this)
    if (for.training) {
      p <- this$EncodeForTraining(p)
    }
    return (p)
  })
setMethodS3("setParams", "Covariance", conflict="quiet",
  function(this, p, for.training=FALSE, ...) {
    if (for.training) {
      p <- this$DecodeForTraining(p)
    }
    p.plain <- this$UndecorateNames(p=p)
    this$setParamsPlain(p=p.plain)
    return (invisible(this))
  })

# Virtual Field: paramsPlain
# Gives a vector of parameter values, whose names are NOT decorated by the id
# of this Covariance object.
#
# Args:
#   p: A (named) vector of new parameter values (we ONLY use ones which are
#      named, and whose names match up with names of parameters.)
#
# Returns:
#   The parameters for this covariance function, but with names undecorated by
#   its id.
#
# NOTE: accessor method is responsibility of subclasses.
setMethodS3("setParamsPlain", "Covariance", conflict="quiet",
  function(this, p, ...) {
    # First, we need a vector where every parameter name has a value: either
    # the new value from 'p', or if none was given, the current value.
    p.good.names <- this$getParamsPlain()
    to.change <- names(p)[which(names(p) %in% names(p.good.names))]
    p.good.names[to.change] <- p[to.change]
    # paramsPlainImplementation requires a vector with a value for every
    # parameter; we took care of this above.
    this$paramsPlainImplementation(p=p.good.names)
    p.clamped <- this$ClampedParamVals(p=p.good.names)
    p.names <- names(p.clamped)
    clamped <- which(p.good.names[p.names] != p.clamped[p.names])
    this$paramsPlainImplementation(p=p.clamped)
    if (any(clamped)) {
      culprits <- paste(sep='', collapse=' ', '"', p.names[clamped], '"')
      warning(paste("These parameters had to be clamped:\n", culprits, "\n"))
    }
    return (invisible(this))
  })

# Virtual Field: lower
# Lower bounds for the parameter values.
#
# Args:
#   L: A (named) vector of new lower bounds (we ONLY use ones which are named,
#      and whose names match up with names of parameters.)
#   for.training:  If TRUE, we ignore "constant" parameters (i.e., where
#      lower=upper) and return the *log* of any "scale" parameters (such as ell
#      or sigma.f for the SE model).
#
# Returns:
#   The lower bounds for the parameters for this covariance function.
setMethodS3("getLower", "Covariance", conflict="quiet",
  function(this, for.training=FALSE, ...) {
    L <- this$getLowerPlain()
    names(L) <- this$getParamNames()
    if (for.training) {
      L <- this$EncodeForTraining(L)
    }
    return (L)
})
setMethodS3("setLower", "Covariance", conflict="quiet",
  function(this, L, for.training=FALSE, ...) {
    if (for.training) {
      L <- DecodeForTraining(L)
    }
    L.plain <- this$UndecorateNames(L)
    this$setLowerPlain(L=L.plain)
    return (invisible(this))
})

# Virtual Field: upper
# Upper bounds for the parameter values.
#
# Args:
#   U: A (named) vector of new upper bounds (we ONLY use ones which are named,
#      and whose names match up with names of parameters.)
#   for.training:  If TRUE, we ignore "constant" parameters (i.e., where
#      lower=upper) and return the *log* of any "scale" parameters (such as ell
#      or sigma.f for the SE model).
#
# Returns:
#   The upper bounds for the parameters for this covariance function.
setMethodS3("getUpper", "Covariance", conflict="quiet",
  function(this, for.training=FALSE, ...) {
    U <- this$getUpperPlain()
    names(U) <- this$getParamNames()
    if (for.training) {
      U <- this$EncodeForTraining(U)
    }
    return (U)
  })
setMethodS3("setUpper", "Covariance", conflict="quiet",
  function(this, U, for.training=FALSE, ...) {
    if (for.training) {
      U <- this$DecodeForTraining(U)
    }
    U.plain <- this$UndecorateNames(U)
    this$setUpperPlain(U=U.plain)
    return (invisible(this))
  })

#-------------------------------------------------------------------------------
# (Covariance) PUBLIC METHODS:

clone.Covariance <- function(this, ...) {
  # Deep clone this Covariance object (clones LazyMatrix objects too).
  #
  # Returns:
  #   A deep clone of the Covariance object.
  Cov <- clone.Object(this)
  return (Cov)
}

setMethodS3("EncodeForTraining", "Covariance", conflict="quiet",
  function(this, values, ...) {
    # Takes a named (i.e., with the 'id' tag already prepended) vector of
    #    parameters (or lower or upper bounds); transforms "scale"-type
    #    parameters into logspace and change their names.
    #
    # Args:
    #   values:  A named vector of parameter values (with ID already
    #      prepended).
    #
    # Returns:
    #   'values', but with scale-type parameters changed to their log (and
    #   appropriately renamed)
    logspace.names <- this$PrependId(this$logspaceNames)
    i.log <- which(names(values) %in% logspace.names)
    values[i.log] <- log(values[i.log])
    names(values)[i.log] <- paste(sep="", names(values)[i.log], LogspaceTag())
    return (values)
  })

setMethodS3("FixConstParam", "Covariance", conflict="quiet",
  function(this, p.name, p.value, ...) {
    # Sets the parameter with plain-name 'p.name' to the value 'p.value', along
    # with the lower and upper bounds (making it effectively constant).
    #
    # Args:
    #   p.name:  The (undecorated) name of the parameter to change.
    #   p.value:  The new constant value for the named parameter.
    #
    # Returns:
    #   'this', but used for its side-effect.
    p.old <- this$paramsPlain[p.name]
    p <- p.value
    names(p) <- p.name
    # Make the parameter into a constant:
    this$upperPlain <- p.old
    this$lowerPlain <- p.old
    # The order in which we move the values is chosen to preserve the ordering,
    # (upper>param>lower), and thereby avoid generating warnings:
    if (p.value > p.old) {  # Going up
      this$upperPlain <- p
      this$paramsPlain <- p
      this$lowerPlain <- p
    }
    if (p.value < p.old) {  # Going down
      this$lowerPlain <- p
      this$paramsPlain <- p
      this$upperPlain <- p
    }
    return (invisible(this))
  })

setMethodS3("KInIn", "Covariance", conflict="quiet",
  function(this, d, ...) {
    # Gives the covariance matrix, based on the current parameter values, for
    # the supplied dataset.
    #
    # Args:
    #   d:  A Dataset object encapsulating the data to train on.
    #
    # Returns:
    #   The K-matrix for this dataset.
    return (this$K.specific(X=d$X))
  })

setMethodS3("KInOut", "Covariance", conflict="quiet",
  function(this, d, X.out=d$X, ...) {
    # Gives the covariance matrix, based on the current parameter values, for
    # the supplied dataset and prediction points.
    #
    # Args:
    #   d:  A Dataset object encapsulating the data to train on.
    #   X.out:  A matrix (with d$d columns) of X-locations where we want
    #      predictions.
    #
    # Returns:
    #   The K-matrix between the X-values for this dataset, and the X-values
    #   where we want to predict.

    # Recompute iff necessary:
    return (this$K.specific(X=d$X, X.out=X.out))
  })

setMethodS3("KOutIn", "Covariance", conflict="quiet",
  function(this, d, X.out=d$X, ...) {
    # Gives the covariance matrix, based on the current parameter values, for
    # the supplied prediction points and dataset.
    #
    # Args:
    #   d:  A Dataset object encapsulating the data to train on.
    #   X.out:  A matrix (with d$d columns) of X-locations where we want
    #      predictions.
    #
    # Returns:
    #   The K-matrix between the X-values where we want to predict, and the
    #   X-values for this dataset.
    return (t(this$KInOut(d=d, X.out=X.out)))
  })

setMethodS3("KOutOut", "Covariance", conflict="quiet",
  function(this, X.out, ...) {
    # Gives the covariance matrix, based on the current parameter values, for
    # the supplied dataset.
    #
    # Args:
    #   X.out:  A matrix (with d$d columns) of X-locations where we want
    #      predictions.
    #
    # Returns:
    #   The K-matrix for this dataset.

    # Recompute iff necessary:
    return (this$K.specific(X=X.out))
  })

setMethodS3("KInInDeriv", "Covariance", conflict="quiet",
  function(this, d, param, ...) {
    # Calculate element-wise derivative of the input-input covariance matrix,
    # with respect to Dataset 'd', and the parameter whose (decorated) name is
    # 'param'.
    #
    # Args:
    #   d:  The Dataset whose X-values we are training on.
    #   param:  The name of the parameter with respect to which we're
    #      differentiating .
    #
    # Returns:
    #   The element-by-element derivative of the total K-matrix for Dataset
    #   'd', with respect to the parameter named 'param'.
    #
    # Notes:
    #   As this is the superclass, all we do here is to check whether this is
    #   our problem (i.e., whether this object actually has a parameter named
    #   'param').  If not, just return a zero matrix.  Otherwise, we let the
    #   subclasses handle the actual computation.
    if (param %in% this$paramNames) {
      names(param) <- param
      param.plain <- names(this$UndecorateNames(p=param))
      K.deriv <- this$KDerivImplementation(d=d, param=param.plain)
    } else {
      K.deriv <- matrix(0, nrow=d$n, ncol=d$n)
    }
    return (K.deriv)
  })

setMethodS3("PrependId", "Covariance", conflict="quiet",
  function(this, strings, ...) {
    return (paste(sep='.', this$getId(), strings))
  })

print.Covariance <- function(this, indent=0, ...) {
  # Prints the Covariance object in a pleasing way.
  #
  # Returns:
  #   Used for its side-effect.
  tab <- paste(collapse='', rep(' ', indent))
  cat(sprintf("%sCovariance, id=%-20s (%s)\n", tab, Wrap(this$id, "'"),
      class(this)[1]))
  PrintParams(lower=this$lowerPlain, params=this$paramsPlain,
    upper=this$upperPlain, indent=indent, ...)
  return (invisible(this))
}

#-------------------------------------------------------------------------------
# (Covariance) PUBLIC "HELPER" FUNCTIONS:

LogspaceTag <- function() {
  # A tag we append to variable names to indicate they have been transformed
  # into logspace.
  #
  # Returns:
  #   A string: this tag.
  return (".__LOG__")
}

DecodeForTraining <- function(values, ...) {
  # Takes a named vector of parameters (or lower or upper bounds); transforms
  #    any tagged as "logspace" back to linear (changing names
  #    appropriately).
  #
  # Args:
  #   values:  A named vector of parameter values
  #
  # Returns:
  #   'values', but with log-transformed parameters changed back to linear
  #   scale (and appropriately renamed)
  log.pattern <- paste(sep='', LogspaceTag(), '$')
  i.log <- grep(pattern=log.pattern, x=names(values))
  values[i.log] <- exp(values[i.log])
  names(values)[i.log] <- sub(pattern=log.pattern, replacement='',
    x=names(values)[i.log])
  return (values)
}

#-------------------------------------------------------------------------------
# (Covariance) PRIVATE METHODS:

setMethodS3("ClampedParamVals", "Covariance", private=TRUE, conflict="quiet",
  # Tells what the current parameter values *would* be, if clamped to lie
  # within the lower and upper bounds.  (Naturally, this is intended for
  # *internal* use, by every function which allows the parameter values to
  # change: they are expected to use *this* function so the parameter values
  # *always* stay within bounds.)
  #
  # Args:
  #   warn: Whether or not to print a warning in case params are outside
  #     bounds.
  #
  # Returns:
  #   A named vector of the new parameter values (with plain names).
  function(this, warn=TRUE, p=this$getParamsPlain(), ...) {
    lower <- this$getLowerPlain()
    upper <- this$getUpperPlain()
    p.new <- ClampNamed(x=p, lower=lower, upper=upper)
    p.p <- p
    if (warn && !isTRUE(all.equal(p.new[names(p)], p[names(p)]))) {
      culprits <- paste(sep='', collapse=' ', '"',
        names(p)[which(p != p.new)], '"')
      warning(paste("These parameters were outside the bounds:\n",
          culprits, "\n"))
    }
    return (p.new[names(p)])
  })

setMethodS3("ClampParams", "Covariance", private=TRUE, conflict="quiet",
  function(this, warn=TRUE, ...) {
    # Clamps the parameter values to fall between the upper and lower
    # boundaries; also, invalidates K if they change and optionally prints a
    # warning.
    #
    # Args:
    #   warn: If TRUE, we print a warning if anything actually changes.
    #
    # Returns:
    #   Object whose parameters got clamped, but used for its side-effects.

    # The following *should* be OK, as long as neither setParams nor
    # setParamsPlain calls ClampParams...
    this$setParamsPlain(p=this$ClampedParamVals(warn=warn))
    return (invisible(this))
  })

setMethodS3("UndecorateNames", "Covariance", private=TRUE, conflict="quiet",
  function(this, p, ...) {
    # Returns entries in 'p' which match 'this$id', but with the id stripped
    # out of the names vector.
    #
    # Args:
    #   p: A (named) vector of values corresponding to the parameters of this
    #      Covariance (they could be the parameter values themselves, or
    #      perhaps upper or lower bounds).
    #
    # Returns:
    #   'p' stripped down so it only has the relevant entries with plain names.

    # Which entries in 'p' correspond to our parameters?
    relevant.indices <- which(names(p) %in% this$getParamNames())
    p.relevant <- p[relevant.indices]
    # Strip out our ID string
    id.prefix <- paste(sep='', this$getId(), '.')
    names(p.relevant) <- gsub(id.prefix, '', names(p.relevant))
    return (p.relevant)
  })

setMethodS3("PushUpperBounds", "Covariance", private=TRUE, conflict="quiet",
  function(this, U.min, warn.on.cross=TRUE, ...) {
    # Adjusts upper boundaries to "make way" for the new lower bounds in
    # 'U.min', giving warnings if any adjustments are actually necessary.
    #
    # Args:
    #   U.min:  New lower bounds for corresponding parameters of this Covariance,
    #      where names(U.min) are *already* undecorated.
    #   warn.on.cross:  If TRUE, we will warn if the bounds are actually moved
    #
    # Returns:
    #   Plain-named version of 'U.min' with only the relevant entries (though
    #   this method is mainly used for its side-effect).

    # Which ones will be affected?  (If none, just skip to the finish.)
    to.change <- names(U.min)[which(names(U.min) %in% this$paramNamesPlain)]
    if (length(to.change) < 1) {
      return (c())
    }
    # Remember old values so we can check later whether they've changed
    upper.old <- this$upperPlain
    # If the bounds cross, set them both to the new value
    bounds.crossed <- U.min[to.change] > upper.old[to.change]
    if (any(bounds.crossed)) {
      this$setUpperPlain(U=pmax(upper.old[to.change], U.min[to.change]))
      if (warn.on.cross) {
        culprits <- paste(sep='', collapse=' ', '"',
          to.change[which(bounds.crossed)], '"')
        warning(paste("Upper bounds forcibly moved for these parameters:\n",
            culprits, "\n"))
      }
    }
    return (U.min[to.change])
  })

setMethodS3("PushLowerBounds", "Covariance", private=TRUE, conflict="quiet",
  function(this, L.max, warn.on.cross=TRUE, ...) {
    # Adjusts lower boundaries to "make way" for the new upper bounds in 'L.max',
    # giving warnings if any adjustments are actually necessary.
    #
    # Args:
    #   L.max: New upper bounds for corresponding parameters of this Covariance,
    #      where names(L.max) are *already* undecorated.
    #   warn.on.cross:  If TRUE, we will warn if the bounds are actually moved
    #
    # Returns:
    #   Plain-named version of 'L.max' with only the relevant entries (though
    #   this method is mainly used for its side-effect).

    # Which ones will be affected?  (If none, just skip to the finish.)
    to.change <- names(L.max)[which(names(L.max) %in% this$paramNamesPlain)]
    if (length(to.change) < 1) {
      return (c())
    }
    # Remember old values so we can check later whether they've changed
    lower.old <- this$lowerPlain
    # If the bounds cross, set them both to the new value
    bounds.crossed <- L.max[to.change] < lower.old[to.change]
    if (any(bounds.crossed)) {
      this$setLowerPlain(L=pmin(lower.old[to.change], L.max[to.change]))
      if (warn.on.cross) {
        culprits <- paste(sep='', collapse=' ', '"',
          to.change[which(bounds.crossed)], '"')
        warning(paste("Lower bounds forcibly moved for these parameters:\n",
            culprits, "\n"))
      }
    }
    return (L.max[to.change])
  })

################################################################################
# SUBCLASS:                     CovarianceNoise
#
# Independent pixels, governed only by their (common) variance.
#
# Methods:
#   KDerivImplementation:  The element-wise derivative of KInIn.

setConstructorS3("CovarianceNoise",
  function(id="noise", sigma=NA, sigma.bounds=NA, ...) {
    # Constructs a CovarianceNoise object with the given values and boundaries.
    #
    # Args:
    #   sigma:  The magnitude of the noise.
    #   sigma.bounds:  Range of acceptable values of sigma
    #
    # Returns:
    #   CovarianceNoise object with the given parameter values.

    # Ideas here are the same as for the CovarianceSE constructor, but simpler
    # because there is one fewer parameter.
    sigma.good <- InitializeBoundedQuantity(ok.range=c(0, Inf),
      quantity=sigma, bounds=sigma.bounds, logspace=TRUE)

    # Construct the CovarianceNoise object:
    extend(Covariance(..., id=id), "CovarianceNoise",
      .sigma        = sigma.good$quantity,
      .sigma.bounds = sigma.good$bounds)
  })

# Virtual Field: logspaceNames (read-only)
# A vector of names indicating which parameters should be optimized in
# logspace.
#
# Returns:
#   Names of parameters to be optimized in logspace.
setMethodS3("getLogspaceNames", "CovarianceNoise", conflict="quiet",
  function(this, ...) {
    return (c("sigma"))
  })

# Virtual Field: paramNamesPlain (read-only)
# Gives the "basenames" (i.e. names undecorated by the id string) of the
# parameters.
#
# Returns:
#   The basenames of the parameters.
setMethodS3("getParamNamesPlain", "CovarianceNoise", conflict="quiet",
  function(this, ...) {
    return (c("sigma"))
  })

# Virtual Field: paramsPlain
# Gives a vector of parameter values, whose names are NOT decorated by the id
# of this Covariance object.
#
# SEE ALSO: paramsPlain for superclass Covariance
#
# Args:
#   p: A (named) vector of new parameter values (we ONLY use ones which are
#      named, and whose names match up with names of parameters).  Note that
#      'p' must have a value for every parameter in this Covariance; ensuring
#      this is the responsibility of the calling function.
#
# Returns:
#   The parameters for this covariance function, but with names undecorated by
#   its id.
setMethodS3("getParamsPlain", "CovarianceNoise", conflict="quiet",
  function(this, ...) {
    p <- c(this$.sigma)
    names(p) <- this$getParamNamesPlain()
    return (p)
  })
setMethodS3("paramsPlainImplementation", "CovarianceNoise", conflict="quiet",
  private=TRUE,
  function(this, p, ...) {
    # Sets any values in 'p' which match our parameter names, without worrying
    # about lower and upper bounds (the function which calls this has the job
    # of worrying about these!).
    this$.sigma <- p["sigma"]
    return (invisible(this))
  })

# Virtual Field: lowerPlain
# Gives a vector of lower bounds for the parameter values, whose names are NOT
# decorated by the id of this Covariance object.
#
# Args:
#   p: A (named) vector of new lower bounds (we ONLY use ones which are named,
#      and whose names match up with names of parameters.)
#
# Returns:
#   The lower bounds for the parameters for this covariance function, but with
#   names undecorated by its id.
setMethodS3("getLowerPlain", "CovarianceNoise", conflict="quiet",
  function(this, ...) {
    L <- c(this$.sigma.bounds[1])
    names(L) <- getParamNamesPlain(this)
    return (L)
})
setMethodS3("setLowerPlain", "CovarianceNoise", conflict="quiet",
  function(this, L, ...) {
    if (length(L) < 1) {
      return(invisible(this))
    }
    L.posdef <- pmax(L, 0)  # Noise cannot be negative

    # Adjust upper bounds to make way for the new values of L
    L.change <- this$PushUpperBounds(U.min=L.posdef)

    L.vals <- this$getLowerPlain()
    L.vals[names(L.change)] <- L.change[names(L.change)]
    this$.sigma.bounds[1] <- L.vals["sigma"]
    this$ClampParams(warn=TRUE)
    return (this)
  })

# Virtual Field: upperPlain
# Gives a vector of upper bounds for the parameter values, whose names are NOT
# decorated by the id of this Covariance object.
#
# Args:
#   U: A (named) vector of new upper bounds (we ONLY use ones which are named,
#      and whose names match up with names of parameters.)
#
# Returns:
#   The upper bounds for the parameters for this covariance function, but with
#   names undecorated by its id.
setMethodS3("getUpperPlain", "CovarianceNoise", conflict="quiet",
  function(this, ...) {
    U <- c(this$.sigma.bounds[2])
    names(U) <- getParamNamesPlain(this)
    return (U)
  })
setMethodS3("setUpperPlain", "CovarianceNoise", conflict="quiet",
  function(this, U, ...) {
    if (length(U) < 1) {
      return(invisible(this))
    }
    U.posdef <- pmax(U, 0)  # SE has no possibly-negative parameters

    # Adjust lower bounds to make way for the new values of U
    U.change <- this$PushLowerBounds(L.max=U.posdef)

    U.vals <- this$getUpperPlain()
    U.vals[names(U.change)] <- U.change[names(U.change)]
    this$.sigma.bounds[2] <- U.vals["sigma"]
    this$ClampParams(warn=TRUE)
    return (this)
  })

#-------------------------------------------------------------------------------
# (CovarianceNoise) PUBLIC METHODS:

setMethodS3("K.specific", "CovarianceNoise", conflict="quiet",
  function(this, X, X.out=NA, ...) {
    # Calculates a covariance matrix for the noise covariance specifically.
    #
    # Args:
    #   X:  X-values for the input points (i.e., where we have data)
    #   X.out:  X-values for the points desired to predict
    #
    # Returns:
    #   The covariance matrix between 'X' and 'X.out', based on the parameter
    #   values in 'this'.
    if (is.na(X.out)) {
      K <- (this$getParamsPlain()["sigma"] ^ 2) * diag(NumPoints(X))
    } else {
      K <- matrix(0, nrow=NumPoints(X.out), ncol=NumPoints(X))
    }
    return (K)
  })

setMethodS3("KDerivImplementation", "CovarianceNoise", conflict="quiet",
  function(this, d, param, ...) {
    # Calculate the element-wise derivative of KInIn, with respect to the
    # parameter whose (plain) name is 'param'.
    #
    # Args:
    #   d:  The Dataset whose X-values determine KInIn.
    #   param:  The (plain) name of the parameter with respect to which we're
    #      differentiating .
    #
    # Returns:
    #   A matrix whose elements are the derivatives of the corresponding
    #   elements in KInIn, with respect to the parameter 'param'.
    if (param == "sigma") {
      K.deriv <- 2 * this$paramsPlain["sigma"] * diag(d$n)
    } else {
      K.deriv <- matrix(0, nrow=d$n, ncol=d$n)
    }
    return (K.deriv)
  })

setMethodS3("Variance", "CovarianceNoise", conflict="quiet",
  # Calculate the Noise variance of the points at X.
  #
  # Args:
  #   X:  The points we want to know the noise variance at.
  #
  # Returns:
  #   A numeric vector of the same length as X, with the corresponding noise
  #   variance.
  function(this, X, ...) {
    return (rep(this$.sigma ^ 2, NumPoints(X)))
  })

################################################################################
# SUBCLASS:                      CovarianceSE
#
# Straight-up squared-exponential covariance.  Governed by two parameters: a
# horizontal and a vertical lengthscale.
#
# Methods:
#   KDerivImplementation:  The element-wise derivative of KInIn.

setConstructorS3("CovarianceSE", function(..., id="SE",
    ell=NA, sigma.f=NA, ell.bounds=NA, sigma.f.bounds=NA) {
    # Constructs a CovarianceSE object with the given values and boundaries.
    #
    # Args:
    #   ell: A characteristic ('horizontal') lengthscale over which function
    #        values are correlated.
    #   sigma.f: The "vertical" lengthscale.
    #   ell.bounds: Range of acceptable values of ell
    #   sigma.f.bounds: Range of acceptable values of sigma.f
    #
    # Returns:
    #   CovarianceSE object with the given parameter values.
    ell.good <- InitializeBoundedQuantity(ok.range=c(0, Inf),
      quantity=ell, bounds=ell.bounds, logspace=TRUE)
    sigma.f.good <- InitializeBoundedQuantity(ok.range=c(0, Inf),
      quantity=sigma.f, bounds=sigma.f.bounds, logspace=TRUE)

    # Construct the CovarianceSE object:
    extend(Covariance(..., id=id), "CovarianceSE",
      .ell            = ell.good$quantity,
      .ell.bounds     = ell.good$bounds,
      .sigma.f        = sigma.f.good$quantity,
      .sigma.f.bounds = sigma.f.good$bounds)
})

# Virtual Field: logspaceNames (read-only)
# A vector of names indicating which parameters should be optimized in
# logspace.
#
# Returns:
#   Names of parameters to be optimized in logspace.
setMethodS3("getLogspaceNames", "CovarianceSE", conflict="quiet",
  function(this, ...) {
    return (c("ell", "sigma.f"))
  })

# Virtual Field: paramNamesPlain (read-only)
# Gives the "basenames" (i.e. names undecorated by the id string) of the
# parameters.
#
# Returns:
#   The basenames of the parameters.
setMethodS3("getParamNamesPlain", "CovarianceSE", conflict="quiet",
  function(this, ...) {
    return (c("ell", "sigma.f"))
  })

# Virtual Field: paramsPlain
# Gives a vector of parameter values, whose names are NOT decorated by the id
# of this Covariance object.
#
# SEE ALSO: paramsPlain for superclass Covariance
#
# Args:
#   p: A (named) vector of new parameter values (we ONLY use ones which are
#      named, and whose names match up with names of parameters.)
#
# Returns:
#   The parameters for this covariance function, but with names undecorated by
#   its id.
setMethodS3("getParamsPlain", "CovarianceSE", conflict="quiet",
  function(this, ...) {
    p <- c(this$.ell, this$.sigma.f)
    names(p) <- getParamNamesPlain(this)
    return (p)
})
setMethodS3("paramsPlainImplementation", "CovarianceSE", conflict="quiet",
  private=TRUE,
  function(this, p, ...) {
    # Sets any values in 'p' which match our parameter names, without worrying
    # about lower and upper bounds (the function which calls this has the job
    # of worrying about these!).
    p.old <- this$getParamsPlain()
    to.change <- names(p)[which(names(p) %in% names(p.old))]
    p.old[to.change] <- p[to.change]
    this$.ell <- p.old["ell"]
    this$.sigma.f <- p.old["sigma.f"]
    return (invisible(this))
  })

# Virtual Field: lowerPlain
# Gives a vector of lower bounds for the parameter values, whose names are NOT
# decorated by the id of this Covariance object.
#
# Args:
#   p: A (named) vector of new parameter values (we ONLY use ones which are
#      named, and whose names match up with names of parameters.)
#
# Returns:
#   The lower bounds for the parameters for this covariance function, but with
#   names undecorated by its id.
setMethodS3("getLowerPlain", "CovarianceSE", conflict="quiet",
  function(this, ...) {
    L <- c(this$.ell.bounds[1], this$.sigma.f.bounds[1])
    names(L) <- getParamNamesPlain(this)
    return (L)
})
setMethodS3("setLowerPlain", "CovarianceSE", conflict="quiet",
  function(this, L, ...) {
    if (length(L) < 1) {
      return(invisible(this))
    }

    L.posdef <- pmax(L, 0)  # SE has no possibly-negative parameters

    # Adjust upper bounds to make way for the new values of L
    L.change <- this$PushUpperBounds(U.min=L.posdef)

    L.vals <- this$getLowerPlain()
    L.vals[names(L.change)] <- L.change[names(L.change)]
    this$.ell.bounds[1] <- L.vals["ell"]
    this$.sigma.f.bounds[1] <- L.vals["sigma.f"]
    this$ClampParams(warn=TRUE)
    return (this)
  })

# Virtual Field: upperPlain
# Gives a vector of upper bounds for the parameter values, whose names are NOT
# decorated by the id of this Covariance object.
#
# Args:
#   U: A (named) vector of new upper bounds (we ONLY use ones which are named,
#      and whose names match up with names of parameters.)
#
# Returns:
#   The upper bounds for the parameters for this covariance function, but with
#   names undecorated by its id.
setMethodS3("getUpperPlain", "CovarianceSE", conflict="quiet",
  function(this, ...) {
    U <- c(this$.ell.bounds[2], this$.sigma.f.bounds[2])
    names(U) <- getParamNamesPlain(this)
    return (U)
  })
setMethodS3("setUpperPlain", "CovarianceSE", conflict="quiet",
  function(this, U, ...) {
    if (length(U) < 1) {
      return(invisible(this))
    }
    U.posdef <- pmax(U, 0)  # SE has no possibly-negative parameters

    # Adjust lower bounds to make way for the new values of U
    U.change <- this$PushLowerBounds(L.max=U.posdef)

    U.vals <- this$getUpperPlain()
    U.vals[names(U.change)] <- U.change[names(U.change)]
    this$.ell.bounds[2] <- U.vals["ell"]
    this$.sigma.f.bounds[2] <- U.vals["sigma.f"]
    this$ClampParams(warn=TRUE)
    return (this)
  })

setMethodS3("K.specific", "CovarianceSE", conflict="quiet",
  function(this, X, X.out=X, ...) {
    # Calculates a covariance matrix for the SE covariance specifically.
    #
    # Args:
    #   X:  X-values for the input points (i.e., where we have data)
    #   X.out:  X-values for the points desired to predict
    #
    # Returns:
    #   The covariance matrix between 'X' and 'X.out', based on the parameter
    #   values in 'this'.
    X.dist <- DistanceMatrix(X=X, X.out=X.out)
    p <- this$getParamsPlain()
    return (p["sigma.f"] ^ 2 * exp(-0.5 * (X.dist / p["ell"]) ^ 2))
  })

#-------------------------------------------------------------------------------
# (CovarianceSE) PUBLIC METHODS:

setMethodS3("KDerivImplementation", "CovarianceSE", conflict="quiet",
  function(this, d, param, ...) {
    # Calculate the element-wise derivative of KInIn, with respect to the
    # parameter whose (plain) name is 'param'.
    #
    # Args:
    #   d:  The Dataset whose X-values determine KInIn.
    #   param:  The (plain) name of the parameter with respect to which we're
    #      differentiating .
    #
    # Returns:
    #   A matrix whose elements are the derivatives of the corresponding
    #   elements in KInIn, with respect to the parameter 'param'.
    p <- this$paramsPlain
    if (param == "ell") {
      K.deriv <- this$KInIn(d=d) * (DistanceMatrix(X=d$X) ^ 2) / (p["ell"] ^ 3)
    } else if (param == "sigma.f") {
      K.deriv <- 2 * this$KInIn(d=d) / p["sigma.f"]
    } else {
      K.deriv <- matrix(0, nrow=d$n, ncol=d$n)
    }
    return (K.deriv)
  })

setMethodS3("Variance", "CovarianceSE", conflict="quiet",
  # Calculate the SE variance of the points at X.
  #
  # Args:
  #   X:  The points we want to know the SE variance at.
  #
  # Returns:
  #   A numeric vector of the same length as X, with the corresponding SE
  #   variance.
  function(this, X, ...) {
    return (rep(this$.sigma.f ^ 2, NumPoints(X)))
  })

################################################################################
# SUBCLASS:                  CovarianceSELocalized
#
# Squared-exponential covariance (like CovarianceSE), but localized to a
# particular region by masking off sigma.f.
#
# Methods:
#   KDerivImplementation:  The element-wise derivative of KInIn.

setConstructorS3("CovarianceSELocalized", function(..., id="SE",
    ell=NA, sigma.f=NA, X.L=NA, X.R=NA, ell.bounds=NA, sigma.f.bounds=NA,
    X.L.bounds=NA, X.R.bounds=NA) {
    # Constructs a CovarianceSELocalized object with the given values and
    # boundaries.
    #
    # Args:
    #   ell:  A characteristic ('horizontal') lengthscale over which function
    #      values are correlated.
    #   sigma.f:  The "vertical" lengthscale.
    #   X.L:  The left boundary of the region
    #   X.R:  The right boundary of the region
    #   ell.bounds:  Range of acceptable values of ell
    #   sigma.f.bounds:  Range of acceptable values of sigma.f
    #   X.L.bounds:  Range of acceptable values for X.L
    #   X.R.bounds:  Range of acceptable values for X.R
    #
    # Returns:
    #   CovarianceSELocalized object with the given parameter values.

    ell.good <- InitializeBoundedQuantity(ok.range=c(0, Inf),
      quantity=ell, bounds=ell.bounds, logspace=TRUE)
    sigma.f.good <- InitializeBoundedQuantity(ok.range=c(0, Inf),
      quantity=sigma.f, bounds=sigma.f.bounds, logspace=TRUE)
    X.L.good <- InitializeBoundedQuantity(quantity=X.L, bounds=X.L.bounds)
    X.R.good <- InitializeBoundedQuantity(quantity=X.R, bounds=X.R.bounds)

    # Construct the CovarianceSELocalized object:
    extend(Covariance(..., id=id), "CovarianceSELocalized",
      .ell            = ell.good$quantity,
      .ell.bounds     = ell.good$bounds,
      .sigma.f        = sigma.f.good$quantity,
      .sigma.f.bounds = sigma.f.good$bounds,
      .X.L            = X.L.good$quantity,
      .X.L.bounds     = X.L.good$bounds,
      .X.R            = X.R.good$quantity,
      .X.R.bounds     = X.R.good$bounds)
})

#-------------------------------------------------------------------------------
# (CovarianceSELocalized) PUBLIC "HELPER" FUNCTIONS:

erf <- function(x) {
  # The standard function 'erf', based on the CDF of the standard normal
  # distribution.
  #
  # Args:
  #   x:  The x-value for erf(x).
  #
  # Returns:
  #   erf(x)
  2 * pnorm(sqrt(2) * x) - 1
}

standard.mask <- function(x) {
  # A sigmoid mask which goes from (-Inf, -1) to (+Inf, +1), with a slope of 1
  # in the middle.
  #
  # Args:
  #   x:  The coordinate to evaluate the mask at.
  #
  # Returns:
  #   The standard sigmoid evaluated at x.
  return (as.vector(erf(sqrt(pi) * x * 0.5)))
}

standard.mask.deriv <- function(x) {
  # Derivative of the standard sigmoid mask 'standard.mask'.
  #
  # Args:
  #   x:  The coordinate to evaluate the mask derivative at.
  #
  # Returns:
  #   The derivative of the standard sigmoid, evaluated at x.
  return (as.vector(exp(-0.25 * pi * (x ^ 2))))
}

localize.mask <- function(X, X.L, X.R, ell) {
  # Mask values in X according to whether they're between X.L and X.R, with
  # transition lengthscale ell.
  #
  # Args:
  #   X:  The values to mask.
  #   X.L:  The left boundary of the region.
  #   X.L:  The right boundary of the region.
  #   ell:  The transition lengthscale at the edges of the region.

  upper.mask <- standard.mask((X - X.L) / ell)
  lower.mask <- standard.mask((X - X.R) / ell)
  return (as.vector(0.5 * (upper.mask - lower.mask)))
}

#-------------------------------------------------------------------------------
# (CovarianceSELocalized) PUBLIC VIRTUAL FIELDS:

# Virtual Field: logspaceNames (read-only)
# A vector of names indicating which parameters should be optimized in
# logspace.
#
# Returns:
#   Names of parameters to be optimized in logspace.
setMethodS3("getLogspaceNames", "CovarianceSELocalized", conflict="quiet",
  function(this, ...) {
    return (c("ell", "sigma.f"))
  })

# Virtual Field: paramNamesPlain (read-only)
# Gives the "basenames" (i.e. names undecorated by the id string) of the
# parameters.
#
# Returns:
#   The basenames of the parameters.
setMethodS3("getParamNamesPlain", "CovarianceSELocalized", conflict="quiet",
  function(this, ...) {
    return (c("ell", "sigma.f", "X.L", "X.R"))
  })

# Virtual Field: paramsPlain
# Gives a vector of parameter values, whose names are NOT decorated by the id
# of this Covariance object.
#
# SEE ALSO: paramsPlain for superclass Covariance
#
# Args:
#   p: A (named) vector of new parameter values (we ONLY use ones which are
#      named, and whose names match up with names of parameters.)
#
# Returns:
#   The parameters for this covariance function, but with names undecorated by
#   its id.
setMethodS3("getParamsPlain", "CovarianceSELocalized", conflict="quiet",
  function(this, ...) {
    p <- c(this$.ell, this$.sigma.f, this$.X.L, this$.X.R)
    names(p) <- getParamNamesPlain(this)
    return (p)
})
setMethodS3("paramsPlainImplementation", "CovarianceSELocalized", conflict="quiet",
  private=TRUE,
  function(this, p, ...) {
    # Sets any values in 'p' which match our parameter names, without worrying
    # about lower and upper bounds (the function which calls this has the job
    # of worrying about these!).
    p.old <- this$getParamsPlain()
    to.change <- names(p)[which(names(p) %in% names(p.old))]
    p.old[to.change] <- p[to.change]
    this$.ell <- p.old["ell"]
    this$.sigma.f <- p.old["sigma.f"]
    this$.X.L <- p.old["X.L"]
    this$.X.R <- p.old["X.R"]
    return (invisible(this))
  })

# Virtual Field: lowerPlain
# Gives a vector of lower bounds for the parameter values, whose names are NOT
# decorated by the id of this Covariance object.
#
# Args:
#   p: A (named) vector of new parameter values (we ONLY use ones which are
#      named, and whose names match up with names of parameters.)
#
# Returns:
#   The lower bounds for the parameters for this covariance function, but with
#   names undecorated by its id.
setMethodS3("getLowerPlain", "CovarianceSELocalized", conflict="quiet",
  function(this, ...) {
    L <- c(this$.ell.bounds[1], this$.sigma.f.bounds[1], this$.X.L.bounds[1],
      this$.X.R.bounds[1])
    names(L) <- getParamNamesPlain(this)
    return (L)
})
setMethodS3("setLowerPlain", "CovarianceSELocalized", conflict="quiet",
  function(this, L, ...) {
    if (length(L) < 1) {
      return(invisible(this))
    }

    posdef.names <- c('ell', 'sigma.f')
    L[posdef.names] <- pmax(L[posdef.names], 0)

    # Adjust upper bounds to make way for the new values of L
    L.change <- this$PushUpperBounds(U.min=L)

    L.vals <- this$getLowerPlain()
    L.vals[names(L.change)] <- L.change[names(L.change)]
    this$.ell.bounds[1] <- L.vals["ell"]
    this$.sigma.f.bounds[1] <- L.vals["sigma.f"]
    this$.X.L.bounds[1] <- L.vals["X.L"]
    this$.X.R.bounds[1] <- L.vals["X.R"]
    this$ClampParams(warn=TRUE)
    return (this)
  })

# Virtual Field: upperPlain
# Gives a vector of upper bounds for the parameter values, whose names are NOT
# decorated by the id of this Covariance object.
#
# Args:
#   U: A (named) vector of new upper bounds (we ONLY use ones which are named,
#      and whose names match up with names of parameters.)
#
# Returns:
#   The upper bounds for the parameters for this covariance function, but with
#   names undecorated by its id.
setMethodS3("getUpperPlain", "CovarianceSELocalized", conflict="quiet",
  function(this, ...) {
    U <- c(this$.ell.bounds[2], this$.sigma.f.bounds[2], this$.X.L.bounds[2],
      this$.X.R.bounds[2])
    names(U) <- getParamNamesPlain(this)
    return (U)
  })
setMethodS3("setUpperPlain", "CovarianceSELocalized", conflict="quiet",
  function(this, U, ...) {
    if (length(U) < 1) {
      return(invisible(this))
    }

    posdef.names <- c('ell', 'sigma.f')
    U[posdef.names] <- pmax(U[posdef.names], 0)

    # Adjust lower bounds to make way for the new values of U
    U.change <- this$PushLowerBounds(L.max=U)

    U.vals <- this$getUpperPlain()
    U.vals[names(U.change)] <- U.change[names(U.change)]
    this$.ell.bounds[2] <- U.vals["ell"]
    this$.sigma.f.bounds[2] <- U.vals["sigma.f"]
    this$.X.L.bounds[2] <- U.vals["X.L"]
    this$.X.R.bounds[2] <- U.vals["X.R"]
    this$ClampParams(warn=TRUE)
    return (this)
  })

setMethodS3("K.specific", "CovarianceSELocalized", conflict="quiet",
  function(this, X, X.out=X, ...) {
    # Calculates a covariance matrix for the SELocalized covariance
    # specifically.
    #
    # Args:
    #   X:  X-values for the input points (i.e., where we have data)
    #   X.out:  X-values for the points desired to predict
    #
    # Returns:
    #   The covariance matrix between 'X' and 'X.out', based on the parameter
    #   values in 'this'.
    X.dist <<- DistanceMatrix(X=X, X.out=X.out)
    p <- this$getParamsPlain()
    mask <<- localize.mask(X=X, X.L=p["X.L"], X.R=p["X.R"], ell=p["ell"])
    mask.out <<- localize.mask(X=X.out, X.L=p["X.L"], X.R=p["X.R"], ell=p["ell"])
    return (outer(mask.out, mask) *
      (p["sigma.f"] ^ 2) * exp(-0.5 * (X.dist / p["ell"]) ^ 2))
  })

#-------------------------------------------------------------------------------
# (CovarianceSELocalized) PUBLIC METHODS:

setMethodS3("KDerivImplementation", "CovarianceSELocalized", conflict="quiet",
  function(this, d, param, ...) {
    # Calculate the element-wise derivative of KInIn, with respect to the
    # parameter whose (plain) name is 'param'.
    #
    # Args:
    #   d:  The Dataset whose X-values determine KInIn.
    #   param:  The (plain) name of the parameter with respect to which we're
    #      differentiating .
    #
    # Returns:
    #   A matrix whose elements are the derivatives of the corresponding
    #   elements in KInIn, with respect to the parameter 'param'.
    p <- this$paramsPlain
    K.unmasked <- K.specific.CovarianceSE(this=this, X=d$X)
    dist.L <- as.vector((d$X - p["X.L"]) / p["ell"])
    dist.R <- as.vector((d$X - p["X.R"]) / p["ell"])
    mask <- localize.mask(X=d$X, X.L=p["X.L"], X.R=p["X.R"], ell=p["ell"])
    if (param == "ell") {
      d.mask.d.ell <- -(0.5 / p["ell"]) * (
        dist.L * standard.mask.deriv(dist.L) -
        dist.R * standard.mask.deriv(dist.R))
      mask.deriv.matrix <<- outer(mask, d.mask.d.ell)
      mask.deriv.matrix <- mask.deriv.matrix + t(mask.deriv.matrix)
      K.deriv <- (mask.deriv.matrix * K.unmasked
        + this$KInIn(d=d) * (DistanceMatrix(X=d$X) ^ 2) / (p["ell"] ^ 3))
    } else if (param == "sigma.f") {
      K.deriv <- 2 * this$KInIn(d=d) / p["sigma.f"]
    } else if (param == "X.L") {
      d.mask.d.X.L <- -(0.5 / p["ell"]) * standard.mask.deriv(dist.L)
      mask.deriv.matrix <- outer(mask, d.mask.d.X.L)
      mask.deriv.matrix <- mask.deriv.matrix + t(mask.deriv.matrix)
      K.deriv <- mask.deriv.matrix * K.unmasked
    } else if (param == "X.R") {
      d.mask.d.X.R <-  (0.5 / p["ell"]) * standard.mask.deriv(dist.R)
      mask.deriv.matrix <- outer(mask, d.mask.d.X.R)
      mask.deriv.matrix <- mask.deriv.matrix + t(mask.deriv.matrix)
      K.deriv <- mask.deriv.matrix * K.unmasked
    } else {
      K.deriv <- matrix(0, nrow=d$n, ncol=d$n)
    }
    return (K.deriv)
  })

setMethodS3("Variance", "CovarianceSELocalized", conflict="quiet",
  # Calculate the SE variance of the points at X.
  #
  # Args:
  #   X:  The points we want to know the SE variance at.
  #
  # Returns:
  #   A numeric vector of the same length as X, with the corresponding SE
  #   variance.
  function(this, X, ...) {
    p <- this$paramsPlain
    mask <- localize.mask(X=X, X.L=p["X.L"], X.R=p["X.R"], ell=p["ell"])
    return ((p["sigma.f"] * mask) ^ 2)
  })

################################################################################
# SUBCLASS:                  CovarianceSEVaryingEll
#
# Straight-up squared-exponential covariance with ell depending on X (and
# constant sigma.f).  ell is defined by the spline interpolation of a set of
# (X, ell) pairs.
#
# Methods:
#   KDerivImplementation:  The element-wise derivative of KInIn.
#
# NOTES:
#   For now, the parameters are not mutable.

setConstructorS3("CovarianceSEVaryingEll", function(..., id="SEVaryingEll",
    X.ell=NA, ell=NA, sigma.f=NA) {
    # Constructs a CovarianceSEVaryingEll object with the given parameter
    # values.
    #
    # Args:
    #   X.ell:  numeric vector of X-values corresponding to 'ell'
    #   ell: Samples of ell(X) (a characteristic lengthscale over which
    #      function values are correlated) at the locations in X.ell.
    #   sigma.f: The "vertical" lengthscale.
    #
    # Returns:
    #   CovarianceSEVaryingEll object with the given parameter values.

    # Construct the CovarianceSEVaryingEll object:
    extend(Covariance(..., id=id), "CovarianceSEVaryingEll",
      .X.ell   = X.ell,
      .ell     = ell,
      .sigma.f = sigma.f)
})

# Virtual Field: logspaceNames (read-only)
# A vector of names indicating which parameters should be optimized in
# logspace.
#
# Returns:
#   Names of parameters to be optimized in logspace.
setMethodS3("getLogspaceNames", "CovarianceSEVaryingEll", conflict="quiet",
  function(this, ...) {
    p.names <- this$getParamNamesPlain()
    return (p.names[-grep("X.ell", p.names)])
  })

# Virtual Field: paramNamesPlain (read-only)
# Gives the "basenames" (i.e. names undecorated by the id string) of the
# parameters.
#
# Returns:
#   The basenames of the parameters.
setMethodS3("getParamNamesPlain", "CovarianceSEVaryingEll", conflict="quiet",
  function(this, ...) {
    p.names <- c(
      paste(sep='', 'sigma.f.', 1:length(this$.sigma.f)),
      paste(sep='', 'X.ell.', 1:length(this$.X.ell)),
      paste(sep='', 'ell.', 1:length(this$.ell)))
    return (p.names)
  })

# Virtual Field: paramsPlain
# Gives a vector of parameter values, whose names are NOT decorated by the id
# of this Covariance object.
#
# SEE ALSO: paramsPlain for superclass Covariance
#
# Args:
#   p: A (named) vector of new parameter values (we ONLY use ones which are
#      named, and whose names match up with names of parameters.)
#
# Returns:
#   The parameters for this covariance function, but with names undecorated by
#   its id.
setMethodS3("getParamsPlain", "CovarianceSEVaryingEll", conflict="quiet",
  function(this, ...) {
    p <- c(this$.sigma.f, this$.X.ell, this$.ell)
    n <- getParamNamesPlain(this)
    names(p) <- n
    return (p)
})
setMethodS3("paramsPlainImplementation", "CovarianceSEVaryingEll", conflict="quiet",
  private=TRUE,
  function(this, p, ...) {
    warning(
      "At present, param values in 'CovarianceSEVaryingEll' are not mutable.\n")
    ## Sets any values in 'p' which match our parameter names, without worrying
    ## about lower and upper bounds (the function which calls this has the job
    ## of worrying about these!).
    #p.old <- this$getParamsPlain()
    #to.change <- names(p)[which(names(p) %in% names(p.old))]
    #p.old[to.change] <- p[to.change]
    return (invisible(this))
  })

# Virtual Field: lowerPlain
# Gives a vector of lower bounds for the parameter values, whose names are NOT
# decorated by the id of this Covariance object.
#
# Args:
#   p: A (named) vector of new parameter values (we ONLY use ones which are
#      named, and whose names match up with names of parameters.)
#
# Returns:
#   The lower bounds for the parameters for this covariance function, but with
#   names undecorated by its id.
setMethodS3("getLowerPlain", "CovarianceSEVaryingEll", conflict="quiet",
  function(this, ...) {
    return (this$getParamsPlain(...))
})
setMethodS3("setLowerPlain", "CovarianceSEVaryingEll", conflict="quiet",
  function(this, L, ...) {
    warning(
      "At present, lower bounds in 'CovarianceSEVaryingEll' are not mutable.\n")
    return (invisible(this))
  })

# Virtual Field: upperPlain
# Gives a vector of upper bounds for the parameter values, whose names are NOT
# decorated by the id of this Covariance object.
#
# Args:
#   U: A (named) vector of new upper bounds (we ONLY use ones which are named,
#      and whose names match up with names of parameters.)
#
# Returns:
#   The upper bounds for the parameters for this covariance function, but with
#   names undecorated by its id.
setMethodS3("getUpperPlain", "CovarianceSEVaryingEll", conflict="quiet",
  function(this, ...) {
    return (this$getParamsPlain(...))
  })
setMethodS3("setUpperPlain", "CovarianceSEVaryingEll", conflict="quiet",
  function(this, U, ...) {
    warning(
      "At present, lower bounds in 'CovarianceSEVaryingEll' are not mutable.\n")
    return (invisible(this))
  })

#-------------------------------------------------------------------------------
# (CovarianceSEVaryingEll) PUBLIC METHODS:

setMethodS3("ell", "CovarianceSEVaryingEll", conflict="quiet",
  function(this, X, ...) {
    # Evaluate ell(X) at the applied X-points (by spline-interpolating)
    return (exp(spline(xout=X, x=this$.X.ell, y=log(this$.ell))$y))
    return (spline(xout=X, x=this$.X.ell, y=this$.ell)$y)
  })

setMethodS3("sigma.f", "CovarianceSEVaryingEll", conflict="quiet",
  function(this, X, ...) {
    # This code repeats a single-valued sigma if necessary, so it can handle
    # both single-sigma and varying-sigma.
    sigma.vals <- data.frame(x=this$.X.ell, y=this$.sigma.f)
    # Evaluate sigma.f(X) at the applied X-points (by spline-interpolating)
    sigma.at.X <- with(sigma.vals, spline(x=x, y=y, xout=X)$y)
    return (sigma.at.X)
  })

setMethodS3("K.specific", "CovarianceSEVaryingEll", conflict="quiet",
  function(this, X, X.out=X, ...) {
    # Calculates a covariance matrix for the SE covariance specifically.
    #
    # Args:
    #   X:  X-values for the input points (i.e., where we have data)
    #   X.out:  X-values for the points desired to predict
    #
    # Returns:
    #   The covariance matrix between 'X' and 'X.out', based on the parameter
    #   values in 'this'.
    X.dist <- DistanceMatrix(X=X, X.out=X.out)
    # Check and warn if user wants to extrapolate
    if (min(X.out) < min(X) || max(X.out) > max(X)) {
      warning(
        "Extrapolation is highly inadvisable with CovarianceSEVaryingEll!\n")
    }
    # Get params and calculate the matrix
    ell     <- this$ell(X=X    )
    ell.out <- this$ell(X=X.out)
    sigma     <- this$sigma.f(X=X    )
    sigma.out <- this$sigma.f(X=X.out)
    sum.ell.sq <- outer(ell.out ^ 2, ell ^ 2, '+')
    K <- (outer(sigma.out, sigma) 
      * sqrt(2.0 * outer(ell.out, ell) / sum.ell.sq)
      * exp(-(X.dist ^ 2) / sum.ell.sq))
    gc()
    return (K)
  })

setMethodS3("KDerivImplementation", "CovarianceSEVaryingEll", conflict="quiet",
  function(this, d, param, ...) {
    # Calculate the element-wise derivative of KInIn, with respect to the
    # parameter whose (plain) name is 'param'.
    #
    # Args:
    #   d:  The Dataset whose X-values determine KInIn.
    #   param:  The (plain) name of the parameter with respect to which we're
    #      differentiating .
    #
    # Returns:
    #   A matrix whose elements are the derivatives of the corresponding
    #   elements in KInIn, with respect to the parameter 'param'.
    warning("This function not yet implemented for CovarianceSEVaryingEll.\n")
    return (NULL)
  })

setMethodS3("Variance", "CovarianceSEVaryingEll", conflict="quiet",
  # Calculate the SE variance of the points at X.
  #
  # Args:
  #   X:  The points we want to know the SE variance at.
  #
  # Returns:
  #   A numeric vector of the same length as X, with the corresponding SE
  #   variance.
  function(this, X, ...) {
    return (rep(mean(this$.sigma.f) ^ 2, NumPoints(X)))
  })

################################################################################
# SUBCLASS:                   CovarianceSEAniso2D
#
# Anisotropic squared-exponential covariance for 2D data, with an extra
# parameter governing the angle.
#
# Methods:
#   KDerivImplementation:  The element-wise derivative of KInIn.

setConstructorS3("CovarianceSEAniso2D",
  function(..., id="Aniso2D", ell.1=NA, ell.2=NA, theta.1=NA, sigma.f=NA,
    ell.1.bounds=NA, ell.2.bounds=NA, theta.1.bounds=NA, sigma.f.bounds=NA) {
    pos.def.range <- c(0, Inf)
    ell.1.good <- InitializeBoundedQuantity(ok.range=pos.def.range,
      quantity=ell.1, bounds=ell.1.bounds, logspace=TRUE)
    ell.2.good <- InitializeBoundedQuantity(ok.range=pos.def.range,
      quantity=ell.2, bounds=ell.2.bounds, logspace=TRUE)
    theta.1.good <- InitializeBoundedQuantity(ok.range=pi * c(-1, 1),
      quantity=theta.1, bounds=theta.1.bounds)
    sigma.f.good <- InitializeBoundedQuantity(ok.range=pos.def.range,
      quantity=sigma.f, bounds=sigma.f.bounds, logspace=TRUE)

    extend(Covariance(..., id=id), "CovarianceSEAniso2D",
      .ell.1          = ell.1.good$quantity,
      .ell.1.bounds   = ell.1.good$bounds,
      .ell.2          = ell.2.good$quantity,
      .ell.2.bounds   = ell.2.good$bounds,
      .theta.1        = theta.1.good$quantity,
      .theta.1.bounds = theta.1.good$bounds,
      .sigma.f        = sigma.f.good$quantity,
      .sigma.f.bounds = sigma.f.good$bounds)
  })

#-------------------------------------------------------------------------------
# (CovarianceSEAniso2D) PUBLIC VIRTUAL FIELDS:

# Virtual Field: logspaceNames (read-only)
# A vector of names indicating which parameters should be optimized in
# logspace.
#
# Returns:
#   Names of parameters to be optimized in logspace.
setMethodS3("getLogspaceNames", "CovarianceSEAniso2D", conflict="quiet",
  function(this, ...) {
    return (c("ell.1", "ell.2", "sigma.f"))
  })

# Virtual Field: paramNamesPlain (read-only)
# Gives the "basenames" (i.e. names undecorated by the id string) of the
# parameters.
#
# Returns:
#   The basenames of the parameters.
setMethodS3("getParamNamesPlain", "CovarianceSEAniso2D", conflict="quiet",
  function(this, ...) {
    return (c("ell.1", "ell.2", "theta.1", "sigma.f"))
  })

# Virtual Field: paramsPlain
# Gives a vector of parameter values, whose names are NOT decorated by the id
# of this Covariance object.
#
# SEE ALSO: paramsPlain for superclass Covariance
#
# Args:
#   p: A (named) vector of new parameter values (we ONLY use ones which are
#      named, and whose names match up with names of parameters.)
#
# Returns:
#   The parameters for this covariance function, but with names undecorated by
#   its id.
setMethodS3("getParamsPlain", "CovarianceSEAniso2D", conflict="quiet",
  function(this, ...) {
    p <- c(this$.ell.1, this$.ell.2, this$.theta.1, this$.sigma.f)
    names(p) <- getParamNamesPlain(this)
    return (p)
})
setMethodS3("paramsPlainImplementation", "CovarianceSEAniso2D", conflict="quiet",
  private=TRUE,
  function(this, p, ...) {
    # Sets any values in 'p' which match our parameter names, without worrying
    # about lower and upper bounds (the function which calls this has the job
    # of worrying about these!).
    p.old <- this$getParamsPlain()
    to.change <- names(p)[which(names(p) %in% names(p.old))]
    p.old[to.change] <- p[to.change]
    this$.ell.1 <- p.old["ell.1"]
    this$.ell.2 <- p.old["ell.2"]
    this$.theta.1 <- p.old["theta.1"]
    this$.sigma.f <- p.old["sigma.f"]
    return (invisible(this))
  })

# Virtual Field: lowerPlain
# Gives a vector of lower bounds for the parameter values, whose names are NOT
# decorated by the id of this Covariance object.
#
# Args:
#   p: A (named) vector of new parameter values (we ONLY use ones which are
#      named, and whose names match up with names of parameters.)
#
# Returns:
#   The lower bounds for the parameters for this covariance function, but with
#   names undecorated by its id.
setMethodS3("getLowerPlain", "CovarianceSEAniso2D", conflict="quiet",
  function(this, ...) {
    L <- c(this$.ell.1.bounds[1], this$.ell.2.bounds[1],
      this$.theta.1.bounds[1], this$.sigma.f.bounds[1])
    names(L) <- getParamNamesPlain(this)
    return (L)
})
setMethodS3("setLowerPlain", "CovarianceSEAniso2D", conflict="quiet",
  function(this, L, ...) {
    # Adjust upper bounds to make way for the new values of L
    L.change <- this$PushUpperBounds(U.min=L)

    L.vals <- this$getLowerPlain()
    L.vals[names(L.change)] <- L.change[names(L.change)]
    this$.ell.1.bounds[1] <- L.vals["ell.1"]
    this$.ell.2.bounds[1] <- L.vals["ell.2"]
    this$.theta.1.bounds[1] <- L.vals["theta.1"]
    this$.sigma.f.bounds[1] <- L.vals["sigma.f"]
    this$ClampParams(warn=TRUE)
    return (this)
  })

# Virtual Field: upperPlain
# Gives a vector of upper bounds for the parameter values, whose names are NOT
# decorated by the id of this Covariance object.
#
# Args:
#   U: A (named) vector of new upper bounds (we ONLY use ones which are named,
#      and whose names match up with names of parameters.)
#
# Returns:
#   The upper bounds for the parameters for this covariance function, but with
#   names undecorated by its id.
setMethodS3("getUpperPlain", "CovarianceSEAniso2D", conflict="quiet",
  function(this, ...) {
    U <- c(this$.ell.1.bounds[2], this$.ell.2.bounds[2],
      this$.theta.1.bounds[2], this$.sigma.f.bounds[2])
    names(U) <- getParamNamesPlain(this)
    return (U)
  })
setMethodS3("setUpperPlain", "CovarianceSEAniso2D", conflict="quiet",
  function(this, U, ...) {
    # Adjust lower bounds to make way for the new values of U
    U.change <- this$PushLowerBounds(L.max=U)

    U.vals <- this$getUpperPlain()
    U.vals[names(U.change)] <- U.change[names(U.change)]
    this$.ell.1.bounds[2] <- U.vals["ell.1"]
    this$.ell.2.bounds[2] <- U.vals["ell.2"]
    this$.theta.1.bounds[2] <- U.vals["theta.1"]
    this$.sigma.f.bounds[2] <- U.vals["sigma.f"]
    this$ClampParams(warn=TRUE)
    return (this)
  })

setMethodS3("K.specific", "CovarianceSEAniso2D", conflict="quiet",
  function(this, X, X.out=X, ...) {
    # Calculates a covariance matrix for the SE covariance specifically.
    #
    # Args:
    #   X:  X-values for the input points (i.e., where we have data)
    #   X.out:  X-values for the points desired to predict
    #
    # Returns:
    #   The covariance matrix between 'X' and 'X.out', based on the parameter
    #   values in 'this'.
    p <- this$getParamsPlain()
    ell <- c(p["ell.1"], p["ell.2"])

    # Rotate and scale the X-values
    cos.1 <- cos(p["theta.1"])
    sin.1 <- sin(p["theta.1"])
    R <- matrix(c(cos.1, sin.1, -sin.1, cos.1), ncol=2)
    X.rot     <- (X     %*% R) / matrix(
      rep(ell, each=NumPoints(X    )), ncol=2)
    X.out.rot <- (X.out %*% R) / matrix(
      rep(ell, each=NumPoints(X.out)), ncol=2)

    X.dist <- DistanceMatrix(X=X.rot, X.out=X.out.rot)
    return (p["sigma.f"] ^ 2 * exp(-0.5 * (X.dist ^ 2)))
  })

#-------------------------------------------------------------------------------
# (CovarianceSEAniso2D) PUBLIC METHODS:

setMethodS3("KDerivImplementation", "CovarianceSEAniso2D", conflict="quiet",
  function(this, d, param, ...) {
    # Calculate the element-wise derivative of KInIn, with respect to the
    # parameter whose (plain) name is 'param'.
    #
    # Args:
    #   d:  The Dataset whose X-values determine KInIn.
    #   param:  The (plain) name of the parameter with respect to which we're
    #      differentiating .
    #
    # Returns:
    #   A matrix whose elements are the derivatives of the corresponding
    #   elements in KInIn, with respect to the parameter 'param'.
    Delta.1 <- outer(d$X[, 1], d$X[, 1], '-')
    Delta.2 <- outer(d$X[, 2], d$X[, 2], '-')
    p <- this$paramsPlain
    c.theta <- cos(p["theta.1"])
    s.theta <- sin(p["theta.1"])
    if (param == "ell.1") {
      K.deriv <- (this$KInIn(d=d) / p["ell.1"]) * (
        (c.theta * Delta.1 - s.theta * Delta.2) / p["ell.1"]) ^ 2
    } else if (param == "ell.2") {
      K.deriv <- (this$KInIn(d=d) / p["ell.2"]) * (
        (s.theta * Delta.1 + c.theta * Delta.2) / p["ell.2"]) ^ 2
    } else if (param == "sigma.f") {
      K.deriv <- 2 * this$KInIn(d=d) / p["sigma.f"]
    } else if (param == "theta.1") {
      K.deriv <- (this$KInIn(d=d) 
        * (1 / (p["ell.1"] ^ 2) - 1 / (p["ell.2"] ^ 2)) 
        * (0.5 * sin(2 * p["theta.1"]) * (Delta.1 ^ 2 - Delta.2 ^ 2)
          + cos(2 * p["theta.1"]) * Delta.1 * Delta.2)
        )
    } else {
      K.deriv <- matrix(0, nrow=d$n, ncol=d$n)
    }
    return (K.deriv)
  })

setMethodS3("Variance", "CovarianceSEAniso2D", conflict="quiet",
  # Calculate the SE variance of the points at X.
  #
  # Args:
  #   X:  The points we want to know the SE variance at.
  #
  # Returns:
  #   A numeric vector of the same length as X, with the corresponding SE
  #   variance.
  function(this, X, ...) {
    return (rep(this$.sigma.f ^ 2, NumPoints(X)))
  })

################################################################################
# CLASS:                             Model
#
# A collection of covariance structures which can be trained on datasets and
# make predictions.
#
# Helper Functions:
#   LogML:  The logarithm of the (M)arginal (L)ikelihood for this model, given
#      a Dataset.
#   GradLogML:  The gradient (w.r.t. the parameter values) of the logarithm of
#      the (M)arginal (L)ikelihood for this model, given a Dataset.
#
# Virtual Fields (R = read only):
# (R) contributionIds:  The id strings of the Covariance objects which
#        contribute to this Model.
# ( ) id:  A character string used to identify this model.
# ( ) lower:  A named vector of lower bounds on parameters governing this
#     Model.
# ( ) params:  A named vector of parameters governing this Model.
# (R) signalIds:  The ids of contributions which are considered to be "signal",
#        i.e., not noise.
# ( ) upper:  A named vector of upper bounds on parameters governing this
#     Model.
# (R) varyingParamNames:  Names of parameters with nonzero range between upper
#     and lower bounds.
#
# Methods:
#   AddCovariance:  Add a new Covariance object to this Model.
#   clone:  Deep clone of this Model object (clones Covariance's as well).
#   Forget:  Clear all large matrices to save memory.
#   Freeze:  Make all parameter values constant.
#   L:  Lower-triangular Cholesky decomposition of the covariance matrix
#      (useful for generating random draws).
#   NamedCovariance:  Retrieve a clone of the first contributing Covariance
#      object with the given id.
#   PlotBubblingSurfaces2D:  Plot smoothly varying random surfaces to visualize
#      the uncertainty.
#   PlotCovariance2D:  Plot the covariance of points in 2D with respect to
#      another point; also plot the derivatives w.r.t. each parameter.
#   PosteriorInterval:  bounds on the uncertainty (in terms of the standard
#      deviation) in the prediction at a given point.
#   PosteriorMean:  The optimal prediction of the underlying function at a set
#      of points.
#   PosteriorStandardDeviation:  pointwise sigma (for *transformed* values) at
#      each prediction point.
#   PredictionMatrix:  A matrix relating function values at the output points
#      to function values at the input points.
#   print:  Prettied-up summary of this Model object.
#   SetNoiseBounds:  Add a 'noise' contribution with the given bounds.
#   Summary:  A data.frame summarizing predictions for multiple model subsets.
#   Train:  Optimize this Model's parameters to describe the given data.
#
# This class should work just fine, as long as
#   a) we are training on all the datapoints together (i.e., not breaking them
#      up into subregions to divide-and-conquer), and
#   b) this$params returns a vector which is amenable to simple optimization
#      routines (i.e., none of the parameters require special treatment).
# If either of these conditions fail, a specialized subclass should be created
# instead.  Note that *both* these conditions fail for *both* scenarios
# considered in our Journal of Applied Crystallography paper, despite the fact
# that I wrote this software to perform the analysis for that paper.  That's
# OK; having these classes still makes it very much easier to build the
# specialized functions I need.

setConstructorS3("Model",
  function(id="") {
    # Constructs a Model object with the given ID.
    #
    # Args:
    #   id:  A string which labels this model.
    #
    # Returns:
    #   A Model with the given ID.

    extend(Object(), "Model",
      .id             = id,
      .K.chol         = LazyMatrix(),
      .L              = LazyMatrix(),
      .last.trained.d = NA,
      .contributions  = list()
      )
  })

#-------------------------------------------------------------------------------
# (Model) PUBLIC "HELPER" FUNCTIONS:

LogML <- function(par=model$getParams(for.training=TRUE), model, d,
  update.params=TRUE) {
  # Sets the parameter values to 'par' for Model 'model', and returns the log
  # of the (M)arginal (L)ikelihood for describing Dataset 'd'.
  #
  # Args:
  #   par:  The parameter values to test.
  #   model:  The Model object we're optimizing.
  #   d:  The Dataset we're training on.
  #   update.params:  logical; if TRUE, we should change model's params to the
  #      values in par.
  #
  # Returns:
  #   The log of the marginal likelihood (also has a side-effect of setting
  #   model$params <- par, unless update.params is FALSE (although this would
  #   create an extra copy of model, and thus probably be less efficient).
  if (!update.params) {
    model <- clone(model)
  }
  model$setParams(p=DecodeForTraining(par))
  Y <- d$xformedDpts
  # The following calculation is based on Equation 5.8 in Rasmussen and
  # Williams, "Gaussian Processes for Machine Learning".
  term.data.fit   <- -0.5 * t(Y) %*% model$KInv(d) %*% Y
  term.complexity <- -0.5 * model$LogDetK(d)
  term.num.dpts   <- -0.5 * d$n * log(2 * pi)
  return (term.data.fit + term.complexity + term.num.dpts)
}

LogML1D <- function(value, name, model, d, update.params=TRUE) {
  # Calculate the log of the marginal likelihood for Model 'model' on Dataset
  # 'd', if model parameter given by 'name' is replaced by 'value'.
  #
  # Args:
  #   value:  The new value of the parameter
  #   name:  The name of the parameter
  #   model:  The Model object we're optimizing.
  #   d:  The Dataset we're training on.
  #   update.params:  logical; if TRUE, we should change model's params to the
  #      given value.
  #
  # Returns:
  #   The log of the marginal likelihood.
  names(value) <- name
  return (LogML(par=value, model=model, d=d, update.params=update.params))
}

GradLogML <- function(par=model$getParams(for.training=TRUE), model, d,
  update.params=TRUE) {
  # Sets the parameter values to 'par' for Model 'model', and returns the
  # gradient (w.r.t. par) of the log of the (M)arginal (L)ikelihood for
  # describing Dataset 'd'.
  #
  # Args:
  #   par:  The parameter values to test.
  #   model:  The Model object we're optimizing.
  #   d:  The Dataset we're training on.
  #   update.params:  logical; if TRUE, we should change model's params to the
  #      given value.
  #
  # Returns:
  #   The gradient of the log of the marginal likelihood (also has a
  #   side-effect of setting model$params <- par).
  if (!update.params) {
    model <- clone(model)
  }
  model$setParams(p=DecodeForTraining(par))
  Y <- d$xformedDpts
  # The following calculations are based on Equation 5.9 in Rasmussen and
  # Williams, "Gaussian Processes for Machine Learning".
  alpha <- model$KInv(d) %*% Y
  mat.for.grad <- alpha %*% t(alpha) - model$KInv(d)
  var.names <- names(model$getParams(for.training=TRUE))
  good.names <- names(par)[which(names(par) %in% var.names)]
  grad <- c()
  for (p.n in good.names) {
    grad[p.n] <- 0.5 * SmartTrace(model$KDeriv(d=d, param=p.n), mat.for.grad)
  }
  return (grad)
}

#-------------------------------------------------------------------------------
# (Model) PUBLIC VIRTUAL FIELDS:

# Virtual Field: contributionIds
# A list of id's for all contributions in this Model.
#
# Returns:
#   The id's for this Model's contributing Covariance objects.
setMethodS3("getContributionIds", "Model", conflict="quiet",
  function(this, ...) {
    id.list <- c()
    for (covar in this$.contributions) {
      id.list <- c(id.list, covar$id)
    }
    return (id.list)
  })

# Virtual Field: id
# A character string identifying this Model object.
#
# Args:
#   id: the string to change the id to.
#
# Returns:
#   The id of this Model object.
setMethodS3("getId", "Model", conflict="quiet",
  function(this, ...) {
    return (this$.id)
  })
setMethodS3("setId", "Model", conflict="quiet",
  function(this, id, ...) {
    this$.id <- id
    return (this)
  })

# Virtual Field: params
# A named vector of parameters governing this Model.
#
# Args:
#   p: A (named) vector of new parameter values (we ONLY use ones which are
#      named, and whose names match up with names of parameters.)
#   for.training:  If TRUE, we ignore "constant" parameters (i.e., where
#      lower=upper) and return the *log* of any "scale" parameters (such as ell
#      or sigma.f for the SE model).
#
# Returns:
#   A vector with values for ell and sigma.f
setMethodS3("getParams", "Model", conflict="quiet",
  function(this, for.training=FALSE, ...) {
    p <- c()
    for (covar in this$.contributions) {
      p <- c(p, covar$getParams(for.training=for.training))
    }
    if (for.training) {
      unlog.params <- DecodeForTraining(p)
      i.vary <- which(names(unlog.params) %in% this$getVaryingParamNames())
      p <- p[i.vary]
    }
    return (p)
  })
setMethodS3("setParams", "Model", conflict="quiet",
  function(this, p, for.training=FALSE, ...) {
    for (covar in this$.contributions) {
      covar$setParams(p=p, for.training=for.training)
    }
    return (invisible(this))
  })

# Virtual Field: lower
# Lower bounds for the parameter values.
#
# Args:
#   L: A (named) vector of new lower bounds (we ONLY use ones which are named,
#      and whose names match up with names of parameters.)
#   for.training:  If TRUE, we ignore "constant" parameters (i.e., where
#      lower=upper) and return the *log* of any "scale" parameters (such as ell
#      or sigma.f for the SE model).
#
# Returns:
#   The lower bounds for the parameters for this model.
setMethodS3("getLower", "Model", conflict="quiet",
  function(this, for.training=FALSE, ...) {
    L <- c()
    for (covar in this$.contributions) {
      L <- c(L, covar$getLower(for.training=for.training))
    }
    return (L)
  })
setMethodS3("setLower", "Model", conflict="quiet",
  function(this, L, for.training=FALSE, ...) {
    for (covar in this$.contributions) {
      covar$setLower(L=L, for.training=for.training)
    }
    return (invisible(this))
  })

# Virtual Field: upper
# Upper bounds for the parameter values.
#
# Args:
#   U: A (named) vector of new upper bounds (we ONLY use ones which are named,
#      and whose names match up with names of parameters.)
#   for.training:  If TRUE, we ignore "constant" parameters (i.e., where
#      lower=upper) and return the *log* of any "scale" parameters (such as ell
#      or sigma.f for the SE model).
#
# Returns:
#   The upper bounds for the parameters for this model.
setMethodS3("getUpper", "Model", conflict="quiet",
  function(this, for.training=FALSE, ...) {
    U <- c()
    for (covar in this$.contributions) {
      U <- c(U, covar$getUpper(for.training=for.training))
    }
    return (U)
  })
setMethodS3("setUpper", "Model", conflict="quiet",
  function(this, U, for.training=FALSE, ...) {
    for (covar in this$.contributions) {
      covar$setUpper(U=U, for.training=for.training)
    }
    return (invisible(this))
  })

# Virtual Field: signalIds (read-only)
# A list of id's for all non-noise contributions in this Model.
#
# Returns:
#   The id's for this Model's non-noise contributing Covariance objects.
setMethodS3("getSignalIds", "Model", conflict="quiet",
  function(this, ...) {
    # The ids of contributions which are considered to be "signal", i.e., not
    # noise.
    return (this$getContributionIds()[
        which(this$getContributionIds() != 'noise')])
  })

# Virtual Field: varyingParamNames
# Names of the parameters which are not constant
#
# Args:
#   U: A (named) vector of new upper bounds (we ONLY use ones which are named,
#      and whose names match up with names of parameters.)
#
# Returns:
#   The upper bounds for the parameters for this model.
setMethodS3("getVaryingParamNames", "Model", conflict="quiet",
  function(this, ...) {
    U <- this$getUpper()
    L <- this$getLower()
    p.names <- names(U)
    i.const <- which(U[p.names] == L[p.names])
    if (length(i.const) > 0) {
      return (p.names[-i.const])
    } else {
      return (p.names)
    }
  })


#-------------------------------------------------------------------------------
# (Model) PUBLIC METHODS:

setMethodS3("AddCovariance", "Model", conflict="quiet",
  function(this, covariance, on.duplicate.id="rename", ...) {
    # Add another Covariance structure to this model (it contributes
    # additively).
    #
    # Args:
    #   covariance:  A Covariance object to be cloned and added to this model.
    #   on.duplicate.id:  Character string, one of ("rename", "replace"),
    #      directing how to handle a new contribution whose ID is the same as
    #      an existing one's.
    #
    # Returns:
    #   Used for its side-effect.

    # We CLONE it so we OWN it.  (Don't want anyone else to fiddle with the
    # Covariance object, EXCEPT this Model object.)
    our.covar <- clone(covariance)

    # Make sure this contribution has a unique id:
    if (on.duplicate.id == "rename") {
      label <- 1
      unique.id <- our.covar$id
      old.ids <- this$contributionIds
      while (unique.id %in% old.ids) {
        unique.id <- paste(sep='', our.covar$id, '.', label)
        label <- label + 1
      }
      our.covar$id <- unique.id
    } else if (on.duplicate.id == "replace") {
      this$.contributions[[our.covar$id]] <- NULL
    }

    # We want to be able to refer to this contribution by its id.
    new.contribution <- list(our.covar)
    names(new.contribution) <- our.covar$id
    this$.contributions <- c(this$.contributions, new.contribution)
    return (invisible(this))
  })

clone.Model <- function(this, ...) {
  # Deep clone this Model object (clones Covariance's too).
  #
  # Returns:
  #   A deep clone of the model object.
  M <- clone.Object(this)
  M$.contributions <- list()
  for (covar in this$.contributions) {
    covar.clone <- list(clone(covar))
    names(covar.clone) <- covar$id
    M$.contributions <- c(M$.contributions, covar.clone)
  }
  M$.K.chol <- clone(M$.K.chol)
  M$.L <- clone(M$.L)
  return (M)
}

setMethodS3("Forget", "Model", conflict="quiet",
  function(this, ...) {
    this$.K.chol <- LazyMatrix()
    this$.L <- LazyMatrix()
    this$.last.trained.d <- NA
    gc()
  })

setMethodS3("Freeze", "Model", conflict="quiet",
  function(this, p.names=names(this$params), ...) {
    good.i <- which(p.names %in% names(this$params))
    for (p.name in p.names[good.i]) {
      this$lower <- this$params[p.name]
      this$upper <- this$params[p.name]
    }
  })

setMethodS3("L", "Model", conflict="quiet",
  function(this, d, X.out=d$X, contributions=this$getSignalIds(), ...) {
    # Compute the lower-triangular Cholesky decomposition of the covariance
    # matrix (useful for generating random draws).
    # 
    # Args:
    #   d:  The Dataset to evaluate the Model on.
    #   X.out:  matrix; the X-points where we want to predict the function.
    #   contributions:  character vector; id's of the contributing Covariances
    #      we want to predict (other contributions are considered noise);
    #      default is every contribution not named 'noise'.
    #
    # Returns:
    #   The lower-triangular Cholesky decomposition of the covariance matrix.
    this$ComputeL(d=d, X.out=X.out, contributions=contributions)
    return (this$.L$M)
  })

setMethodS3("NamedCovariance", "Model", conflict="quiet",
  function(this, id, ...) {
    # Retrieve a clone of the first contributing Covariance objects with the
    #   given id.
    #
    # Args:
    #   id:  The id (character) of the covariance to retrieve and clone.
    #
    # Returns:
    #   A clone of the Covariance object with the given id (or else NA if no
    #   such object exists).
    for (covar in this$.contributions) {
      if (covar$id == id) {
        Cov <- clone(covar)
        return (Cov)
      }
    }
    return (NA)
  })

setMethodS3("PlotBubblingSurfaces2D", "Model", conflict="quiet",
  function(this, d, X.out=d$X, contributions=this$getSignalIds(),
    n.surfaces=10, n.times=50, file.name=NULL, ...) {
    # Plot smoothly varying random surfaces to visualize the uncertainty.
    #
    # Args:
    #   d:  The Dataset to evaluate the Model on.
    #   X.out:  matrix; the X-points where we want to predict the function.
    #   contributions:  character vector; id's of the contributing Covariances
    #      we want to predict (other contributions are considered noise).
    #   n.surfaces:  The number of independent draws to take for each bubbler.
    #   n.times:  The final number of interpolated points.
    #   file.name:  Character, the basename of the file to plot to.
    #
    # Returns:
    #   Used for its side-effect (plotting a series of timesteps to PNG files).
    B <- BubblingRandomMatrix(n.pts=NumPoints(X.out), N=n.surfaces,
      n.times=n.times, ...)
    L <- this$L(d=d, X.out=X.out, contributions=contributions)
    Y <- matrix(nrow=NumPoints(X.out), rep(times=n.times,
        this$PosteriorMean(d=d, X.out=X.out, contributions=contributions))) + (
        L %*% B)
    unit <- min(dist(d$X))
    for (i in 1:ncol(Y)) {
      rgl.clear()
      PlotSurface(X=X.out, Y=Y[, i], ...)
      rgl.spheres(x=d$X[, 1], z=d$X[, 2], y=d$dpts, radius=0.2 * unit)
      rgl.snapshot(filename=sprintf("%s_t-%04d.png", file.name, i), top=TRUE)
    }
  })

setMethodS3("PlotCovariance2D", "Model", conflict="quiet",
  function(this, d, file.name=NULL, file.type='png', i=1, ...) {
    # Plot the covariance of points in 2D with respect to another point; also
    # plot the derivatives w.r.t. each parameter.
    #
    # Args:
    #   d:  The Dataset to evaluate the Model on.
    #   file.name:  Character, the name of the file to plot to (file extension
    #      should be automatically appended if necessary)
    #   file.type:  Character, one of ("pdf", "png")
    #   i:  numeric; index of datapoint to use for calculating the covariance
    #
    # Returns:
    #   Used for its side-effects of plotting.
    #
    # Notes:
    #   Lay them out in a grid of two rows.  First column: top is covariance
    #   matrix, bottom is for text (parameter values).  Subsequent columns show
    #   derivative matrices.

    # Set the layout:
    n.rows <- 2
    n.cols <- 1 + ceiling(length(this$params) / n.rows)
    Layout <- grid.layout(nrow = n.rows, ncol = n.cols, 
      widths  = unit(rep(1, n.cols), rep("null", n.cols)), 
      heights = unit(rep(1, n.rows), rep("null", n.rows)))

    # Open the file:
    my.file <- SetupFileInfo(name=file.name, type=file.type)
    if (my.file$type == 'png') {
      if (!require("Cairo")) {
        png(filename=my.file$name, ...)
      } else {
        CairoPNG(filename=my.file$name, ...)
      }
    } else if (my.file$type == 'pdf') {
      cairo_pdf(filename=my.file$name, ...)
    } else {
      stop("Sorry, we can only do pdf or png output for now.")
    }
    LayoutNewGridPage(Layout=Layout)

    # Construct a data.frame for all the matrices
    point.out <- matrix(d$X[i, ], nrow=1)
    K.data <<- data.frame(X.1=d$X[, 1], X.2=d$X[, 2],
      Cov=this$KTotal(d=d)[i, ],
      M=t(this$PredictionMatrix(d=d, X.out=point.out,
          contributions=this$contributionIds[
            which(this$contributionIds != 'noise')]))
      )

    # Generic ggplot setup
    ps.out <- 1.3  # outer point size
    ps.in <- 1.0   # inner point size
    p.base <- (ggplot(data=K.data, aes(x=X.1, y=X.2))
      + scale_colour_gradientn(colours=c('white','red','blue'))
      + geom_point(colour='black', size=ps.out)
      + geom_vline(xintercept=K.data[i, 1])
      + geom_hline(yintercept=K.data[i, 2])
      )

    # Print the covariance matrix in the top left
    p.cov <- (p.base
      + geom_point(aes(colour=Cov), size=ps.in)
      )
    print(p.cov, vp=Subplot(1, 1))
    # Weight matrix in bottom left
    p.weight <- (p.base
      + geom_point(aes(colour=M), size=ps.in)
      + scale_colour_gradient2(colours=c('red', 'white', 'blue'))
      )
    print(p.weight, vp=Subplot(2, 1))

    dev.off()
    return (invisible(this))
  })

setMethodS3("PosteriorMean", "Model", conflict="quiet",
  function(this, d, X.out=d$X, contributions=this$getSignalIds(),
    untransform.result=TRUE, ...) {
    # Computes this Model's optimal prediction of the underlying function's
    # value at every point in 'X.out'.
    #
    # Args:
    #   d:  The Dataset to train the Model on.
    #   X.out:  matrix; the X-points where we want to predict the function.
    #   contributions:  character vector; id's of the contributing Covariances
    #      we want to predict (other contributions are considered noise);
    #      default is every contribution not named 'noise'.
    #   untransform.result:  logical; if TRUE, we transform back to the space
    #      of dpts (as opposed to the space of xformedDpts where training takes
    #      place).
    #
    # Returns:
    #   A numeric vector with optimal predictions at every point in X.out.
    contributions <- this$CheckContributionsAndWarn(contributions)
    M <- this$PredictionMatrix(d=d, X.out=X.out, contributions=contributions)
    result <- M %*% d$xformedDpts
    if (untransform.result) {
      result <- d$Untransform(result)
    }
    return (result)
  })

setMethodS3("PredictionMatrix", "Model", conflict="quiet",
  function(this, d, X.out=d$X, contributions=this$getSignalIds(), ...) {
    # A matrix relating function values at the output points to function values
    # at the input points.
    #
    # Args:
    #   d:  The Dataset to train the Model on.
    #   X.out:  matrix; the X-points where we want to predict the function.
    #   contributions:  character vector; id's of the contributing Covariances
    #      we want to predict (other contributions are considered noise).
    #
    # Returns:
    #   The matrix 'M' such that (Y.pred = M %*% d$xformedDpts).
    K.in.out <- matrix(0, nrow=NumPoints(X.out), ncol=d$n)
    for (c.name in contributions) {
      covar <- this$.contributions[[c.name]]
      covar.K <- covar$KInOut(d=d, X.out=X.out)
      K.in.out <- K.in.out + covar.K
      rm(covar.K)
      gc()
    }
    M <- K.in.out %*% this$KInv(d=d)
    rm(K.in.out)
    gc()
    return (M)
  })

setMethodS3("PosteriorInterval", "Model", conflict="quiet",
  function(this, d, X.out=d$X, num.sd=1, contributions=this$getSignalIds(),
    ...) {
    # Computes bounds on the uncertainty (in terms of the standard deviation)
    # in the prediction at a given point, along with the prediction.
    #
    # Args:
    #   d:  The Dataset to train the Model on.
    #   X.out:  matrix; the X-points where we want to predict the function.
    #   num.sd:  The number of standard deviations from the mean our interval
    #      should include.
    #   contributions:  character vector; id's of the contributing Covariances
    #      we want to predict (other contributions are considered noise);
    #      default is every contribution not named 'noise'.
    #
    # Returns:
    #   A numeric vector with optimal predictions at every point in X.out.
    M <- this$PredictionMatrix(d=d, X.out=X.out, contributions=contributions)
    std.dev <- this$PosteriorStandardDeviation(d=d, X.out=X.out,
      contributions=contributions, ...)
    prediction <- M %*% d$xformedDpts
    return (data.frame(X=X.out,
        mean=d$Untransform(prediction),
        lower=d$Untransform(prediction - num.sd * std.dev),
        upper=d$Untransform(prediction + num.sd * std.dev)))
  })

setMethodS3("PosteriorStandardDeviation", "Model", conflict="quiet",
  function(this, d, X.out=d$X, contributions=this$getSignalIds(), ...) {
    # Computes the posterior "sigma" at a given point.
    #
    # Args:
    #   d:  The Dataset to train the Model on.
    #   X.out:  matrix; the X-points where we want to predict the function.
    #   contributions:  character vector; id's of the contributing Covariances
    #      we want to predict (other contributions are considered noise);
    #      default is every contribution not named 'noise'.
    #
    # Returns:
    #   A numeric vector with the posterior standard deviation at every point
    #   in X.out.
    contributions <- this$CheckContributionsAndWarn(contributions)
    # Calculate the posterior predictive mean.
    N.out <- NumPoints(X.out)
    K.in.out <- matrix(0, nrow=N.out, ncol=d$n)
    variance <- rep(0, N.out)
    for (c.name in contributions) {
      covar <- this$.contributions[[c.name]]
      K.in.out <- K.in.out + covar$KInOut(d=d, X.out=X.out)
      variance <- variance + covar$Variance(X=X.out)
      gc()
    }
    M <- this$PredictionMatrix(d=d, X.out=X.out, contributions=contributions)
    # The 'rowSums' bit is a fancy way to calculate only the diagonal elements
    # of the matrix product (K.in.out %*% this$KInv(d=d) %*% t(K.in.out)).
    # We force the variance to be non-negative; this is necessary because
    # roundoff errors in the matrix multiplication can cause it to go negative
    # by negligible amounts.
    variance <- pmax(0, variance - rowSums(M * K.in.out))
    std.dev <- sqrt(variance)
    return (std.dev)
  })

print.Model <- function(this, indent=0, ...) {
  # Pretty-prints information for this Model object.
  #
  # Returns:
  #   Used for its side-effect.
  tab <- Spaces(num=indent)
  cat(sprintf("%s%s, id='%s'\n", tab, class(this)[1], this$id))
  cat(sprintf("%s%sCONTRIBUTING COVARIANCES:\n", tab, Spaces(2)))
  for (covar in this$.contributions) {
    cat(sprintf("%s%sid=%-20s (%s)\n", tab, Spaces(4), Wrap(covar$id, "'"),
        class(covar)[1]))
  }
  PrintParams(lower=this$lower, upper=this$upper, params=this$params,
    indent=indent)
  return (invisible(this))
}

setMethodS3("SetNoiseBounds", "Model", conflict="quiet",
  function(this, sigma.vals, ...) {
    # Set the range of values for the noise level in this model.
    #
    # Args:
    #   sigma.vals:  A numeric vector, whose range() sets the range of values
    #   for the noise level.
    #
    # Returns:
    #   Used for its side-effect.
    this$AddCovariance(CovarianceNoise(id="noise", sigma.bounds=sigma.vals),
      on.duplicate.id="replace")
  })

setMethodS3("Train", "Model", conflict="quiet",
  function(this, d, force.retrain=FALSE, ...) {
    # Optimize this Model's parameters so they describe the given data.
    #
    # Args:
    #   d:  Dataset; the data which our parameters should describe.
    #
    # Returns:
    #   Used for its side-effect (i.e., changing the values of this Model's
    #   parameters to reflect the given Dataset).
    if (force.retrain || !d$Same(d=this$.last.trained.d,
        compare=c("X", "xformedDpts", "noiseVar"))) {
      old.params <- this$params
      lower <- this$getLower(for.training=TRUE)
      upper <- this$getUpper(for.training=TRUE)
      params <- this$getParams(for.training=TRUE)
      if (length(params) > 1) {
        this$.opt <- optim(method="L-BFGS-B", control=list(fnscale=-1),
          par=params, lower=lower, upper=upper,
          fn=LogML, gr=GradLogML,
          # Extra parameters needed by 'fn' and 'gr':
          model=this, d=d)
        if (this$.opt$convergence > 50) {  # L-BFGS-B warning or error
          warning(paste(sep='',
              "L-BFGS-B optimization ran into trouble; params NOT changed:\n    ",
              this$.opt$message))
          cat("By the way, the params were as follows:\n        ", this$.opt$par, "\n")
          this$params <- old.params
        } else {
          this$.last.trained.d <- clone(d)
        }
      } else if (length(params) == 1) {
        opt <- optimize(f=LogML1D, maximum=TRUE,
          lower=lower, upper=upper,
          model=this, d=d, name=names(params))
        best <- opt$maximum
        names(best) <- names(params)
        this$params <- best
      }
    }
    return (invisible(this))
  })

#-------------------------------------------------------------------------------
# (Model) PRIVATE METHODS:

setMethodS3("CheckContributionsAndWarn", "Model", private=TRUE, conflict="quiet",
  function(this, contributions, ...) {
    # Check the list of names in 'contributions' to see which we actually have,
    # returning a validated list of names, and warning if any were requested
    # which do not exist, or if no valid contributions remain.
    #
    # Args:
    #   contributions:  A character vector with the names of the contributions
    #      which we should check.
    #
    # Returns:
    #   A character vector of the unique contribution names which actually
    #   exist.

    # Check which requested contributions are present in the model.
    contributions <- unique(contributions)
    existing <- which(contributions %in% this$contributionIds)
    bad.names <- contributions[-existing]
    if (length(bad.names) > 0) {
      culprits <- paste(sep='', collapse=' ', '"', bad.names, '"')
      warning(paste("The following contribution IDs do not exist:\n",
          culprits, "\n"))
    }
    contributions <- contributions[existing]

    # Make sure we have at least one contribution!
    if (length(contributions) < 0) {
      warning("No nonzero contributions actually supplied!\n")
    }

    return (contributions)
  })

setMethodS3("ComputeL", "Model", private=TRUE, conflict="quiet",
  function(this, d, X.out=d$X, contributions, ...) {
    # Compute the Cholesky decomposition of the model's current covariance
    # matrix for datapoints 'd'.
    #
    # Args:
    #   d:  The Dataset for evaluating this model.
    #   X.out:  matrix; the X-points where we want to predict the function.
    #   contributions:  character vector; id's of the contributing Covariances
    #      we want to predict (other contributions are considered noise);
    #      default is every contribution not named 'noise'.
    #
    # Returns:
    #   The Cholesky decomposition of the output points' covariance matrix for
    #   this model and dataset.
    ingredients <- list(X=d$X, X.out=X.out, noiseVar=d$noiseVar,
      params=this$params)
    if (this$.L$NeedToRecalculate(ingredients=ingredients)) {
      this$ComputeKChol(d=d)
      K.out.out <- d$noiseVar * diag(NumPoints(X.out))
      K.in.out <- matrix(0, nrow=NumPoints(X.out), ncol=NumPoints(d$X))
      for (c.id in contributions) {
        covar <- this$.contributions[[c.id]]
        K.out.out <- K.out.out + covar$KOutOut(X.out=X.out)
        K.in.out <- K.in.out + covar$KInOut(d=d, X.out=X.out)
      }
      K.posterior <- K.out.out - (K.in.out %*% chol2inv(this$.K.chol$M) %*%
        t(K.in.out))
      this$.L$StoreMatrix(M=t(chol(K.posterior)), ingredients=ingredients)
    }
    return (invisible(this))
  })

setMethodS3("ComputeKChol", "Model", private=TRUE, conflict="quiet",
  function(this, d, ...) {
    # Compute the Cholesky decomposition of the model's current covariance
    # matrix for datapoints 'd'.
    #
    # Args:
    #   d:  The Dataset for evaluating this model.
    #
    # Returns:
    #   The Cholesky decomposition of the total covariance matrix for this
    #   model.
    ingredients <- list(X=d$X, noiseVar=d$noiseVar, params=this$params)
    if (this$.K.chol$NeedToRecalculate(ingredients=ingredients)) {
      K.tot <- this$KTotal(d=d)
      K.chol <- DebugIfError(chol.default, K.tot)
      this$.K.chol$StoreMatrix(M=K.chol, ingredients=ingredients)
    }
    return (invisible(this))
  })

setMethodS3("KDeriv", "Model", private=TRUE, conflict="quiet",
  function(this, d, param, ...) {
    # Computes the element-by-element derivative of the total K-matrix for
    # Dataset 'd', with respect to the parameter named 'param'.
    #
    # Args:
    #   d:  The Dataset for evaluating this model.
    #   param:  The name of the parameter with respect to which we're
    #      differentiating .
    #
    # Returns:
    #   The element-by-element derivative of the total K-matrix for Dataset
    #   'd', with respect to the parameter named 'param'.
    K.deriv <- matrix(0, nrow=d$n, ncol=d$n)
    in.logspace <- (length(grep(pattern=LogspaceTag(), x=param)) > 0)
    for (covar in this$.contributions) {
      if (in.logspace) {  # Decode the param name, if necessary
        param <- sub(pattern=LogspaceTag(), replacement="", x=param)
      }
      d.covar_d.param <- covar$KInInDeriv(d=d, param=param)
      if (in.logspace) {  # d/d(log(x)) = x*d/d(x)
        d.covar_d.param <- d.covar_d.param * this$params[param]
      }
      K.deriv <- K.deriv + d.covar_d.param
    }
    return (K.deriv)
  })

setMethodS3("K", "Model", private=TRUE, conflict="quiet",
  function(this, X, X.out=X, contributions, ...) {
    # Calculate a covariance matrix from X to X.out, including only the named
    # contributions.
    #
    # Args:
    #   X:  2-column numeric matrix; the input points (where we have data).
    #   X.out:  2-column numeric matrix; The output points (where we wish to
    #      make predictions).
    #   contributions:  character vector; a list of names of covariances to
    #      include.
    #
    # Returns:
    #   The covariance matrix from X to X.out.
    K <- matrix(0, nrow=NumPoints(X.out), ncol=NumPoints(X))
    # Add covariance matrix for each additive contribution:
    for (covar in this$.contributions) {
      if (covar$id %in% contributions) {
        K <- K + covar$K.specific(X=X, X.out=X.out)
      }
    }
    return (K)
  })

setMethodS3("KTotal", "Model", private=TRUE, conflict="quiet",
  function(this, d, ...) {
    # Compute the total covariance matrix for this model with respect to the
    # data in Dataset 'd'.
    #
    # Args:
    #   d:  The Dataset for evaluating this model.
    #
    # Returns:
    #   The total covariance matrix for this model.
    K.tot <- matrix(0, nrow=d$n, ncol=d$n)
    for (covar in this$.contributions) {
      K.tot <- K.tot + covar$KInIn(d=d)
      gc()
    }
    # Cap it off with the noise associated with the data:
    diag(K.tot) <- diag(K.tot) + d$noiseVar
    return (K.tot)
  })

setMethodS3("KInv", "Model", private=TRUE, conflict="quiet",
  function(this, d, ...) {
    # The inverse total covariance matrix for this model, with respect to the
    # given dataset.
    #
    # Args:
    #   d:  The Dataset for evaluating this model.
    #
    # Returns:
    #   The inverse total covariance matrix for this model.
    this$ComputeKChol(d=d)
    return (chol2inv(this$.K.chol$M))
  })

setMethodS3("LogDetK", "Model", private=TRUE, conflict="quiet",
  function(this, d, ...) {
    # Compute the logarithm of the determinant of the model's total covariance
    # matrix for the points in Dataset 'd'.
    #
    # Args:
    #   d:  Dataset at whose points we evaluate the covariance matrix.
    #
    # Returns:
    #   The logarithm of the determinant of the model's covariance matrix.
    this$ComputeKChol(d=d)
    return (2 * sum(log(diag(this$.K.chol$M))))
  })

################################################################################
# SECTION:                    Standalone Functions
################################################################################

