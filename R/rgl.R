#' Plot a triangulated surface using rgl
#'
#' This function takes a set of function values (\code{Y}) evaluated at
#' corresponding points in the 2D plane (\code{X}) and plots the corresponding
#' surface.  It uses a triangulation of the points, which can be supplied (as
#' \code{tri}) or can also be calculated automatically.  If dY is given, it also
#' plots the +- 1 sigma surfaces translucently.  If \code{Y.scale} is not
#' supplied, this function will guess the scale.
#'
#' @param X  2-column matrix holding the X-points for data.
#' @param Y  numeric vector of datapoints
#' @param dY  optional numeric vector of uncertainties in datapoints
#' @param tri  3-column matrix of indices into X, whose rows are triangles.
#' @param new.window  Whether to open a new window
#' @param Y.scale  The factor for scaling Y when plotting.
#'
#' @export
PlotSurface <- function(X, Y, dY=NA, tri=NA, new.window=FALSE,
  Y.scale=max(dist(X)) / diff(range(Y)), ...) {
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


