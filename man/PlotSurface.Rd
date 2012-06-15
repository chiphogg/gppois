\name{PlotSurface}
\alias{PlotSurface}
\title{Plot a triangulated surface using rgl}
\usage{
  PlotSurface(X, Y, dY = NA, tri = NA, new.window = FALSE,
    Y.scale = max(dist(X))/diff(range(Y)), ...)
}
\arguments{
  \item{X}{2-column matrix holding the X-points for data.}

  \item{Y}{numeric vector of datapoints}

  \item{dY}{optional numeric vector of uncertainties in
  datapoints}

  \item{tri}{3-column matrix of indices into X, whose rows
  are triangles.}

  \item{new.window}{Whether to open a new window}

  \item{Y.scale}{The factor for scaling Y when plotting.}
}
\description{
  This function takes a set of function values (\code{Y})
  evaluated at corresponding points in the 2D plane
  (\code{X}) and plots the corresponding surface.  It uses
  a triangulation of the points, which can be supplied (as
  \code{tri}) or can also be calculated automatically.  If
  dY is given, it also plots the +- 1 sigma surfaces
  translucently.  If \code{Y.scale} is not supplied, this
  function will guess the scale.
}

