\name{FindGapPoints}
\alias{FindGapPoints}
\title{Flag points which are "in the gap"}
\usage{
  FindGapPoints(X, X.out, gap.thresh = NA)
}
\arguments{
  \item{X}{d-column matrix whose rows represent points
  where we have data.}

  \item{X.out}{d-column matrix whose rows represent points
  where we want to predict.}

  \item{gap.thresh}{The threshhold distance which
  determines \sQuote{gap} points, i.e., \sQuote{non-gap}
  points are closer than \code{gap.thresh} to any
  datapoint.}
}
\value{
  A logical matrix, of length NumPoints(X.out), TRUE if the
  corresponding point in X.out is a "gap point".
}
\description{
  This function takes two sets of points: \code{X} is
  usually points where we have data, and \code{X.out} is
  usually gridded points where we evaluate the Gaussian
  Process.  It looks for points in \code{X.out} which are
  "gap" points, i.e., far enough away from every point in
  \code{X}.  The threshhold is \code{gap.thresh}, which
  defaults to the largest nearest-neighbor distance in X.
}

