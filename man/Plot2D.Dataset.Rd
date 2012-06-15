\name{Plot2D.Dataset}
\alias{Plot2D.Dataset}
\title{Scatterplot for 2D dataset}
\arguments{
  \item{max.points}{The maximum number of points to show
  (these are randomly sampled from the available points).}

  \item{...}{Not used.}
}
\value{
  Used for its side-effect.
}
\description{
  Use rgl to plot 2D data.  If there are too many
  datapoints, plotting and interacting become very slow.
  Thus, this function allows to limit the number plotted
  using the \code{max.points} parameter.
}
\seealso{
  \code{\link{Dataset}}
}

