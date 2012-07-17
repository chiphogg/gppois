\name{Dataset}
\alias{Dataset}
\title{A wrapper class for data being analyzed}
\arguments{
  \item{id}{(character) An id which identifies this
  Dataset.}

  \item{data}{(\code{data.frame}) The raw data.}

  \item{X.names}{(character vector) The names of the
  columns which represent covariates.  (NB: the length of
  this vector implicitly defines the dimensionality of the
  data: if you give 2 column names, you have a 2D
  dataset.)}

  \item{noise.var}{(named numeric vector) Desired noise
  levels for each column (i.e., a threshhold for precision:
  we don't care about differences smaller than the
  corresponding scale).}

  \item{column}{(character) The name of the quantity to
  select for analysis.}

  \item{poisson.names}{(character vector) The column names
  which represent Poisson data.}

  \item{tol.factor}{(numeric) The relative tolerance
  factor: by default, we set the noise variance for each
  column to \code{tol.factor * sd(col)}.}

  \item{data.offset}{(named numeric vector) An optional
  offset which is subtracted from each column.  If an
  unnamed value is included, it is treated as a default.
  The "default default" is to subtract off the mean for
  each column.}

  \item{...}{Not used.}
}
\description{
  This is the constructor for a \code{Dataset} object:
  basically, a wrapper for data.frame objects with some
  added knowledge that helps out the Model.  For example:
  \itemize{ \item It knows which columns are independent
  (\sQuote{x}-type) variables, and which columns represent
  datapoints \item It knows the dimensionality of the data,
  and how many datapoints there are \item It keeps track of
  both the raw datapoints, and \dQuote{transformed}
  versions which are easier to analyze with Gaussian
  processes.  \item It knows }
}
\section{Data transformations}{
  There are two "types" of data transformation: baseline
  adjustments (for Gaussian data), and the Anscombe
  transform (for turning Poisson into Gaussian).

  \describe{ \item{\strong{Baseline adjustments}}{
  Currently, this amounts to subtracting a constant from
  the data (defaults to the mean).  However, future
  versions should be more flexible: let the user specify
  anything up to quadratic using \dQuote{const},
  \dQuote{lin}, or \dQuote{quad}, and keep track of the
  coefficients (so we can apply it to regions outside or
  between the datapoints).} \item{\strong{Anscombe
  transform}}{ See \code{\link{Anscombe}}.  This
  transformation takes data with independent Poisson noise,
  and turns it into independent \emph{Gaussian} noise --
  with a \strong{constant} variance of 1/4.  Basically,
  this ability puts the \dQuote{pois} in \dQuote{gppois}.
  } }
}

