\name{SetNoiseBounds.Model}
\alias{SetNoiseBounds.Model}
\title{Uncertainty about the noise level}
\arguments{
  \item{this}{The Model object.}

  \item{sigma.vals}{A numeric vector, such that
  range(sigma.vals) sets the range of values for the noise
  level.}

  \item{...}{Not used.}
}
\description{
  The noise level is not necessarily known precisely a
  priori.  If it is uncertain, this function lets you set a
  range of possible values.
}
\details{
  Keep in mind that the noise is a \dQuote{scale-type}
  parameter, so we are uncertain about its \emph{order of
  magnitude}.  Thus, it's a good idea to pass boundaries
  like, say, \code{c(1e-7, 1e-3)}, rather than boundaries
  which have the same order of magnitude.  (Unless of
  course you do know the order of magnitude ahead of time!)
}
\seealso{
  \code{\link{Model}}
}

