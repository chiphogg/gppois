\name{AddCovariance.Model}
\alias{AddCovariance}
\alias{AddCovariance.Model}
\title{Add a new Covariance to this Model}
\arguments{
  \item{covariance}{A Covariance object to be cloned and
  added to this model.}

  \item{on.duplicate.id}{(character) One of
  (\dQuote{rename}, \dQuote{replace}), directing how to
  handle a new contribution whose ID is the same as an
  existing one's.}

  \item{...}{Not used.}
}
\description{
  Add another Covariance structure to this model (it
  contributes additively). We CLONE it so we OWN it.
  (Don't want anyone else to fiddle with the Covariance
  object, EXCEPT this Model object.)
}
\seealso{
  \code{\link{Model}}

  \code{\link{Covariance}}
}

