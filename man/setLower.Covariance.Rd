\name{setLower.Covariance}
\alias{Covariance$lower}
\alias{getLower.Covariance}
\alias{setLower.Covariance}
\title{Lower bounds for parameters using full names}
\arguments{
  \item{L}{A (named) vector of new lower bounds on
  parameter values (we \emph{only} use ones which are
  named, and whose names match up with names of
  parameters.)}

  \item{for.training}{(logical) If TRUE, we return values
  more suitable for optimization: "scale"-type parameters
  are given as logarithms (which amounts to using the
  Jeffreys prior), and "constant" params (i.e., where the
  lower bound equals the upper bound) are omitted
  entirely.}

  \item{...}{Not used.}
}
\description{
  Gets or sets the lower bounds for one or more parameters,
  using the full name (ID + basename) of each parameter.
  See section on Undecorated Names in help page for
  \code{\link{paramsPlain}} for further commentary on name
  decoration.
}
\seealso{
  \code{\link{getParamsPlain}}

  \code{\link{Covariance}}
}

