\name{getParams.Covariance}
\alias{Covariance$params}
\alias{getParams.Covariance}
\alias{setParams.Covariance}
\title{Parameter values for this Covariance}
\arguments{
  \item{p}{(named numeric vector) New parameter values.
  (N.B: we \emph{only} use ones which are named, and whose
  names match up with names of parameters.)}

  \item{for.training}{(logical) If TRUE, we return values
  more suitable for optimization: "scale"-type parameters
  are given as logarithms (which amounts to using the
  Jeffreys prior), and "constant" params (i.e., where the
  lower bound equals the upper bound) are omitted
  entirely.}

  \item{...}{Not used.}
}
\value{
  (named numeric vector) Parameter values for this
  Covariance.
}
\description{
  The current values of the parameters which govern this
  Covariance.
}
\section{Optimization mode}{
  By default, we return all parameters.  However, this
  would if the caller is an \emph{optimization} routine,
  there are at least two important drawbacks.  First, if
  any parameters are fixed, it's wasteful (and potentially
  hazardous) to pass these to the optimization routine.
  Second, and more seriously, it assumes a flat prior on
  \emph{all} parameters, even "scale"-type parameters.
  This causes the optimization to fail outright even in
  many simple cases.

  To circumvent these problems, we provide an
  \emph{optimization mode}, which returns only non-constant
  parameters, and puts "scale"-type parameters in logspace.
  (The latter corresponds to the Jeffreys prior: flat in
  \emph{log-space}, rather than real-space.  It represents
  uncertainty about the \emph{order of magnitude} of the
  parameter.)
}

\section{Handling crossed boundaries}{
  The lower bound, upper bound, and value of every
  parameter must \emph{always} be properly ordered:\cr
  lower <= param <= upper\cr Sometimes a proposed move
  might violate that condition.  The way it's handled
  depends on whether it's a \emph{parameter} that was
  moved, or a \emph{boundary}: \itemize{ \item Moved
  \emph{parameter}\cr The boundaries are hard.  The value
  will be clamped at the closest allowable value, and a
  warning will be given.  (Note that this implies a
  parameter can be set "constant" by setting its upper
  bound equal to its lower bound.)  \item Moved
  \emph{boundary}\cr A moved boundary can "push" the
  parameter value, or even the other boundary.  e.g., if
  you have (lower=3, param=4, upper=5) and set upper to 2,
  the final state will be (lower=2, param=2, upper=2).  }
}
\seealso{
  \code{\link{getParamsPlain.Covariance}}

  \code{\link{Covariance}}
}

