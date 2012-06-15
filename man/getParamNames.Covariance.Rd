\name{getParamNames.Covariance}
\alias{Covariance$paramNames}
\alias{getParamNames.Covariance}
\title{ID-decorated names of Covariance parameters}
\arguments{
  \item{...}{Not used.}
}
\value{
  Parameter names in the form "id.basename", where basename
  is the "undecorated" parameter name.
}
\description{
  This function returns the names of the Covariance
  parameters, decorated with the Covariance ID as a prefix.
  This decoration helps distinguish the parameters if a
  model has multiple copies of a given type of covariance.
  For example, if we used both coarse and fine
  lengthscales, we might have parameters named "coarse.ell"
  and "fine.ell" instead of two parameters ambiguously
  named "ell".
}
\seealso{
  \code{\link{getParams.Covariance}}

  \code{\link{getParams.Model}}

  \code{\link{Covariance}}
}

