\name{Dataset$quantity}
\alias{Dataset$quantity}
\alias{getQuantity}
\alias{getQuantity.Dataset}
\alias{quantity}
\alias{setQuantity}
\alias{setQuantity.Dataset}
\title{Currently selected quantity (column)}
\arguments{
  \item{...}{Not used.}
}
\value{
  (character) The name of the quantity which is currently
  selected.
}
\description{
  A Dataset consists of \code{d} covariates and one or more
  observed variables. Only one of these variables can be
  analyzed at a time; \code{quantity} gives the name of the
  variable currently selected.
}
\details{
  Usually accesed using \code{this$quantity} syntax.
}
\seealso{
  \code{\link{Dataset}}
}

