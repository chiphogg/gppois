# R.oo methods defined using 'setMethodS3' automatically declare their
# generics.  This is usually nice -- but sometimes not nice.  For example:
#   * "pure virtual" methods (defined only in subclasses)
#     - e.g., the Variance function in Covariance subclasses
#   * related methods in different classes with the same (generic) name
#     - e.g., getParams in Model and Covariance
#
# The former should be documented in the superclass file.  The latter are
# documented herein.

#' Lower bounds for parameters
#'
#' Objects governed by parameters also have boundaries on those parameters.
#' These R.oo functions retrieve and/or set these lower bounds.
#'
#' Boundaries push obstacles aside when they are changed.  For example, suppose
#' a parameter and its bounds are set as (lower=3, param=4, upper=5).  Setting
#' lower to 6 will result in values (lower=6, param=6, upper=6): the other
#' values are "dragged along".
#'
#' @name lower
#' @aliases getLower setLower
#' @export getLower setLower
#'
#' @seealso \code{\link{upper}}
#' @seealso \code{\link{params}}
#' @seealso \code{\link{Covariance$upper}}
#' @seealso \code{\link{Model$upper}}
NULL

#' Upper bounds for parameters
#'
#' Objects governed by parameters also have boundaries on those parameters.
#' These R.oo functions retrieve and/or set these upper bounds.
#'
#' Boundaries push obstacles aside when they are changed.  For example, suppose
#' a parameter and its bounds are set as (lower=3, param=4, upper=5).  Setting
#' lower to 6 will result in values (lower=6, param=6, upper=6): the other
#' values are "dragged along".
#'
#' @name upper
#' @aliases getUpper setUpper
#' @export getUpper setUpper
#'
#' @seealso \code{\link{lower}}
#' @seealso \code{\link{params}}
#' @seealso \code{\link{Covariance$upper}}
#' @seealso \code{\link{Model$upper}}
NULL

#' Parameter values
#'
#' These R.oo functions retrieve and/or set the values of parameters which
#' govern an object (i.e., typically a \code{\link{Covariance}}, or a
#' \code{\link{Model}} which is a collection of Covariances).
#'
#' When parameters are moved, they cannot push the boundaries aside.  For
#' example, suppose a parameter and its bounds are set as 
#' (lower=3, param=4, upper=5).  Setting param to 6 will result in values
#' (lower=3, param=5, upper=5): it is clamped to lie within the bounds.
#'
#' @name params
#' @aliases getParams setParams
#' @export getParams setParams
#'
#' @seealso \code{\link{lower}}
#' @seealso \code{\link{upper}}
#' @seealso \code{\link{Covariance$params}}
#' @seealso \code{\link{Model$params}}
NULL

