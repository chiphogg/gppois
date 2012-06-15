#############################################################################/**
# @RdocClass LazyMatrix
#
# Wrapper to avoid recomputing matrices
#
# \description{
#   This is the \emph{constructor} for a \code{LazyMatrix} object:
#   @get "title".
#
#   Note that this class doesn't do any actual computation!  It just stores the
#   results, and also tells other code whether or not it needs to.
#
#   Here is the class hierarchy: 
#   @classhierarchy
#
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods
# }
#
# @author
#*/###########################################################################
setConstructorS3("LazyMatrix",
  function(...) {
    extend(Object(), "LazyMatrix",
      .last.ingredients = list(),
      .matrix           = matrix()
      )
  })


#' The previously-computed matrix
#'
#' Accesses the matrix which was previously computed.
#'
#' @name getM.LazyMatrix
#' @aliases M getM
#' @S3method getM LazyMatrix
#'
#' @param ... Not used.
#'
#' @usage Dataset$M
#'
#' @return The previously-computed matrix.
#'
#' @seealso \code{\link{LazyMatrix}}
setMethodS3("getM", "LazyMatrix", conflict="quiet",
  function(this, ...) {
    return (this$.matrix)
  })

#' Decides whether we need to recompute
#'
#' Decides whether or not we would need to recompute the matrix, by checking
#' whether the supplied "ingredients" are equal to the ones previously supplied.
#'
#' @S3method NeedToRecalculate LazyMatrix
#' @name NeedToRecalculate.LazyMatrix
#'
#' @param ... Not used.
#'
#' @usage LazyMatrix$NeedToRecalculate
#'
#' @return TRUE iff we will need to recompute this matrix because the
#'    ingredients are different.
#'
#' @seealso \code{\link{LazyMatrix}}
setMethodS3("NeedToRecalculate", "LazyMatrix", conflict="quiet",
  function(this, ingredients, ...) {
    # Args:
    #   ingredients:  A named list of quantities this matrix depends on.
    #
    # Returns:

    # First, make sure we have all the old ingredients, and no new ones
    same.names <- (
      all(names(ingredients) %in% names(this$.last.ingredients)) &&
      all(names(this$.last.ingredients) %in% names(ingredients))
      )
    if (!same.names) {
      return (TRUE)
    }

    # Make sure the contents of corresponding ingredients are identical.
    for (name in names(ingredients)) {
      if (!identical(ingredients[[name]], this$.last.ingredients[[name]])) {
        return (TRUE)
      }
    }

    return (FALSE)
  })

#' Stores a calculated matrix
#'
#' Stores the supplied matrix in memory, and remembers the ingredients used to
#' calculate it.
#'
#' @S3method StoreMatrix LazyMatrix
#' @name StoreMatrix.LazyMatrix
#'
#' @param M The matrix to store.
#' @param ... Not used.
#'
#' @seealso \code{\link{LazyMatrix}}
setMethodS3("StoreMatrix", "LazyMatrix", conflict="quiet",
  function(this, M, ingredients, ...) {
    # Stores the supplied matrix in memory, along with the ingredients used to
    # calculate it.
    #
    # Args:
    #   M:  The matrix to store.
    #   ingredients:  A list of quantities used to compute 'M'.
    #
    # Returns:
    #   Used for its side-effect (i.e., storing these quantities).
    this$.matrix <- M
    this$.last.ingredients <- ingredients
    return (invisible(this))
  })

