#' Anscombe transform for Poisson-noised data
#'
#' The Anscombe transform takes i.i.d. Poisson-noised data into i.i.d. normal
#' with constant variance (to a very good approximation).  Gaussian-noised data
#' is much easier to fit with Gaussian Processes.  The inverse Anscombe
#' transform takes functions from this transformed space back into the original
#' space.
#'
#' @rdname Anscombe
#' @aliases Anscombe AnscombeInverse
#'
#' @param Y.p  Data (presumably Poisson-noised) to be transformed
#' @param Y.g  Data (presumably Anscombe-transformed) to be untransformed to
#'    the original Poisson space.
#' @export
#' @return The transformed (\code{Anscombe}) or untransformed
#'    (\code{AnscombeInverse}) data
Anscombe <- function(Y.p) {
  return (sqrt(Y.p + 3.0 / 8.0))
}

#' @rdname Anscombe
AnscombeInverse <- function(Y.g) {
  return (Y.g ^ 2 - (1 / 8))
}

