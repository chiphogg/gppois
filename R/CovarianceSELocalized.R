############################################################################/**
# @RdocClass CovarianceSELocalized
#
# @title "Localized Squared-exponential Covariance"
#
# \description{
#   A squared-exponential covariance which is limited in spatial extent.
#   In addition to the two usual parameters (i.e., the horizontal and vertical
#   lengthscales), there are two \dQuote{boundary} parameters, \code{X.L} and
#   \code{X.R}.  \code{sigma.f} transitions smoothly to zero outside the region
#   between \code{X.L} and \code{X.R}.  The transition lengthscale is
#   \code{ell} (if it were any smaller, it could introduce sub-\code{ell}
#   features).
#
#   @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#   \item{id}{(character) A string to identify this covariance object.}
#   \item{ell}{(numeric) A characteristic horizontal scale for features in
#      functions being modeled.}
#   \item{sigma.f}{(numeric) A characteristic vertical scale for features in
#      functions being modeled.}
#   \item{X.L}{(numeric) The left boundary of the localized region.}
#   \item{X.R}{(numeric) The right boundary of the localized region.}
#   \item{ell.bounds}{(numeric) The range of values which \code{ell} might
#      assume.}
#   \item{sigma.f.bounds}{(numeric) The range of values which \code{sigma.f}
#      might assume.}
#   \item{X.L.bounds}{(numeric) The range of values which \code{X.L} might
#      assume.}
#   \item{X.R.bounds}{(numeric) The range of values which \code{X.R} might
#      assume.}
#   \item{...}{Not used.}
# }
#
# \section{Covariance Parameters}{
#   This section lists the fit parameters corresponding to this type of
#   Covariance.  Any parameters marked as \dQuote{(Scale parameter)} will be
#   optimized in log-space, consistent with the Jeffreys prior.
#
#   \describe{
#     \item{ell}{(Scale parameter) The horizontal feature lengthscale.}
#     \item{sigma.f}{(Scale parameter) The vertical feature lengthscale.}
#     \item{X.L}{The left boundary of the localized region.}
#     \item{X.R}{The right boundary of the localized region.}
#   }
# }
#
# \section{Fields and Methods}{
#  @allmethods
#
# }
#
# @author
#*/###########################################################################
setConstructorS3("CovarianceSELocalized", function(..., id="SE",
    ell=NA, sigma.f=NA, X.L=NA, X.R=NA, ell.bounds=NA, sigma.f.bounds=NA,
    X.L.bounds=NA, X.R.bounds=NA) {
    ell.good <- InitializeBoundedQuantity(ok.range=c(0, Inf),
      quantity=ell, bounds=ell.bounds, logspace=TRUE)
    sigma.f.good <- InitializeBoundedQuantity(ok.range=c(0, Inf),
      quantity=sigma.f, bounds=sigma.f.bounds, logspace=TRUE)
    X.L.good <- InitializeBoundedQuantity(quantity=X.L, bounds=X.L.bounds)
    X.R.good <- InitializeBoundedQuantity(quantity=X.R, bounds=X.R.bounds)
    # Construct the CovarianceSELocalized object:
    extend(Covariance(..., id=id), "CovarianceSELocalized",
      .ell            = ell.good$quantity,
      .ell.bounds     = ell.good$bounds,
      .sigma.f        = sigma.f.good$quantity,
      .sigma.f.bounds = sigma.f.good$bounds,
      .X.L            = X.L.good$quantity,
      .X.L.bounds     = X.L.good$bounds,
      .X.R            = X.R.good$quantity,
      .X.R.bounds     = X.R.good$bounds)
})

standard.mask <- function(x) {
  # A sigmoid mask which goes from (-Inf, -1) to (+Inf, +1), with a slope of 1
  # in the middle.
  #
  # Args:
  #   x:  The coordinate to evaluate the mask at.
  #
  # Returns:
  #   The standard sigmoid evaluated at x.
  return (as.vector(erf(sqrt(pi) * x * 0.5)))
}

standard.mask.deriv <- function(x) {
  # Derivative of the standard sigmoid mask 'standard.mask'.
  #
  # Args:
  #   x:  The coordinate to evaluate the mask derivative at.
  #
  # Returns:
  #   The derivative of the standard sigmoid, evaluated at x.
  return (as.vector(exp(-0.25 * pi * (x ^ 2))))
}

localize.mask <- function(X, X.L, X.R, ell) {
  # Mask values in X according to whether they're between X.L and X.R, with
  # transition lengthscale ell.
  #
  # Args:
  #   X:  The values to mask.
  #   X.L:  The left boundary of the region.
  #   X.L:  The right boundary of the region.
  #   ell:  The transition lengthscale at the edges of the region.
  upper.mask <- standard.mask((X - X.L) / ell)
  lower.mask <- standard.mask((X - X.R) / ell)
  return (as.vector(0.5 * (upper.mask - lower.mask)))
}

#' Names of "scale"-type parameters
#'
#' A character vector of names, indicating which parameters are considered to be
#' "scale" parameters.  (Read \dQuote{Optimization mode} section of
#' \code{\link{getParams.Covariance}} to see what this means.)
#'
#' @name getLogspaceNames.CovarianceSELocalized
#' @aliases CovarianceSELocalized$logspaceNames getLogspaceNames.CovarianceSELocalized
#' @S3method getLogspaceNames CovarianceSELocalized
#' @export getLogspaceNames getLogspaceNames.CovarianceSELocalized
#'
#' @param ... Not used.
#'
#' @usage CovarianceSELocalized$logspaceNames
#'
#' @return Names of parameters to be optimized in logspace.
#'
#' @seealso \code{\link{CovarianceSELocalized}}
setMethodS3("getLogspaceNames", "CovarianceSELocalized", conflict="quiet",
  function(this, ...) {
    return (c("ell", "sigma.f"))
  })

#' Basenames of parameters
#'
#' Gives the "basenames" (i.e. names undecorated by the id string) of the
#' parameters.
#'
#' @name getParamNamesPlain.CovarianceSELocalized
#' @aliases CovarianceSELocalized$paramNamesPlain getParamNamesPlain.CovarianceSELocalized
#' @S3method getParamNamesPlain CovarianceSELocalized
#' @export getParamNamesPlain getParamNamesPlain.CovarianceSELocalized
#'
#' @param ... Not used.
#'
#' @usage CovarianceSELocalized$paramNamesPlain
#'
#' @return The basenames of the parameters.
#'
#' @seealso \code{\link{CovarianceSELocalized}}
setMethodS3("getParamNamesPlain", "CovarianceSELocalized", conflict="quiet",
  function(this, ...) {
    return (c("ell", "sigma.f", "X.L", "X.R"))
  })

#' Parameter values with plain names
#'
#' Gives a vector of parameter values, whose names are NOT decorated by the id
#' of this Covariance object.
#'
#' @name getParamsPlain.CovarianceSELocalized
#' @aliases CovarianceSELocalized$paramsPlain getParamsPlain.CovarianceSELocalized
#' @S3method getParamsPlain CovarianceSELocalized
#' @export getParamsPlain getParamsPlain.CovarianceSELocalized
#'
#' @param ... Not used.
#'
#' @usage CovarianceSELocalized$paramsPlain
#'
#' @return The parameters for this covariance function, but with names
#'    undecorated by its id.
#'
#' @seealso \code{\link{setParamsPlain.Covariance}}
#' @seealso \code{\link{CovarianceSELocalized}}
setMethodS3("getParamsPlain", "CovarianceSELocalized", conflict="quiet",
  function(this, ...) {
    p <- c(this$.ell, this$.sigma.f, this$.X.L, this$.X.R)
    names(p) <- getParamNamesPlain(this)
    return (p)
})
setMethodS3("paramsPlainImplementation", "CovarianceSELocalized", conflict="quiet",
  private=TRUE,
  function(this, p, ...) {
    # Sets any values in 'p' which match our parameter names, without worrying
    # about lower and upper bounds (the function which calls this has the job
    # of worrying about these!).
    p.old <- this$getParamsPlain()
    to.change <- names(p)[which(names(p) %in% names(p.old))]
    p.old[to.change] <- p[to.change]
    this$.ell <- p.old["ell"]
    this$.sigma.f <- p.old["sigma.f"]
    this$.X.L <- p.old["X.L"]
    this$.X.R <- p.old["X.R"]
    return (invisible(this))
  })

#' Lower bounds for params, with plain names
#'
#' Gives a vector of lower bounds for the parameter values, whose names are NOT
#' decorated by the id of this Covariance object.
#'
#' @name getLowerPlain.CovarianceSELocalized
#' @aliases CovarianceSELocalized$lowerPlain 
#' @aliases getLowerPlain.CovarianceSELocalized
#' @aliases setLowerPlain.CovarianceSELocalized
#' @S3method getLowerPlain CovarianceSELocalized
#' @export getLowerPlain getLowerPlain.CovarianceSELocalized
#'
#' @param L A (named) vector of new lower bounds (we ONLY use ones which are
#'    named, and whose names match up with names of parameters.)
#' @param ... Not used.
#'
#' @usage CovarianceSELocalized$lowerPlain
#'
#' @return The lower bounds for the parameters for this covariance function, but
#'    with names undecorated by its id.
#'
#' @seealso \code{\link{getUpperPlain.CovarianceSELocalized}}
#' @seealso \code{\link{CovarianceSELocalized}}
setMethodS3("getLowerPlain", "CovarianceSELocalized", conflict="quiet",
  function(this, ...) {
    L <- c(this$.ell.bounds[1], this$.sigma.f.bounds[1], this$.X.L.bounds[1],
      this$.X.R.bounds[1])
    names(L) <- getParamNamesPlain(this)
    return (L)
})
setMethodS3("setLowerPlain", "CovarianceSELocalized", conflict="quiet",
  function(this, L, ...) {
    if (length(L) < 1) {
      return(invisible(this))
    }
    posdef.names <- c('ell', 'sigma.f')
    L[posdef.names] <- pmax(L[posdef.names], 0)
    # Adjust upper bounds to make way for the new values of L
    L.change <- this$PushUpperBounds(U.min=L)

    L.vals <- this$getLowerPlain()
    L.vals[names(L.change)] <- L.change[names(L.change)]
    this$.ell.bounds[1] <- L.vals["ell"]
    this$.sigma.f.bounds[1] <- L.vals["sigma.f"]
    this$.X.L.bounds[1] <- L.vals["X.L"]
    this$.X.R.bounds[1] <- L.vals["X.R"]
    this$ClampParams(warn=TRUE)
    return (this)
  })

#' Upper bounds for params, with plain names
#'
#' Gives a vector of upper bounds for the parameter values, whose names are NOT
#' decorated by the id of this Covariance object.
#'
#' @name getUpperPlain.CovarianceSELocalized
#' @aliases CovarianceSELocalized$upperPlain
#' @aliases getUpperPlain.CovarianceSELocalized
#' @aliases setUpperPlain.CovarianceSELocalized
#' @S3method getUpperPlain CovarianceSELocalized
#' @export getUpperPlain getUpperPlain.CovarianceSELocalized
#'
#' @param U A (named) vector of new upper bounds (we ONLY use ones which are
#'    named, and whose names match up with names of parameters.)
#' @param ... Not used.
#'
#' @usage CovarianceSELocalized$upperPlain
#'
#' @return The upper bounds for the parameters for this covariance function, but
#'    with names undecorated by its id.
#'
#' @seealso \code{\link{getLowerPlain.CovarianceSELocalized}}
#' @seealso \code{\link{CovarianceSELocalized}}
setMethodS3("getUpperPlain", "CovarianceSELocalized", conflict="quiet",
  function(this, ...) {
    U <- c(this$.ell.bounds[2], this$.sigma.f.bounds[2], this$.X.L.bounds[2],
      this$.X.R.bounds[2])
    names(U) <- getParamNamesPlain(this)
    return (U)
  })
setMethodS3("setUpperPlain", "CovarianceSELocalized", conflict="quiet",
  function(this, U, ...) {
    if (length(U) < 1) {
      return(invisible(this))
    }
    posdef.names <- c('ell', 'sigma.f')
    U[posdef.names] <- pmax(U[posdef.names], 0)
    # Adjust lower bounds to make way for the new values of U
    U.change <- this$PushLowerBounds(L.max=U)

    U.vals <- this$getUpperPlain()
    U.vals[names(U.change)] <- U.change[names(U.change)]
    this$.ell.bounds[2] <- U.vals["ell"]
    this$.sigma.f.bounds[2] <- U.vals["sigma.f"]
    this$.X.L.bounds[2] <- U.vals["X.L"]
    this$.X.R.bounds[2] <- U.vals["X.R"]
    this$ClampParams(warn=TRUE)
    return (this)
  })

#' Localized Squared-Exponential Covariance matrix
#'
#' Calculates a covariance matrix for the localized squared-exponential
#' covariance function.
#'
#' @S3method K.specific CovarianceSELocalized
#' @export K.specific K.specific.CovarianceSELocalized
#' @name K.specific.CovarianceSELocalized
#'
#' @param X  X-values for the input points (i.e., where we have data)
#' @param X.out  X-values for the points desired to predict
#' @param ... Not used.
#'
#' @return The covariance matrix taking \code{X} into \code{X.out}, based on the
#'    parameter values in \code{this}.
#'
#' @seealso \code{\link{CovarianceSELocalized}}
setMethodS3("K.specific", "CovarianceSELocalized", conflict="quiet",
  function(this, X, X.out=X, ...) {
    X.dist <<- DistanceMatrix(X=X, X.out=X.out)
    p <- this$getParamsPlain()
    mask <<- localize.mask(X=X, X.L=p["X.L"], X.R=p["X.R"], ell=p["ell"])
    mask.out <<- localize.mask(X=X.out, X.L=p["X.L"], X.R=p["X.R"], ell=p["ell"])
    return (outer(mask.out, mask) *
      (p["sigma.f"] ^ 2) * exp(-0.5 * (X.dist / p["ell"]) ^ 2))
  })

#' Element-wise derivatives of Covariance matrix
#'
#' Calculate the element-wise derivative of \code{KInIn}, with respect to the
#' parameter whose (plain) name is \code{param}.
#'
#' @S3method KDerivImplementation CovarianceSELocalized
#' @export KDerivImplementation KDerivImplementation.CovarianceSELocalized
#' @name KDerivImplementation.CovarianceSELocalized
#'
#' @param d  The Dataset whose X-values determine KInIn.
#' @param param  The (plain) name of the parameter with respect to which we're
#'    differentiating.
#' @param ... Not used.
#'
#' @return A matrix whose elements are the derivatives of the corresponding
#'    elements in KInIn, with respect to the parameter \code{param}.
#'
#' @seealso \code{\link{CovarianceSELocalized}}
setMethodS3("KDerivImplementation", "CovarianceSELocalized", conflict="quiet",
  function(this, d, param, ...) {
    p <- this$paramsPlain
    K.unmasked <- K.specific.CovarianceSE(this=this, X=d$X)
    dist.L <- as.vector((d$X - p["X.L"]) / p["ell"])
    dist.R <- as.vector((d$X - p["X.R"]) / p["ell"])
    mask <- localize.mask(X=d$X, X.L=p["X.L"], X.R=p["X.R"], ell=p["ell"])
    if (param == "ell") {
      d.mask.d.ell <- -(0.5 / p["ell"]) * (
        dist.L * standard.mask.deriv(dist.L) -
        dist.R * standard.mask.deriv(dist.R))
      mask.deriv.matrix <<- outer(mask, d.mask.d.ell)
      mask.deriv.matrix <- mask.deriv.matrix + t(mask.deriv.matrix)
      K.deriv <- (mask.deriv.matrix * K.unmasked
        + this$KInIn(d=d) * (DistanceMatrix(X=d$X) ^ 2) / (p["ell"] ^ 3))
    } else if (param == "sigma.f") {
      K.deriv <- 2 * this$KInIn(d=d) / p["sigma.f"]
    } else if (param == "X.L") {
      d.mask.d.X.L <- -(0.5 / p["ell"]) * standard.mask.deriv(dist.L)
      mask.deriv.matrix <- outer(mask, d.mask.d.X.L)
      mask.deriv.matrix <- mask.deriv.matrix + t(mask.deriv.matrix)
      K.deriv <- mask.deriv.matrix * K.unmasked
    } else if (param == "X.R") {
      d.mask.d.X.R <-  (0.5 / p["ell"]) * standard.mask.deriv(dist.R)
      mask.deriv.matrix <- outer(mask, d.mask.d.X.R)
      mask.deriv.matrix <- mask.deriv.matrix + t(mask.deriv.matrix)
      K.deriv <- mask.deriv.matrix * K.unmasked
    } else {
      K.deriv <- matrix(0, nrow=d$n, ncol=d$n)
    }
    return (K.deriv)
  })

#' Localized SE variance at each point
#'
#' Calculate the localized SE variance of the points at X: i.e., the a priori
#' uncertainty at each point.
#'
#' @S3method Variance CovarianceSELocalized
#' @export Variance Variance.CovarianceSELocalized
#' @name Variance.CovarianceSELocalized
#'
#' @param X  The points we want to know the localized SE variance at.
#' @param ... Not used.
#'
#' @return A numeric vector of the same length as X, with the corresponding
#'    localized SE variance.
#'
#' @seealso \code{\link{CovarianceSELocalized}}
setMethodS3("Variance", "CovarianceSELocalized", conflict="quiet",
  function(this, X, ...) {
    p <- this$paramsPlain
    mask <- localize.mask(X=X, X.L=p["X.L"], X.R=p["X.R"], ell=p["ell"])
    return ((p["sigma.f"] * mask) ^ 2)
  })

