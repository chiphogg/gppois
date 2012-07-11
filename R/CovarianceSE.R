############################################################################/**
# @RdocClass CovarianceSE
#
# @title "(S)quared-(E)xponential Covariance"
#
# \description{
#   The standard squared-exponential covariance.  Governed by two parameters: a
#   horizontal and a vertical lengthscale.
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
#   \item{ell.bounds}{(numeric) The range of values which \code{ell} might
#      assume.}
#   \item{sigma.f.bounds}{(numeric) The range of values which \code{sigma.f}
#      might assume.}
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
setConstructorS3("CovarianceSE", function(..., id="SE",
    ell=NA, sigma.f=NA, ell.bounds=NA, sigma.f.bounds=NA) {
    ell.good <- InitializeBoundedQuantity(ok.range=c(0, Inf),
      quantity=ell, bounds=ell.bounds, logspace=TRUE)
    sigma.f.good <- InitializeBoundedQuantity(ok.range=c(0, Inf),
      quantity=sigma.f, bounds=sigma.f.bounds, logspace=TRUE)
    # Construct the CovarianceSE object:
    extend(Covariance(..., id=id), "CovarianceSE",
      .ell            = ell.good$quantity,
      .ell.bounds     = ell.good$bounds,
      .sigma.f        = sigma.f.good$quantity,
      .sigma.f.bounds = sigma.f.good$bounds)
})

#' Names of "scale"-type parameters
#'
#' A character vector of names, indicating which parameters are considered to be
#' "scale" parameters.  (Read \dQuote{Optimization mode} section of
#' \code{\link{getParams.Covariance}} to see what this means.)
#'
#' @name getLogspaceNames.CovarianceSE
#' @aliases CovarianceSE$logspaceNames getLogspaceNames.CovarianceSE
#' @S3method getLogspaceNames CovarianceSE
#' @export getLogspaceNames getLogspaceNames.CovarianceSE
#'
#' @param ... Not used.
#'
#' @usage CovarianceSE$logspaceNames
#'
#' @return Names of parameters to be optimized in logspace.
#'
#' @seealso \code{\link{CovarianceSE}}
setMethodS3("getLogspaceNames", "CovarianceSE", conflict="quiet",
  function(this, ...) {
    return (c("ell", "sigma.f"))
  })

#' Basenames of parameters
#'
#' Gives the "basenames" (i.e. names undecorated by the id string) of the
#' parameters.
#'
#' @name getParamNamesPlain.CovarianceSE
#' @aliases CovarianceSE$paramNamesPlain getParamNamesPlain.CovarianceSE
#' @S3method getParamNamesPlain CovarianceSE
#' @export getParamNamesPlain getParamNamesPlain.CovarianceSE
#'
#' @param ... Not used.
#'
#' @usage CovarianceSE$paramNamesPlain
#'
#' @return The basenames of the parameters.
#'
#' @seealso \code{\link{CovarianceSE}}
setMethodS3("getParamNamesPlain", "CovarianceSE", conflict="quiet",
  function(this, ...) {
    return (c("ell", "sigma.f"))
  })

#' Parameter values with plain names
#'
#' Gives a vector of parameter values, whose names are NOT decorated by the id
#' of this Covariance object.
#'
#' @name getParamsPlain.CovarianceSE
#' @aliases CovarianceSE$paramsPlain getParamsPlain.CovarianceSE
#' @S3method getParamsPlain CovarianceSE
#' @export getParamsPlain getParamsPlain.CovarianceSE
#'
#' @param ... Not used.
#'
#' @usage CovarianceSE$paramsPlain
#'
#' @return The parameters for this covariance function, but with names
#'    undecorated by its id.
#'
#' @seealso \code{\link{setParamsPlain.Covariance}}
#' @seealso \code{\link{CovarianceSE}}
setMethodS3("getParamsPlain", "CovarianceSE", conflict="quiet",
  function(this, ...) {
    p <- c(this$.ell, this$.sigma.f)
    names(p) <- getParamNamesPlain(this)
    return (p)
})
setMethodS3("paramsPlainImplementation", "CovarianceSE", conflict="quiet",
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
    return (invisible(this))
  })

#' Lower bounds for params, with plain names
#'
#' Gives a vector of lower bounds for the parameter values, whose names are NOT
#' decorated by the id of this Covariance object.
#'
#' @name getLowerPlain.CovarianceSE
#' @aliases CovarianceSE$lowerPlain getLowerPlain.CovarianceSE setLowerPlain.CovarianceSE
#' @S3method getLowerPlain CovarianceSE
#' @export getLowerPlain getLowerPlain.CovarianceSE
#'
#' @param L A (named) vector of new lower bounds (we ONLY use ones which are
#'    named, and whose names match up with names of parameters.)
#' @param ... Not used.
#'
#' @usage CovarianceSE$lowerPlain
#'
#' @return The lower bounds for the parameters for this covariance function, but
#'    with names undecorated by its id.
#'
#' @seealso \code{\link{getUpperPlain.CovarianceSE}}
#' @seealso \code{\link{CovarianceSE}}
setMethodS3("getLowerPlain", "CovarianceSE", conflict="quiet",
  function(this, ...) {
    L <- c(this$.ell.bounds[1], this$.sigma.f.bounds[1])
    names(L) <- getParamNamesPlain(this)
    return (L)
})
setMethodS3("setLowerPlain", "CovarianceSE", conflict="quiet",
  function(this, L, ...) {
    if (length(L) < 1) {
      return(invisible(this))
    }
    L.posdef <- pmax(L, 0)  # SE has no possibly-negative parameters
    # Adjust upper bounds to make way for the new values of L
    L.change <- PushUpperBounds(this, U.min=L.posdef)
    L.vals <- this$getLowerPlain()
    L.vals[names(L.change)] <- L.change[names(L.change)]
    this$.ell.bounds[1] <- L.vals["ell"]
    this$.sigma.f.bounds[1] <- L.vals["sigma.f"]
    ClampParams(this, warn=TRUE)
    return (this)
  })

#' Upper bounds for params, with plain names
#'
#' Gives a vector of upper bounds for the parameter values, whose names are NOT
#' decorated by the id of this Covariance object.
#'
#' @name getUpperPlain.CovarianceSE
#' @aliases CovarianceSE$upperPlain getUpperPlain.CovarianceSE setUpperPlain.CovarianceSE
#' @S3method getUpperPlain CovarianceSE
#' @export getUpperPlain getUpperPlain.CovarianceSE
#'
#' @param U A (named) vector of new upper bounds (we ONLY use ones which are
#'    named, and whose names match up with names of parameters.)
#' @param ... Not used.
#'
#' @usage CovarianceSE$upperPlain
#'
#' @return The upper bounds for the parameters for this covariance function, but
#'    with names undecorated by its id.
#'
#' @seealso \code{\link{getLowerPlain.CovarianceSE}}
#' @seealso \code{\link{CovarianceSE}}
setMethodS3("getUpperPlain", "CovarianceSE", conflict="quiet",
  function(this, ...) {
    U <- c(this$.ell.bounds[2], this$.sigma.f.bounds[2])
    names(U) <- getParamNamesPlain(this)
    return (U)
  })
setMethodS3("setUpperPlain", "CovarianceSE", conflict="quiet",
  function(this, U, ...) {
    if (length(U) < 1) {
      return(invisible(this))
    }
    U.posdef <- pmax(U, 0)  # SE has no possibly-negative parameters
    # Adjust lower bounds to make way for the new values of U
    U.change <- PushLowerBounds(this, L.max=U.posdef)

    U.vals <- this$getUpperPlain()
    U.vals[names(U.change)] <- U.change[names(U.change)]
    this$.ell.bounds[2] <- U.vals["ell"]
    this$.sigma.f.bounds[2] <- U.vals["sigma.f"]
    ClampParams(this, warn=TRUE)
    return (this)
  })

#' Squared-Exponential Covariance matrix
#'
#' Calculates a covariance matrix for the squared-exponential covariance
#' function.
#'
#' @S3method K.specific CovarianceSE
#' @export K.specific K.specific.CovarianceSE
#' @name K.specific.CovarianceSE
#'
#' @param X  X-values for the input points (i.e., where we have data)
#' @param X.out  X-values for the points desired to predict
#' @param ... Not used.
#'
#' @return The covariance matrix taking \code{X} into \code{X.out}, based on the
#'    parameter values in \code{this}.
#'
#' @seealso \code{\link{CovarianceSE}}
setMethodS3("K.specific", "CovarianceSE", conflict="quiet",
  function(this, X, X.out=X, ...) {
    X.dist <- DistanceMatrix(X=X, X.out=X.out)
    p <- this$getParamsPlain()
    return (p["sigma.f"] ^ 2 * exp(-0.5 * (X.dist / p["ell"]) ^ 2))
  })

#' Element-wise derivatives of Covariance matrix
#'
#' Calculate the element-wise derivative of \code{KInIn}, with respect to the
#' parameter whose (plain) name is \code{param}.
#'
#' @S3method KDerivImplementation CovarianceSE
#' @export KDerivImplementation KDerivImplementation.CovarianceSE
#' @name KDerivImplementation.CovarianceSE
#'
#' @param d  The Dataset whose X-values determine KInIn.
#' @param param  The (plain) name of the parameter with respect to which we're
#'    differentiating.
#' @param ... Not used.
#'
#' @return A matrix whose elements are the derivatives of the corresponding
#'    elements in KInIn, with respect to the parameter \code{param}.
#'
#' @seealso \code{\link{CovarianceSE}}
setMethodS3("KDerivImplementation", "CovarianceSE", conflict="quiet",
  function(this, d, param, ...) {
    p <- this$paramsPlain
    if (param == "ell") {
      K.deriv <- this$KInIn(d=d) * (DistanceMatrix(X=d$X) ^ 2) / (p["ell"] ^ 3)
    } else if (param == "sigma.f") {
      K.deriv <- 2 * this$KInIn(d=d) / p["sigma.f"]
    } else {
      K.deriv <- matrix(0, nrow=d$n, ncol=d$n)
    }
    return (K.deriv)
  })

#' SE variance at each point
#'
#' Calculate the SE variance of the points at X: i.e., the a priori uncertainty
#' at each point.
#'
#' @S3method Variance CovarianceSE
#' @export Variance Variance.CovarianceSE
#' @name Variance.CovarianceSE
#'
#' @param X  The points we want to know the SE variance at.
#' @param ... Not used.
#'
#' @return A numeric vector of the same length as X, with the corresponding
#'    SE variance.
#'
#' @seealso \code{\link{CovarianceSE}}
setMethodS3("Variance", "CovarianceSE", conflict="quiet",
  function(this, X, ...) {
    return (rep(this$.sigma.f ^ 2, NumPoints(X)))
  })

