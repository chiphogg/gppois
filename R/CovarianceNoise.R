############################################################################/**
# @RdocClass CovarianceNoise
#
# @title "Covariance describing random noise"
#
# \description{
#   This subclass lets us treat noise in a unified way within our Model.
#
#   @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#   \item{id}{(character) A string to identify this covariance object.}
#   \item{sigma}{(numeric) The default value of the noise.}
#   \item{sigma.bounds}{(numeric) The range of values which \code{sigma} might assume.}
#   \item{...}{Not used.}
# }
#
# \section{Covariance Parameters}{
#   This section lists the fit parameters corresponding to this type of
#   Covariance.  Any parameters marked as \dQuote{(Scale parameter)} will be
#   optimized in log-space, consistent with the Jeffreys prior.
#
#   \describe{
#     \item{sigma}{(Scale parameter) The magnitude (standard deviation) of the
#        noise.}
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
setConstructorS3("CovarianceNoise",
  function(id="noise", sigma=NA, sigma.bounds=NA, ...) {
    # Ideas here are the same as for the CovarianceSE constructor, but simpler
    # because there is one fewer parameter.
    sigma.good <- InitializeBoundedQuantity(ok.range=c(0, Inf),
      quantity=sigma, bounds=sigma.bounds, logspace=TRUE)
    # Construct the CovarianceNoise object:
    extend(Covariance(..., id=id), "CovarianceNoise",
      .sigma        = sigma.good$quantity,
      .sigma.bounds = sigma.good$bounds)
  })

#' Names of "scale"-type parameters
#'
#' A character vector of names, indicating which parameters are considered to be
#' "scale" parameters.  (Read \dQuote{Optimization mode} section of
#' \code{\link{getParams.Covariance}} to see what this means.)
#'
#' @name getLogspaceNames.CovarianceNoise
#' @aliases CovarianceNoise$logspaceNames getLogspaceNames.CovarianceNoise
#' @S3method getLogspaceNames CovarianceNoise
#'
#' @param ... Not used.
#'
#' @usage CovarianceNoise$logspaceNames
#'
#' @return Names of parameters to be optimized in logspace.
#'
#' @seealso \code{\link{CovarianceNoise}}
setMethodS3("getLogspaceNames", "CovarianceNoise", conflict="quiet",
  function(this, ...) {
    return (c("sigma"))
  })

#' Basenames of parameters
#'
#' Gives the "basenames" (i.e. names undecorated by the id string) of the
#' parameters.
#'
#' @name getParamNamesPlain.CovarianceNoise
#' @aliases CovarianceNoise$paramNamesPlain getParamNamesPlain.CovarianceNoise
#' @S3method getParamNamesPlain CovarianceNoise
#'
#' @param ... Not used.
#'
#' @usage CovarianceNoise$paramNamesPlain
#'
#' @return The basenames of the parameters.
#'
#' @seealso \code{\link{CovarianceNoise}}
setMethodS3("getParamNamesPlain", "CovarianceNoise", conflict="quiet",
  function(this, ...) {
    return (c("sigma"))
  })

#' Parameter values with plain names
#'
#' Gives a vector of parameter values, whose names are NOT decorated by the id
#' of this Covariance object.
#'
#' @name getParamsPlain.CovarianceNoise
#' @aliases CovarianceNoise$paramsPlain getParamsPlain.CovarianceNoise
#' @S3method getParamsPlain CovarianceNoise
#'
#' @param ... Not used.
#'
#' @usage CovarianceNoise$paramsPlain
#'
#' @return The parameters for this covariance function, but with names
#'    undecorated by its id.
#'
#' @seealso \code{\link{setParamsPlain.Covariance}}
#' @seealso \code{\link{CovarianceNoise}}
setMethodS3("getParamsPlain", "CovarianceNoise", conflict="quiet",
  function(this, ...) {
    p <- c(this$.sigma)
    names(p) <- this$getParamNamesPlain()
    return (p)
  })
setMethodS3("paramsPlainImplementation", "CovarianceNoise", conflict="quiet",
  private=TRUE,
  function(this, p, ...) {
    # Sets any values in 'p' which match our parameter names, without worrying
    # about lower and upper bounds (the function which calls this has the job
    # of worrying about these!).
    this$.sigma <- p["sigma"]
    return (invisible(this))
  })

#' Lower bounds for params, with plain names
#'
#' Gives a vector of lower bounds for the parameter values, whose names are NOT
#' decorated by the id of this Covariance object.
#'
#' @name getLowerPlain.CovarianceNoise
#' @aliases CovarianceNoise$lowerPlain getLowerPlain.CovarianceNoise setLowerPlain.CovarianceNoise
#' @S3method getLowerPlain CovarianceNoise
#'
#' @param L A (named) vector of new lower bounds (we ONLY use ones which are
#'    named, and whose names match up with names of parameters.)
#' @param ... Not used.
#'
#' @usage CovarianceNoise$lowerPlain
#'
#' @return The lower bounds for the parameters for this covariance function, but
#'    with names undecorated by its id.
#'
#' @seealso \code{\link{getUpperPlain.CovarianceNoise}}
#' @seealso \code{\link{CovarianceNoise}}
setMethodS3("getLowerPlain", "CovarianceNoise", conflict="quiet",
  function(this, ...) {
    L <- c(this$.sigma.bounds[1])
    names(L) <- getParamNamesPlain(this)
    return (L)
})
setMethodS3("setLowerPlain", "CovarianceNoise", conflict="quiet",
  function(this, L, ...) {
    if (length(L) < 1) {
      return(invisible(this))
    }
    L.posdef <- pmax(L, 0)  # Noise cannot be negative

    # Adjust upper bounds to make way for the new values of L
    L.change <- this$PushUpperBounds(U.min=L.posdef)

    L.vals <- this$getLowerPlain()
    L.vals[names(L.change)] <- L.change[names(L.change)]
    this$.sigma.bounds[1] <- L.vals["sigma"]
    this$ClampParams(warn=TRUE)
    return (this)
  })

#' Upper bounds for params, with plain names
#'
#' Gives a vector of upper bounds for the parameter values, whose names are NOT
#' decorated by the id of this Covariance object.
#'
#' @name getUpperPlain.CovarianceNoise
#' @aliases CovarianceNoise$upperPlain getUpperPlain.CovarianceNoise setUpperPlain.CovarianceNoise
#' @S3method getUpperPlain CovarianceNoise
#'
#' @param U A (named) vector of new upper bounds (we ONLY use ones which are
#'    named, and whose names match up with names of parameters.)
#' @param ... Not used.
#'
#' @usage CovarianceNoise$upperPlain
#'
#' @return The upper bounds for the parameters for this covariance function, but
#'    with names undecorated by its id.
#'
#' @seealso \code{\link{getLowerPlain.CovarianceNoise}}
#' @seealso \code{\link{CovarianceNoise}}
setMethodS3("getUpperPlain", "CovarianceNoise", conflict="quiet",
  function(this, ...) {
    U <- c(this$.sigma.bounds[2])
    names(U) <- getParamNamesPlain(this)
    return (U)
  })
setMethodS3("setUpperPlain", "CovarianceNoise", conflict="quiet",
  function(this, U, ...) {
    if (length(U) < 1) {
      return(invisible(this))
    }
    U.posdef <- pmax(U, 0)  # SE has no possibly-negative parameters

    # Adjust lower bounds to make way for the new values of U
    U.change <- this$PushLowerBounds(L.max=U.posdef)

    U.vals <- this$getUpperPlain()
    U.vals[names(U.change)] <- U.change[names(U.change)]
    this$.sigma.bounds[2] <- U.vals["sigma"]
    this$ClampParams(warn=TRUE)
    return (this)
  })

#' Noise Covariance matrix
#'
#' Calculates a covariance matrix for the noise covariance specifically.
#'
#' @S3method K.specific CovarianceNoise
#' @name K.specific.CovarianceNoise
#'
#' @param X  X-values for the input points (i.e., where we have data)
#' @param X.out  X-values for the points desired to predict
#' @param ... Not used.
#'
#' @note This is very similar to a \code{\link{CovarianceSE}} in the limit as
#'    \code{ell} goes to zero.  However, there is a very important difference!
#'    Consider the case of multiple datapoints at the same X-value.
#'    CovarianceSE assigns a correlation of 1; CovarianceNoise assigns a
#'    correlation of 0.  Caution is required!
#'
#' @return The covariance matrix taking \code{X} into \code{X.out}, based on the
#'    parameter values in \code{this}.
#'
#' @seealso \code{\link{CovarianceNoise}}
setMethodS3("K.specific", "CovarianceNoise", conflict="quiet",
  function(this, X, X.out=NA, ...) {
    if (is.na(X.out) || identical(X, X.out)) {
      K <- (this$getParamsPlain()["sigma"] ^ 2) * diag(NumPoints(X))
    } else {
      K <- matrix(0, nrow=NumPoints(X.out), ncol=NumPoints(X))
    }
    return (K)
  })

#' Element-wise derivatives of Covariance matrix
#'
#' Calculate the element-wise derivative of \code{KInIn}, with respect to the
#' parameter whose (plain) name is \code{param}.
#'
#' @S3method KDerivImplementation CovarianceNoise
#' @name KDerivImplementation.CovarianceNoise
#'
#' @param d  The Dataset whose X-values determine KInIn.
#' @param param  The (plain) name of the parameter with respect to which we're
#'    differentiating.
#' @param ... Not used.
#'
#' @return A matrix whose elements are the derivatives of the corresponding
#'    elements in KInIn, with respect to the parameter \code{param}.
#'
#' @seealso \code{\link{CovarianceNoise}}
setMethodS3("KDerivImplementation", "CovarianceNoise", conflict="quiet",
  function(this, d, param, ...) {
    # Calculate the element-wise derivative of KInIn, with respect to the
    # parameter whose (plain) name is 'param'.
    #
    # Args:
    #   d:  The Dataset whose X-values determine KInIn.
    #   param:  The (plain) name of the parameter with respect to which we're
    #      differentiating .
    #
    # Returns:
    #   A matrix whose elements are the derivatives of the corresponding
    #   elements in KInIn, with respect to the parameter 'param'.
    if (param == "sigma") {
      K.deriv <- 2 * this$paramsPlain["sigma"] * diag(d$n)
    } else {
      K.deriv <- matrix(0, nrow=d$n, ncol=d$n)
    }
    return (K.deriv)
  })

#' Noise variance at each point
#'
#' Calculate the Noise variance of the points at X.
#'
#' @S3method Variance CovarianceNoise
#' @name Variance.CovarianceNoise
#'
#' @param X  The points we want to know the noise variance at.
#' @param ... Not used.
#'
#' @return A numeric vector of the same length as X, with the corresponding
#'    noise variance.
#'
#' @seealso \code{\link{CovarianceNoise}}
setMethodS3("Variance", "CovarianceNoise", conflict="quiet",
  # Calculate the Noise variance of the points at X.
  #
  # @param X  The points we want to know the noise variance at.
  #
  # Returns:
  #   A numeric vector of the same length as X, with the corresponding noise
  #   variance.
  function(this, X, ...) {
    return (rep(this$.sigma ^ 2, NumPoints(X)))
  })

