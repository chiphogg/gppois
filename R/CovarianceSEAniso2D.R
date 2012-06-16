############################################################################/**
# @RdocClass CovarianceSEAniso2D
#
# @title "2D Anisotropic SE Covariance"
#
# \description{
#   A 2D squared-exponential covariance whose eigenvalues of the covariance
#   matrix are not assumed identical.  In other words, features might vary more
#   rapidly in one direction than in the orthogonal direction.
#
#   @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#   \item{id}{(character) A string to identify this covariance object.}
#   \item{ell.1}{(numeric) One characteristic horizontal scale for features in
#      the functions.}
#   \item{ell.2}{(numeric) Another characteristic horizontal scale for features
#      in the functions.}
#   \item{theta.1}{(numeric) The angle of the ell.1 direction.}
#   \item{sigma.f}{(numeric) A characteristic vertical scale for features in
#      functions being modeled.}
#   \item{ell.1.bounds}{(numeric) The range of values which \code{ell.1} might
#      assume.}
#   \item{ell.2.bounds}{(numeric) The range of values which \code{ell.2} might
#      assume.}
#   \item{theta.1.bounds}{(numeric) The range of values which \code{theta.1}
#      might assume.}
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
#     \item{ell.1}{(Scale parameter) One horizontal feature lengthscale.}
#     \item{ell.2}{(Scale parameter) Another horizontal feature lengthscale.}
#     \item{theta.1}{The angle of the ell.1 axis}
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
setConstructorS3("CovarianceSEAniso2D",
  function(..., id="Aniso2D", ell.1=NA, ell.2=NA, theta.1=NA, sigma.f=NA,
    ell.1.bounds=NA, ell.2.bounds=NA, theta.1.bounds=NA, sigma.f.bounds=NA) {
    pos.def.range <- c(0, Inf)
    ell.1.good <- InitializeBoundedQuantity(ok.range=pos.def.range,
      quantity=ell.1, bounds=ell.1.bounds, logspace=TRUE)
    ell.2.good <- InitializeBoundedQuantity(ok.range=pos.def.range,
      quantity=ell.2, bounds=ell.2.bounds, logspace=TRUE)
    theta.1.good <- InitializeBoundedQuantity(ok.range=pi * c(-1, 1),
      quantity=theta.1, bounds=theta.1.bounds)
    sigma.f.good <- InitializeBoundedQuantity(ok.range=pos.def.range,
      quantity=sigma.f, bounds=sigma.f.bounds, logspace=TRUE)
    # Construct the CovarianceSEAniso2D object:
    extend(Covariance(..., id=id), "CovarianceSEAniso2D",
      .ell.1          = ell.1.good$quantity,
      .ell.1.bounds   = ell.1.good$bounds,
      .ell.2          = ell.2.good$quantity,
      .ell.2.bounds   = ell.2.good$bounds,
      .theta.1        = theta.1.good$quantity,
      .theta.1.bounds = theta.1.good$bounds,
      .sigma.f        = sigma.f.good$quantity,
      .sigma.f.bounds = sigma.f.good$bounds)
  })

#' Names of "scale"-type parameters
#'
#' A character vector of names, indicating which parameters are considered to be
#' "scale" parameters.  (Read \dQuote{Optimization mode} section of
#' \code{\link{getParams.Covariance}} to see what this means.)
#'
#' @name getLogspaceNames.CovarianceSEAniso2D
#' @aliases CovarianceSEAniso2D$logspaceNames
#' @aliases getLogspaceNames.CovarianceSEAniso2D
#' @S3method getLogspaceNames CovarianceSEAniso2D
#'
#' @param ... Not used.
#'
#' @usage CovarianceSEAniso2D$logspaceNames
#'
#' @return Names of parameters to be optimized in logspace.
#'
#' @seealso \code{\link{CovarianceSEAniso2D}}
setMethodS3("getLogspaceNames", "CovarianceSEAniso2D", conflict="quiet",
  function(this, ...) {
    return (c("ell.1", "ell.2", "sigma.f"))
  })

#' Basenames of parameters
#'
#' Gives the "basenames" (i.e. names undecorated by the id string) of the
#' parameters.
#'
#' @name getParamNamesPlain.CovarianceSEAniso2D
#' @aliases CovarianceSEAniso2D$paramNamesPlain
#' @aliases getParamNamesPlain.CovarianceSEAniso2D
#' @S3method getParamNamesPlain CovarianceSEAniso2D
#'
#' @param ... Not used.
#'
#' @usage CovarianceSEAniso2D$paramNamesPlain
#'
#' @return The basenames of the parameters.
#'
#' @seealso \code{\link{CovarianceSEAniso2D}}
setMethodS3("getParamNamesPlain", "CovarianceSEAniso2D", conflict="quiet",
  function(this, ...) {
    return (c("ell.1", "ell.2", "theta.1", "sigma.f"))
  })

#' Parameter values with plain names
#'
#' Gives a vector of parameter values, whose names are NOT decorated by the id
#' of this Covariance object.
#'
#' @name getParamsPlain.CovarianceSEAniso2D
#' @aliases CovarianceSEAniso2D$paramsPlain getParamsPlain.CovarianceSEAniso2D
#' @S3method getParamsPlain CovarianceSEAniso2D
#'
#' @param ... Not used.
#'
#' @usage CovarianceSEAniso2D$paramsPlain
#'
#' @return The parameters for this covariance function, but with names
#'    undecorated by its id.
#'
#' @seealso \code{\link{setParamsPlain.Covariance}}
#' @seealso \code{\link{CovarianceSEAniso2D}}
setMethodS3("getParamsPlain", "CovarianceSEAniso2D", conflict="quiet",
  function(this, ...) {
    p <- c(this$.ell.1, this$.ell.2, this$.theta.1, this$.sigma.f)
    names(p) <- getParamNamesPlain(this)
    return (p)
})
setMethodS3("paramsPlainImplementation", "CovarianceSEAniso2D", conflict="quiet",
  private=TRUE,
  function(this, p, ...) {
    # Sets any values in 'p' which match our parameter names, without worrying
    # about lower and upper bounds (the function which calls this has the job
    # of worrying about these!).
    p.old <- this$getParamsPlain()
    to.change <- names(p)[which(names(p) %in% names(p.old))]
    p.old[to.change] <- p[to.change]
    this$.ell.1 <- p.old["ell.1"]
    this$.ell.2 <- p.old["ell.2"]
    this$.theta.1 <- p.old["theta.1"]
    this$.sigma.f <- p.old["sigma.f"]
    return (invisible(this))
  })

#' Lower bounds for params, with plain names
#'
#' Gives a vector of lower bounds for the parameter values, whose names are NOT
#' decorated by the id of this Covariance object.
#'
#' @name getLowerPlain.CovarianceSEAniso2D
#' @aliases CovarianceSEAniso2D$lowerPlain getLowerPlain.CovarianceSEAniso2D
#' @aliases setLowerPlain.CovarianceSEAniso2D
#' @S3method getLowerPlain CovarianceSEAniso2D
#'
#' @param L A (named) vector of new lower bounds (we ONLY use ones which are
#'    named, and whose names match up with names of parameters.)
#' @param ... Not used.
#'
#' @usage CovarianceSEAniso2D$lowerPlain
#'
#' @return The lower bounds for the parameters for this covariance function,
#'    but with names undecorated by its id.
#'
#' @seealso \code{\link{getUpperPlain.CovarianceSEAniso2D}}
#' @seealso \code{\link{CovarianceSEAniso2D}}
setMethodS3("getLowerPlain", "CovarianceSEAniso2D", conflict="quiet",
  function(this, ...) {
    L <- c(this$.ell.1.bounds[1], this$.ell.2.bounds[1],
      this$.theta.1.bounds[1], this$.sigma.f.bounds[1])
    names(L) <- getParamNamesPlain(this)
    return (L)
})
setMethodS3("setLowerPlain", "CovarianceSEAniso2D", conflict="quiet",
  function(this, L, ...) {
    # Adjust upper bounds to make way for the new values of L
    L.change <- this$PushUpperBounds(U.min=L)

    L.vals <- this$getLowerPlain()
    L.vals[names(L.change)] <- L.change[names(L.change)]
    this$.ell.1.bounds[1] <- L.vals["ell.1"]
    this$.ell.2.bounds[1] <- L.vals["ell.2"]
    this$.theta.1.bounds[1] <- L.vals["theta.1"]
    this$.sigma.f.bounds[1] <- L.vals["sigma.f"]
    this$ClampParams(warn=TRUE)
    return (this)
  })

#' Upper bounds for params, with plain names
#'
#' Gives a vector of upper bounds for the parameter values, whose names are NOT
#' decorated by the id of this Covariance object.
#'
#' @name getUpperPlain.CovarianceSEAniso2D
#' @aliases CovarianceSEAniso2D$upperPlain getUpperPlain.CovarianceSEAniso2D
#' @aliases setUpperPlain.CovarianceSEAniso2D
#' @S3method getUpperPlain CovarianceSEAniso2D
#'
#' @param U A (named) vector of new upper bounds (we ONLY use ones which are
#'    named, and whose names match up with names of parameters.)
#' @param ... Not used.
#'
#' @usage CovarianceSEAniso2D$upperPlain
#'
#' @return The upper bounds for the parameters for this covariance function, but
#'    with names undecorated by its id.
#'
#' @seealso \code{\link{getLowerPlain.CovarianceSEAniso2D}}
#' @seealso \code{\link{CovarianceSEAniso2D}}
setMethodS3("getUpperPlain", "CovarianceSEAniso2D", conflict="quiet",
  function(this, ...) {
    U <- c(this$.ell.1.bounds[2], this$.ell.2.bounds[2],
      this$.theta.1.bounds[2], this$.sigma.f.bounds[2])
    names(U) <- getParamNamesPlain(this)
    return (U)
  })
setMethodS3("setUpperPlain", "CovarianceSEAniso2D", conflict="quiet",
  function(this, U, ...) {
    # Adjust lower bounds to make way for the new values of U
    U.change <- this$PushLowerBounds(L.max=U)

    U.vals <- this$getUpperPlain()
    U.vals[names(U.change)] <- U.change[names(U.change)]
    this$.ell.1.bounds[2] <- U.vals["ell.1"]
    this$.ell.2.bounds[2] <- U.vals["ell.2"]
    this$.theta.1.bounds[2] <- U.vals["theta.1"]
    this$.sigma.f.bounds[2] <- U.vals["sigma.f"]
    this$ClampParams(warn=TRUE)
    return (this)
  })

#' Anisotropic 2D SE Covariance matrix
#'
#' Calculates a covariance matrix for the anisotropic 2-dimensional
#' squared-exponential covariance function.
#'
#' @S3method K.specific CovarianceSEAniso2D
#' @name K.specific.CovarianceSEAniso2D
#'
#' @param X  X-values for the input points (i.e., where we have data)
#' @param X.out  X-values for the points desired to predict
#' @param ... Not used.
#'
#'  @return The covariance matrix taking \code{X} into \code{X.out}, based on
#'    the parameter values in \code{this}.
#'
#' @seealso \code{\link{CovarianceSEAniso2D}}
setMethodS3("K.specific", "CovarianceSEAniso2D", conflict="quiet",
  function(this, X, X.out=X, ...) {
    p <- this$getParamsPlain()
    ell <- c(p["ell.1"], p["ell.2"])

    # Rotate and scale the X-values
    cos.1 <- cos(p["theta.1"])
    sin.1 <- sin(p["theta.1"])
    R <- matrix(c(cos.1, sin.1, -sin.1, cos.1), ncol=2)
    X.rot     <- (X     %*% R) / matrix(
      rep(ell, each=NumPoints(X    )), ncol=2)
    X.out.rot <- (X.out %*% R) / matrix(
      rep(ell, each=NumPoints(X.out)), ncol=2)

    X.dist <- DistanceMatrix(X=X.rot, X.out=X.out.rot)
    return (p["sigma.f"] ^ 2 * exp(-0.5 * (X.dist ^ 2)))
  })

#' Element-wise derivatives of Covariance matrix
#'
#' Calculate the element-wise derivative of \code{KInIn}, with respect to the
#' parameter whose (plain) name is \code{param}.
#'
#' @S3method KDerivImplementation CovarianceSEAniso2D
#' @name KDerivImplementation.CovarianceSEAniso2D
#'
#' @param d  The Dataset whose X-values determine KInIn.
#' @param param  The (plain) name of the parameter with respect to which we're
#'    differentiating.
#' @param ... Not used.
#'
#' @return A matrix whose elements are the derivatives of the corresponding
#'    elements in KInIn, with respect to the parameter \code{param}.
#'
#' @seealso \code{\link{CovarianceSEAniso2D}}
setMethodS3("KDerivImplementation", "CovarianceSEAniso2D", conflict="quiet",
  function(this, d, param, ...) {
    Delta.1 <- outer(d$X[, 1], d$X[, 1], '-')
    Delta.2 <- outer(d$X[, 2], d$X[, 2], '-')
    p <- this$paramsPlain
    c.theta <- cos(p["theta.1"])
    s.theta <- sin(p["theta.1"])
    if (param == "ell.1") {
      K.deriv <- (this$KInIn(d=d) / p["ell.1"]) * (
        (c.theta * Delta.1 - s.theta * Delta.2) / p["ell.1"]) ^ 2
    } else if (param == "ell.2") {
      K.deriv <- (this$KInIn(d=d) / p["ell.2"]) * (
        (s.theta * Delta.1 + c.theta * Delta.2) / p["ell.2"]) ^ 2
    } else if (param == "sigma.f") {
      K.deriv <- 2 * this$KInIn(d=d) / p["sigma.f"]
    } else if (param == "theta.1") {
      K.deriv <- (this$KInIn(d=d) 
        * (1 / (p["ell.1"] ^ 2) - 1 / (p["ell.2"] ^ 2)) 
        * (0.5 * sin(2 * p["theta.1"]) * (Delta.1 ^ 2 - Delta.2 ^ 2)
          + cos(2 * p["theta.1"]) * Delta.1 * Delta.2)
        )
    } else {
      K.deriv <- matrix(0, nrow=d$n, ncol=d$n)
    }
    return (K.deriv)
  })

#' Anisotropic 2D SE variance at each point
#'
#' Calculate the Anisotropic 2D SE variance of the points at X: i.e., the
#' a priori uncertainty at each point.
#'
#' @S3method Variance CovarianceSEAniso2D
#' @name Variance.CovarianceSEAniso2D
#'
#' @param X  The points we want to know the anisotropic 2D SE variance at.
#' @param ... Not used.
#'
#' @return A numeric vector of the same length as X, with the corresponding
#'    anisotropic 2D SE variance.
#'
#' @seealso \code{\link{CovarianceSEAniso2D}}
setMethodS3("Variance", "CovarianceSEAniso2D", conflict="quiet",
  function(this, X, ...) {
    return (rep(this$.sigma.f ^ 2, NumPoints(X)))
  })

