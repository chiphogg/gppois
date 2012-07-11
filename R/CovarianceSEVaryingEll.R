############################################################################/**
# @RdocClass CovarianceSEVaryingEll
#
# @title "Nonstationary Squared-Exponential Covariance"
#
# \description{
#   A \emph{nonstationary} squared-exponential covariance (meaning that the
#   parameters are \emph{no longer} independent of the covariates).  It is
#   governed by the two usual parameters (i.e., the horizontal and vertical
#   lengthscales \code{ell} and \code{sigma.f}), but \code{ell} now
#   depends on \code{X}.
#
#   Note that it cannot easily be optimized within the current paradigm!
#   Usually in \pkg{gppois}, a Model is trained on a Dataset by varying the
#   hyperparameters, which are assumed to be small and finite in number.  But
#   now \code{ell} is a continuous function -- there are infinitely many
#   hyperparameters!  Even if we only care about the denoised values of each
#   datapoint, we still have one hyperparameter per datapoint, which is far too
#   many for most optimizers.
#
#   The solution is to optimize it with a different paradigm: the \emph{focus
#   regions} approach (Hogg et al. 2012).  If \code{ell(X)} varies slowly, we
#   can break our function into pieces, perform local fits, and interpolate the
#   resulting values of \code{ell}.  This functionality is not yet integrated
#   into the main package.
#
#   @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#   \item{id}{(character) A string to identify this covariance object.}
#   \item{X.ell}{(numeric) The X-values where \code{ell} is specified.}
#   \item{ell}{(numeric) A series of local values for the characteristic
#      horizontal scale for features in functions being modeled.}
#   \item{sigma.f}{(numeric) A characteristic vertical scale for features in
#      functions being modeled.}
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
# \references{
#    Hogg, C., K. Mullen, and I. Levin (2012). A Bayesian approach for
#    denoising one-dimensional data. Journal of Applied Crystallography, 45(3),
#    pp. 471-481.
# }
#
# @author
#*/###########################################################################
setConstructorS3("CovarianceSEVaryingEll", function(..., id="SEVaryingEll",
    X.ell=NA, ell=NA, sigma.f=NA) {
    # Constructs a CovarianceSEVaryingEll object with the given parameter
    # values.
    #
    # Args:
    #   X.ell:  numeric vector of X-values corresponding to 'ell'
    #   ell: Samples of ell(X) (a characteristic lengthscale over which
    #      function values are correlated) at the locations in X.ell.
    #   sigma.f: The "vertical" lengthscale.
    #
    # Returns:
    #   CovarianceSEVaryingEll object with the given parameter values.

    # Construct the CovarianceSEVaryingEll object:
    extend(Covariance(..., id=id), "CovarianceSEVaryingEll",
      .X.ell   = X.ell,
      .ell     = ell,
      .sigma.f = sigma.f)
})

#' Names of "scale"-type parameters
#'
#' A character vector of names, indicating which parameters are considered to be
#' "scale" parameters.  (Read \dQuote{Optimization mode} section of
#' \code{\link{getParams.Covariance}} to see what this means.)
#'
#' @name getLogspaceNames.CovarianceSEVaryingEll
#' @aliases CovarianceSEVaryingEll$logspaceNames
#' @aliases getLogspaceNames.CovarianceSEVaryingEll
#' @S3method getLogspaceNames CovarianceSEVaryingEll
#' @export getLogspaceNames getLogspaceNames.CovarianceSEVaryingEll
#'
#' @param ... Not used.
#'
#' @usage CovarianceSEVaryingEll$logspaceNames
#'
#' @return Names of parameters to be optimized in logspace.
#'
#' @seealso \code{\link{CovarianceSEVaryingEll}}
setMethodS3("getLogspaceNames", "CovarianceSEVaryingEll", conflict="quiet",
  function(this, ...) {
    p.names <- this$getParamNamesPlain()
    return (p.names[-grep("X.ell", p.names)])
  })

#' Basenames of parameters
#'
#' Gives the "basenames" (i.e. names undecorated by the id string) of the
#' parameters.
#'
#' @name getParamNamesPlain.CovarianceSEVaryingEll
#' @aliases CovarianceSEVaryingEll$paramNamesPlain
#' @aliases getParamNamesPlain.CovarianceSEVaryingEll
#' @S3method getParamNamesPlain CovarianceSEVaryingEll
#' @export getParamNamesPlain getParamNamesPlain.CovarianceSEVaryingEll
#'
#' @param ... Not used.
#'
#' @usage CovarianceSEVaryingEll$paramNamesPlain
#'
#' @return The basenames of the parameters.
#'
#' @seealso \code{\link{CovarianceSEVaryingEll}}
setMethodS3("getParamNamesPlain", "CovarianceSEVaryingEll", conflict="quiet",
  function(this, ...) {
    p.names <- c(
      paste(sep='', 'sigma.f.', 1:length(this$.sigma.f)),
      paste(sep='', 'X.ell.', 1:length(this$.X.ell)),
      paste(sep='', 'ell.', 1:length(this$.ell)))
    return (p.names)
  })

#' Parameter values with plain names
#'
#' Gives a vector of parameter values, whose names are NOT decorated by the id
#' of this Covariance object.
#'
#' @name getParamsPlain.CovarianceSEVaryingEll
#' @aliases CovarianceSEVaryingEll$paramsPlain
#' @aliases getParamsPlain.CovarianceSEVaryingEll
#' @S3method getParamsPlain CovarianceSEVaryingEll
#' @export getParamsPlain getParamsPlain.CovarianceSEVaryingEll
#'
#' @param ... Not used.
#'
#' @usage CovarianceSEVaryingEll$paramsPlain
#'
#' @return The parameters for this covariance function, but with names
#'    undecorated by its id.
#'
#' @seealso \code{\link{setParamsPlain.Covariance}}
#' @seealso \code{\link{CovarianceSEVaryingEll}}
setMethodS3("getParamsPlain", "CovarianceSEVaryingEll", conflict="quiet",
  function(this, ...) {
    p <- c(this$.sigma.f, this$.X.ell, this$.ell)
    n <- getParamNamesPlain(this)
    names(p) <- n
    return (p)
})
setMethodS3("paramsPlainImplementation", "CovarianceSEVaryingEll", conflict="quiet",
  private=TRUE,
  function(this, p, ...) {
    warning(
      "At present, param values in 'CovarianceSEVaryingEll' are not mutable.\n")
    ## Sets any values in 'p' which match our parameter names, without worrying
    ## about lower and upper bounds (the function which calls this has the job
    ## of worrying about these!).
    #p.old <- this$getParamsPlain()
    #to.change <- names(p)[which(names(p) %in% names(p.old))]
    #p.old[to.change] <- p[to.change]
    return (invisible(this))
  })

#' Lower bounds for params, with plain names
#'
#' Gives a vector of lower bounds for the parameter values, whose names are NOT
#' decorated by the id of this Covariance object.
#'
#' @name getLowerPlain.CovarianceSEVaryingEll
#' @aliases CovarianceSEVaryingEll$lowerPlain
#' @aliases getLowerPlain.CovarianceSEVaryingEll
#' @aliases setLowerPlain.CovarianceSEVaryingEll
#' @S3method getLowerPlain CovarianceSEVaryingEll
#' @export getLowerPlain getLowerPlain.CovarianceSEVaryingEll
#'
#' @param L A (named) vector of new lower bounds (we ONLY use ones which are
#'    named, and whose names match up with names of parameters.)
#' @param ... Not used.
#'
#' @usage CovarianceSEVaryingEll$lowerPlain
#'
#' @return The lower bounds for the parameters for this covariance function, but
#'    with names undecorated by its id.
#'
#' @seealso \code{\link{getUpperPlain.CovarianceSEVaryingEll}}
#' @seealso \code{\link{CovarianceSEVaryingEll}}
setMethodS3("getLowerPlain", "CovarianceSEVaryingEll", conflict="quiet",
  function(this, ...) {
    return (this$getParamsPlain(...))
})
setMethodS3("setLowerPlain", "CovarianceSEVaryingEll", conflict="quiet",
  function(this, L, ...) {
    warning(
      "At present, lower bounds in 'CovarianceSEVaryingEll' are not mutable.\n")
    return (invisible(this))
  })

#' Upper bounds for params, with plain names
#'
#' Gives a vector of upper bounds for the parameter values, whose names are NOT
#' decorated by the id of this Covariance object.
#'
#' @name getUpperPlain.CovarianceSEVaryingEll
#' @aliases CovarianceSEVaryingEll$upperPlain
#' @aliases getUpperPlain.CovarianceSEVaryingEll
#' @aliases setUpperPlain.CovarianceSEVaryingEll
#' @S3method getUpperPlain CovarianceSEVaryingEll
#' @export getUpperPlain getUpperPlain.CovarianceSEVaryingEll
#'
#' @param U A (named) vector of new upper bounds (we ONLY use ones which are
#'    named, and whose names match up with names of parameters.)
#' @param ... Not used.
#'
#' @usage CovarianceSEVaryingEll$upperPlain
#'
#' @return The upper bounds for the parameters for this covariance function, but
#'    with names undecorated by its id.
#'
#' @seealso \code{\link{getLowerPlain.CovarianceSEVaryingEll}}
#' @seealso \code{\link{CovarianceSEVaryingEll}}
setMethodS3("getUpperPlain", "CovarianceSEVaryingEll", conflict="quiet",
  function(this, ...) {
    return (this$getParamsPlain(...))
  })
setMethodS3("setUpperPlain", "CovarianceSEVaryingEll", conflict="quiet",
  function(this, U, ...) {
    warning(
      "At present, lower bounds in 'CovarianceSEVaryingEll' are not mutable.\n")
    return (invisible(this))
  })

#' ell(X)
#'
#' Evaluate ell(X) at the applied X-points (by spline-interpolating)
#'
#' @S3method ell CovarianceSEVaryingEll
#' @export ell ell.CovarianceSEVaryingEll
#' @name ell.CovarianceSEVaryingEll
#'
#' @param X The X-values where we want to know \code{ell(X)}.
#' @param ... Not used.
#'
#' @note The spline interpolation takes place in log-space (\code{ell} can never
#'    be negative), but a real-space result is returned.
#'
#' @return \code{ell} evaluated at each point in \code{X}.
#'
#' @seealso \code{\link{CovarianceSEVaryingEll}}
setMethodS3("ell", "CovarianceSEVaryingEll", conflict="quiet",
  function(this, X, ...) {
    # Evaluate ell(X) at the applied X-points (by spline-interpolating)
    return (exp(spline(xout=X, x=this$.X.ell, y=log(this$.ell))$y))
    return (spline(xout=X, x=this$.X.ell, y=this$.ell)$y)
  })

#' sigma.f(X)
#'
#' Evaluate sigma.f(X) at the applied X-points (by spline-interpolating)
#'
#' @S3method sigma.f CovarianceSEVaryingEll
#' @export sigma.f sigma.f.CovarianceSEVaryingEll
#' @name sigma.f.CovarianceSEVaryingEll
#'
#' @param X The X-values where we want to know \code{sigma.f(X)}.
#' @param ... Not used.
#'
#' @note The spline interpolation takes place in log-space (\code{sigma.f} can
#'    never be negative), but a real-space result is returned.
#'
#' @return \code{sigma.f} evaluated at each point in \code{X}.
#'
#' @seealso \code{\link{CovarianceSEVaryingEll}}
setMethodS3("sigma.f", "CovarianceSEVaryingEll", conflict="quiet",
  function(this, X, ...) {
    # This code repeats a single-valued sigma if necessary, so it can handle
    # both single-sigma and varying-sigma.
    sigma.vals <- data.frame(x=this$.X.ell, y=this$.sigma.f)
    # Evaluate sigma.f(X) at the applied X-points (by spline-interpolating)
    sigma.at.X <- with(sigma.vals, spline(x=x, y=y, xout=X)$y)
    return (sigma.at.X)
  })

#' Nonstationary Squared-exponential Covariance matrix
#'
#' Calculates a covariance matrix for the nonstationary squared-exponential
#' covariance function.
#'
#' @S3method K.specific CovarianceSEVaryingEll
#' @export K.specific K.specific.CovarianceSEVaryingEll
#' @name K.specific.CovarianceSEVaryingEll
#'
#' @param X  X-values for the input points (i.e., where we have data)
#' @param X.out  X-values for the points desired to predict
#' @param ... Not used.
#'
#' @return The covariance matrix taking \code{X} into \code{X.out}, based on the
#'    parameter values in \code{this}.
#'
#' @seealso \code{\link{CovarianceSEVaryingEll}}
setMethodS3("K.specific", "CovarianceSEVaryingEll", conflict="quiet",
  function(this, X, X.out=X, ...) {
    X.dist <- DistanceMatrix(X=X, X.out=X.out)
    # Check and warn if user wants to extrapolate
    if (min(X.out) < min(X) || max(X.out) > max(X)) {
      warning(
        "Extrapolation is highly inadvisable with CovarianceSEVaryingEll!\n")
    }
    # Get params and calculate the matrix
    ell     <- this$ell(X=X    )
    ell.out <- this$ell(X=X.out)
    sigma     <- this$sigma.f(X=X    )
    sigma.out <- this$sigma.f(X=X.out)
    sum.ell.sq <- outer(ell.out ^ 2, ell ^ 2, '+')
    K <- (outer(sigma.out, sigma) 
      * sqrt(2.0 * outer(ell.out, ell) / sum.ell.sq)
      * exp(-(X.dist ^ 2) / sum.ell.sq))
    gc()
    return (K)
  })

#' Element-wise derivatives of Covariance matrix
#'
#' Calculate the element-wise derivative of \code{KInIn}, with respect to the
#' parameter whose (plain) name is \code{param}.
#'
#' @S3method KDerivImplementation CovarianceSEVaryingEll
#' @export KDerivImplementation KDerivImplementation.CovarianceSEVaryingEll
#' @name KDerivImplementation.CovarianceSEVaryingEll
#'
#' @param d  The Dataset whose X-values determine KInIn.
#' @param param  The (plain) name of the parameter with respect to which we're
#'    differentiating.
#' @param ... Not used.
#'
#' @return A matrix whose elements are the derivatives of the corresponding
#'    elements in KInIn, with respect to the parameter \code{param}.
#'
#' @seealso \code{\link{CovarianceSEVaryingEll}}
setMethodS3("KDerivImplementation", "CovarianceSEVaryingEll", conflict="quiet",
  function(this, d, param, ...) {
    warning("This function not yet implemented for CovarianceSEVaryingEll.\n")
    return (NULL)
  })

#' Nonstationary SE variance at each point
#'
#' Calculate the nonstationary SE variance of the points at X: i.e., the a
#' priori uncertainty at each point.
#'
#' @S3method Variance CovarianceSEVaryingEll
#' @export Variance Variance.CovarianceSEVaryingEll
#' @name Variance.CovarianceSEVaryingEll
#'
#' @param X  The points we want to know the nonstationary SE variance at.
#' @param ... Not used.
#'
#' @return A numeric vector of the same length as X, with the corresponding
#'    nonstationary SE variance.
#'
#' @seealso \code{\link{CovarianceSEVaryingEll}}
setMethodS3("Variance", "CovarianceSEVaryingEll", conflict="quiet",
  function(this, X, ...) {
    return (rep(mean(this$.sigma.f) ^ 2, NumPoints(X)))
  })

