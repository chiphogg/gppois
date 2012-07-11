#############################################################################/**
# @RdocClass Model
#
# @title "A trainable collection of Covariances"
#
# \description{
#   This is the \emph{constructor} for a \code{Model} object:
#   @get "title".
#
#   The \code{Model} expresses our beliefs or knowledge about a Dataset.  We
#   assume a Dataset can be modeled as the sum of one or more Gaussian Process
#   Covariance functions, each of which is governed by parameters.  The Model
#   can then be trained on a Dataset, a process which selects the parameter
#   values that best describe the data.  The model can then be used to make
#   predictions about the true function -- either at noisy datapoints, or
#   interpolating into data-free regions, or both.
#
#   This class should work just fine, as long as
#     a) we are training on all the datapoints together (i.e., not breaking
#        them up into subregions to divide-and-conquer), and
#     b) this$params returns a vector which is amenable to simple optimization
#        routines (i.e., none of the parameters require special treatment).
#   If either of these conditions fail, a new approach is needed: either a
#   specialized subclass should be created, or the problem should be broken
#   into smaller pieces where these assumptions are good.
#
#   Ironically, \emph{both} these conditions fail for \emph{both} scenarios
#   considered in our Journal of Applied Crystallography paper, despite the
#   fact that I wrote this software to perform the analysis for that paper.
#   I hope to remedy this in a future version.  However, even in the meantime,
#   having these classes still makes it very much easier to build the
#   specialized functions I need.  Moreover, experience has shown that the
#   plain-vanilla class structure is already good enough for other applications
#   unrelated to denoising of scattering curves.
#
#   Here is the class hierarchy:\cr
#   @classhierarchy
#
# }
#
# @synopsis
#
# \arguments{
#   \item{id}{(character) An id which identifies this Model.}
#   \item{...}{Not used.}
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
setConstructorS3("Model",
  function(id="", ...) {
    # Constructs a Model object with the given ID.
    #
    # Args:
    #   id:  A string which labels this model.
    #
    # Returns:
    #   A Model with the given ID.

    extend(Object(), "Model",
      .id             = id,
      .K.chol         = LazyMatrix(),
      .L              = LazyMatrix(),
      .last.trained.d = NA,
      .contributions  = list()
      )
  })

#' Log of the Marginal Likelihood
#'
#' Sets the parameter values to \code{par} for Model \code{model}, and returns
#' the log of the (M)arginal (L)ikelihood for describing Dataset \code{d}.
#' 
#' @param par  The parameter values to use for calculating.
#' @param model  The Model object we're optimizing.
#' @param d  The Dataset we're training on.
#' @param update.params  logical; if TRUE, we should change model's params to the
#'      values in par.
#'
#' @references Rasmussen, C.E. and C.K.I. Williams (2005.) Gaussian Processes
#'   for Machine Learning. The MIT Press.
#'   \url{http://www.gaussianprocess.org/gpml/}
#'
#' @export
#'
#' @return The logarithm of the marginal likelihood; also, has the side-effect
#'    that \code{moedel$params} are set to \code{par} unless
#'    \code{update.params=FALSE}.
LogML <- function(par=model$getParams(for.training=TRUE), model, d,
  update.params=TRUE) {
  if (!update.params) {
    model <- clone(model)
  }
  model$setParams(p=DecodeForTraining(par))
  Y <- d$xformedDpts
  # The following calculation is based on Equation 5.8 in
  # Rasmussen and Williams (2005):
  term.data.fit   <- -0.5 * t(Y) %*% model$KInv(d) %*% Y
  term.complexity <- -0.5 * model$LogDetK(d)
  term.num.dpts   <- -0.5 * d$n * log(2 * pi)
  return (term.data.fit + term.complexity + term.num.dpts)
}

LogML1D <- function(value, name, model, d, update.params=TRUE) {
  # Calculate the log of the marginal likelihood for Model 'model' on Dataset
  # 'd', if model parameter given by 'name' is replaced by 'value'.
  #
  # (This is useful for optimizing when only one parameter actually changes.
  # End users don't need to know about this, so it isn't publicly documented.)
  #
  # Args:
  #   value:  The new value of the parameter
  #   name:  The name of the parameter
  #   model:  The Model object we're optimizing.
  #   d:  The Dataset we're training on.
  #   update.params:  logical; if TRUE, we should change model's params to the
  #      given value.
  #
  # Returns:
  #   The log of the marginal likelihood.
  names(value) <- name
  return (LogML(par=value, model=model, d=d, update.params=update.params))
}

#' Gradient of log marginal likelihood
#'
#' Sets the parameter values to \code{par} for Model \code{model}, and returns
#' the gradient (w.r.t. par) of the log of the (M)arginal (L)ikelihood for
#' describing Dataset \code{d}.  The gradient of LogML is very helpful in
#' optimizing the parameter values.
#'
#' @export
#'
#' @param par  The parameter values to use for calculating.
#' @param model  The Model object we're optimizing.
#' @param d  The Dataset we're training on.
#' @param update.params  logical; if TRUE, we should change model's params to
#'    the given value.
#'
#' @references Rasmussen, C.E. and C.K.I. Williams (2005.) Gaussian Processes
#'   for Machine Learning. The MIT Press.
#'   \url{http://www.gaussianprocess.org/gpml/}
#'
#' @return The gradient of the log of the marginal likelihood (also has a
#'    side-effect of setting \code{model$params <- par}).
GradLogML <- function(par=model$getParams(for.training=TRUE), model, d,
  update.params=TRUE) {
  if (!update.params) {
    model <- clone(model)
  }
  model$setParams(p=DecodeForTraining(par))
  Y <- d$xformedDpts
  # The following calculations are based on Equation 5.9 in
  # Rasmussen and Williams (2005).
  alpha <- model$KInv(d) %*% Y
  mat.for.grad <- alpha %*% t(alpha) - model$KInv(d)
  var.names <- names(model$getParams(for.training=TRUE))
  good.names <- names(par)[which(names(par) %in% var.names)]
  grad <- c()
  for (p.n in good.names) {
    grad[p.n] <- 0.5 * SmartTrace(model$KDeriv(d=d, param=p.n), mat.for.grad)
  }
  return (grad)
}

#' Id's of each contributing Covariance
#'
#' Recall that a Model consists of a collection of Covariance objects.  This
#' function gives their ID's.
#'
#' @name getContributionIds.Model
#' @aliases Model$contributionIds getContributionIds.Model
#' @S3method getContributionIds Model
#' @export getContributionIds getContributionIds.Model
#'
#' @param this The Model object.
#' @param ... Not used.
#'
#' @usage Model$contributionIds
#'
#' @return The id's for this Model's contributing Covariance objects.
#'
#' @seealso \code{\link{Model}}
setMethodS3("getContributionIds", "Model", conflict="quiet",
  function(this, ...) {
    id.list <- c()
    for (covar in this$.contributions) {
      id.list <- c(id.list, covar$id)
    }
    return (id.list)
  })

#' ID string for this Model
#'
#' A character string identifying this Model object.
#'
#' @name getId.Model
#' @aliases Model$id getId.Model setId.Model
#' @S3method getId Model
#' @export getId getId.Model
#'
#' @param this The \code{Model} whose contributions to list.
#' @param id (character) The new ID for \code{this}.
#' @param ... Not used.
#'
#' @usage Model$id
#' @usage Model$id <- id
#'
#' @return The id of \code{this}.
#'
#' @seealso \code{\link{Model}}
setMethodS3("getId", "Model", conflict="quiet",
  function(this, ...) {
    return (this$.id)
  })
setMethodS3("setId", "Model", conflict="quiet",
  function(this, id, ...) {
    this$.id <- id
    return (this)
  })

#' Parameters for the Model
#'
#' A named vector of parameters governing this Model.  Note that the names are
#' \emph{decorated} by prepending the Covariance id; this prevents namespace
#' collisions.
#'
#' @name getParams.Model
#' @aliases Model$params getParams.Model setParams.Model
#' @S3method getParams Model
#' @export getParams getParams.Model
#' @S3method setParams Model
#' @export setParams setParams.Model
#'
#' @param this The Model object.
#' @param p A (named) vector of new parameter values (we ONLY use ones which are
#'      named, and whose names match up with names of parameters.)
#' @param for.training  If TRUE, we ignore "constant" parameters (i.e., where
#'      lower=upper) and return the *log* of any "scale" parameters (such as ell
#'      or sigma.f for the SE model).
#' @param ... Not used.
#'
#' @usage Model$params
#' @usage Model$params <- c(name1=value1, name2=value2, etc.)
#'
#' @return A vector with the current values for each parameter.
#'
#' @seealso \code{\link{getParams.Covariance}} for more about \code{for.training}
#' @seealso \code{\link{getUpper.Model}}
#' @seealso \code{\link{getLower.Model}}
#' @seealso \code{\link{Model}}
setMethodS3("getParams", "Model", conflict="quiet",
  function(this, for.training=FALSE, ...) {
    p <- c()
    for (covar in this$.contributions) {
      p <- c(p, covar$getParams(for.training=for.training))
    }
    if (for.training) {
      unlog.params <- DecodeForTraining(p)
      i.vary <- which(names(unlog.params) %in% this$getVaryingParamNames())
      p <- p[i.vary]
    }
    return (p)
  })
setMethodS3("setParams", "Model", conflict="quiet",
  function(this, p, for.training=FALSE, ...) {
    for (covar in this$.contributions) {
      covar$setParams(p=p, for.training=for.training)
    }
    return (invisible(this))
  })

#' Lower bounds for parameters
#'
#' Lower bounds for the parameter values.
#'
#' @name getLower.Model
#' @aliases Model$lower getLower.Model setLower.Model
#' @S3method getLower Model
#' @export getLower getLower.Model
#' @S3method setLower Model
#' @export setLower setLower.Model
#'
#' @param this The Model object.
#' @param L A (named) vector of new parameter values (we ONLY use ones which are
#'      named, and whose names match up with names of parameters.)
#' @param for.training  If TRUE, we ignore "constant" parameters (i.e., where
#'      lower=upper) and return the *log* of any "scale" parameters (such as ell
#'      or sigma.f for the SE model).
#' @param ... Not used.
#'
#' @usage Model$lower
#' @usage Model$lower <- c(name1=value1, name2=value2, ...)
#'
#' @return The lower bounds for the parameters for this model.
#'
#' @seealso \code{\link{getUpper.Model}}
#' @seealso \code{\link{getParams.Model}}
#' @seealso \code{\link{Model}}
setMethodS3("getLower", "Model", conflict="quiet",
  function(this, for.training=FALSE, ...) {
    L <- c()
    for (covar in this$.contributions) {
      L <- c(L, covar$getLower(for.training=for.training))
    }
    if (for.training) {
      unlog.params <- DecodeForTraining(L)
      i.vary <- which(names(unlog.params) %in% this$getVaryingParamNames())
      L <- L[i.vary]
    }
    return (L)
  })
setMethodS3("setLower", "Model", conflict="quiet",
  function(this, L, for.training=FALSE, ...) {
    for (covar in this$.contributions) {
      covar$setLower(L=L, for.training=for.training)
    }
    return (invisible(this))
  })

#' Upper bounds for parameters
#'
#' Upper bounds for the parameter values.
#'
#' @name getUpper.Model
#' @aliases Model$upper getUpper.Model setUpper.Model
#' @S3method getUpper Model
#' @export getUpper getUpper.Model
#' @S3method setUpper Model
#' @export setUpper setUpper.Model
#'
#' @param this The Model object.
#' @param L A (named) vector of new parameter values (we ONLY use ones which are
#'      named, and whose names match up with names of parameters.)
#' @param for.training  If TRUE, we ignore "constant" parameters (i.e., where
#'      lower=upper) and return the *log* of any "scale" parameters (such as ell
#'      or sigma.f for the SE model).
#' @param ... Not used.
#'
#' @usage Model$upper
#' @usage Model$upper <- c(name1=value1, name2=value2, ...)
#'
#' @return The upper bounds for the parameters for this model.
#'
#' @seealso \code{\link{getLower.Model}}
#' @seealso \code{\link{getParams.Model}}
#' @seealso \code{\link{Model}}
setMethodS3("getUpper", "Model", conflict="quiet",
  function(this, for.training=FALSE, ...) {
    U <- c()
    for (covar in this$.contributions) {
      U <- c(U, covar$getUpper(for.training=for.training))
    }
    if (for.training) {
      unlog.params <- DecodeForTraining(U)
      i.vary <- which(names(unlog.params) %in% this$getVaryingParamNames())
      U <- U[i.vary]
    }
    return (U)
  })
setMethodS3("setUpper", "Model", conflict="quiet",
  function(this, U, for.training=FALSE, ...) {
    for (covar in this$.contributions) {
      covar$setUpper(U=U, for.training=for.training)
    }
    return (invisible(this))
  })

#' ID's of non-noise Covariances
#'
#' A \code{Model} is a collection of \code{Covariances}.  Some subset of these
#' will correspond to the signal of actual interest; the rest are considered
#' noise.  This function returns the list of contributions which are considered
#' to be \dQuote{signal}.
#'
#' @name getSignalIds.Model
#' @aliases Model$signalIds getSignalIds.Model
#' @S3method getSignalIds Model
#' @export getSignalIds getSignalIds.Model
#'
#' @param this The Model object.
#' @param ... Not used.
#'
#' @usage Model$signalIds
#'
#' @return The id's for this Model's non-noise contributing Covariance objects.
#'
#' @seealso \code{\link{Model}}
setMethodS3("getSignalIds", "Model", conflict="quiet",
  function(this, ...) {
    return (this$getContributionIds()[
        which(this$getContributionIds() != 'noise')])
  })

#' Non-constant Model parameters
#'
#' Names of the parameters which are not constant.
#'
#' @name getVaryingParamNames.Model
#' @aliases Model$varyingParamNames getVaryingParamNames.Model
#' @S3method getVaryingParamNames Model
#' @export getVaryingParamNames getVaryingParamNames.Model
#'
#' @param this The Model object.
#' @param ... Not used.
#'
#' @usage Model$varyingParamNames
#'
#' @return The names of all parameters which are not constant.
#'
#' @seealso \code{\link{Model}}
setMethodS3("getVaryingParamNames", "Model", conflict="quiet",
  function(this, ...) {
    U <- this$getUpper()
    L <- this$getLower()
    p.names <- names(U)
    i.const <- which(U[p.names] == L[p.names])
    if (length(i.const) > 0) {
      return (p.names[-i.const])
    } else {
      return (p.names)
    }
  })

#' Add a new Covariance to this Model
#'
#' Add another Covariance structure to this model (it contributes additively).
#' We CLONE it so we OWN it.  (Don't want anyone else to fiddle with the
#' Covariance object, EXCEPT this Model object.)
#'
#' @S3method AddCovariance Model
#' @export AddCovariance AddCovariance.Model
#' @name AddCovariance.Model
#'
#' @param covariance  A Covariance object to be cloned and added to this model.
#' @param on.duplicate.id  (character) One of (\dQuote{rename},
#'    \dQuote{replace}), directing how to handle a new contribution whose ID is
#'    the same as an existing one's.
#' @param ... Not used.
#'
#' @seealso \code{\link{Model}}
#' @seealso \code{\link{Covariance}}
setMethodS3("AddCovariance", "Model", conflict="quiet",
  function(this, covariance, on.duplicate.id="rename", ...) {
    our.covar <- clone(covariance)
    # Make sure this contribution has a unique id:
    if (on.duplicate.id == "rename") {
      label <- 1
      unique.id <- our.covar$id
      old.ids <- this$contributionIds
      while (unique.id %in% old.ids) {
        unique.id <- paste(sep='', our.covar$id, '.', label)
        label <- label + 1
      }
      our.covar$id <- unique.id
    } else if (on.duplicate.id == "replace") {
      this$.contributions[[our.covar$id]] <- NULL
    }

    # We want to be able to refer to this contribution by its id.
    new.contribution <- list(our.covar)
    names(new.contribution) <- our.covar$id
    this$.contributions <- c(this$.contributions, new.contribution)
    return (invisible(this))
  })

#' Deep-clone a Model
#'
#' Deep-clone a Model (i.e., clones Covariance objects too).
#'
#' @method clone Model
#'
#' @param this The Model object to clone.
#' @param ... Not used.
#'
#' @return A deep clone of the Model object.
#'
#' @export
#' @seealso \code{\link{Model}}
clone.Model <- function(this, ...) {
  M <- clone.Object(this)
  M$.contributions <- list()
  for (covar in this$.contributions) {
    covar.clone <- list(clone(covar))
    names(covar.clone) <- covar$id
    M$.contributions <- c(M$.contributions, covar.clone)
  }
  M$.K.chol <- clone(M$.K.chol)
  M$.L <- clone(M$.L)
  return (M)
}

#' Clear precomputed matrices from memory
#'
#' Model objects store computed matrices for easy reuse.  This is usually very
#' convenient, but sometimes we want the memory more than we want the speed!
#' This method tells the Model to forget its results and free up some memory.
#'
#' @S3method Forget Model
#' @export Forget Forget.Model
#' @name Forget.Model
#'
#' @param this The Model object.
#' @param ... Not used.
#'
#' @seealso \code{\link{Model}}
#' @seealso \code{\link{LazyMatrix}}
setMethodS3("Forget", "Model", conflict="quiet",
  function(this, ...) {
    this$.K.chol <- LazyMatrix()
    this$.L <- LazyMatrix()
    this$.last.trained.d <- NA
    gc()
  })

#' Make some parameters constant
#'
#' Makes a subset of parameters constant, by setting the upper and lower bounds
#' equal to the current parameter value.  Defaults to freezing all parameters.
#'
#' @S3method Freeze Model
#' @export Freeze Freeze.Model
#' @name Freeze.Model
#'
#' @param this The Model object.
#' @param p.names  The names of the parameters to freeze.
#' @param ... Not used.
#'
#' @seealso \code{\link{Model}}
setMethodS3("Freeze", "Model", conflict="quiet",
  function(this, p.names=names(this$params), ...) {
    good.i <- which(p.names %in% names(this$params))
    for (p.name in p.names[good.i]) {
      this$lower <- this$params[p.name]
      this$upper <- this$params[p.name]
    }
  })

#' Lower Cholesky root of covariance matrix
#'
#' Compute the lower-triangular Cholesky decomposition of the covariance matrix.
#' This is useful for taking random draws from the posterior.
#'
#' @S3method L Model
#' @export L L.Model
#' @name L.Model
#'
#' @param this The Model object.
#' @param d  The Dataset we're training on.
#' @param X.out  A matrix (with d$d columns) of X-locations where we want
#'    predictions.
#' @param contributions The names of the parameters which are considered
#'    \dQuote{signal}.
#' @param ... Not used.
#'
#' @return The lower-triangular Cholesky decomposition of the covariance matrix
#'    K.  (i.e., the matrix L such that L %*% t(L) = K.)
#'
#' @seealso \code{\link{Model}}
setMethodS3("L", "Model", conflict="quiet",
  function(this, d, X.out=d$X, contributions=this$getSignalIds(), ...) {
    this$ComputeL(d=d, X.out=X.out, contributions=contributions)
    return (this$.L$M)
  })

#' Retrieve one contributing Covariance
#'
#' Retrieve a clone of the first contributing Covariance objects with the given
#' id.
#'
#' @S3method NamedCovariance Model
#' @export NamedCovariance NamedCovariance.Model
#' @name NamedCovariance.Model
#'
#' @param this The Model object.
#' @param id (character) The ID of the Covariance to retrieve.
#' @param ... Not used.
#'
#' @return A clone of the first contributing Covariance objects with the given
#' id.
#'
#' @seealso \code{\link{Model}}
#' @seealso \code{\link{Covariance}}
setMethodS3("NamedCovariance", "Model", conflict="quiet",
  function(this, id, ...) {
    for (covar in this$.contributions) {
      if (covar$id == id) {
        Cov <- clone(covar)
        return (Cov)
      }
    }
    return (NA)
  })

#' Animated uncertainty in a surface
#'
#' This function helps visualize uncertainty using animations.  (See Details.)
#'
#' @S3method PlotBubblingSurfaces2D Model
#' @export PlotBubblingSurfaces2D PlotBubblingSurfaces2D.Model
#' @name PlotBubblingSurfaces2D.Model
#'
#' @param this The Model object.
#' @param d  The Dataset to evaluate the Model on.
#' @param X.out  (matrix) The X-points where we want to predict the function.
#' @param contributions  (character vector) Id's of the contributing Covariances
#'     we want to predict (other contributions are considered noise).
#' @param n.indep  The number of independent draws to take for each datapoint.
#' @param n.times  The final number of animation frames.
#' @param file.name  (character) The basename of the file to plot to.
#' @param ... Not used.
#'
#' @details
#' So, we have a 2D dataset, and a posterior distribution on the underlying true
#' function.  We want to visualize the \emph{uncertainty} in that true function.
#' How?
#'
#' One way to visualize a posterior distribution is to take a number of draws
#' from it.  This is tricky for surfaces; they will tend to overlap and
#' intersect and the plot will be too cluttered.  So instead of plotting the
#' random draws simultaneously, why not show them one at a time?  This does fix
#' the clutter problem, but consecutive draws are unrelated to each other, so
#' now we have a "jumpiness" problem.
#'
#' Both problems can be solved simultaneously, if the draws satisfy these
#' conditions:
#' \enumerate{
#'   \item Consecutive draws are \emph{not} independent; they are correlated,
#'      and the correlation approaches 1 as the time difference approaches 0.
#'   \item However, every draw does have the same \emph{marginal} distribution,
#'      which is equal to the posterior distribution we're trying to visualize.
#' }
#' The first condition makes the motion continuous.  The second makes sure it
#' represents the uncertainty we're actually trying to represent.
#'
#' This software implements a novel solution which satisfies a stronger version
#' of the first condition: specifically, the path of every point on the surface
#' is differentiable infinitely many times.  Thus, we have the smoothest
#' possible animations which actually represent the uncertainty we're trying to
#' show.  I came up with this solution in May 2012; I am currently planning to
#' write it up into a paper after I get back from Japan at the end of June.
#'
#'
#' @seealso \code{\link{BubblingRandomMatrix}}
#' @seealso \code{\link{Model}}
setMethodS3("PlotBubblingSurfaces2D", "Model", conflict="quiet",
  function(this, d, X.out=d$X, contributions=this$getSignalIds(),
    n.indep=10, n.times=50, file.name=NULL, ...) {
    B <- BubblingRandomMatrix(n.pts=NumPoints(X.out), N=n.indep,
      n.times=n.times, ...)
    L <- this$L(d=d, X.out=X.out, contributions=contributions)
    Y <- matrix(nrow=NumPoints(X.out), rep(times=n.times,
        this$PosteriorMean(d=d, X.out=X.out, contributions=contributions))) + (
        L %*% B)
    unit <- min(dist(d$X))
    for (i in 1:ncol(Y)) {
      rgl.clear()
      PlotSurface(X=X.out, Y=Y[, i], ...)
      rgl.spheres(x=d$X[, 1], z=d$X[, 2], y=d$dpts, radius=0.2 * unit)
      rgl.snapshot(filename=sprintf("%s_t-%04d.png", file.name, i), top=TRUE)
    }
  })

#' Best estimate of the true function
#'
#' Computes this Model's optimal prediction of the underlying function's value
#' at every point in 'X.out'.
#'
#' @S3method PosteriorMean Model
#' @export PosteriorMean PosteriorMean.Model
#' @name PosteriorMean.Model
#'
#' @param this The Model object.
#' @param d  The Dataset to train the Model on.
#' @param X.out  (matrix) The X-points where we want to predict the function.
#' @param contributions  (character vector) Id's of the contributing Covariances
#'    we want to predict (other contributions are considered noise); default is
#'    every contribution not named 'noise'.
#' @param untransform.result  logical; if TRUE, we transform back to the space
#'    of \code{\link{dpts}} (as opposed to the space of
#'    \code{\link{xformedDpts}} where training takes place).
#' @param ... Not used.
#'
#' @return A numeric vector with optimal predictions at every point in X.out.
#'
#' @seealso \code{\link{Model}}
setMethodS3("PosteriorMean", "Model", conflict="quiet",
  function(this, d, X.out=d$X, contributions=this$getSignalIds(),
    untransform.result=TRUE, ...) {
    contributions <- this$CheckContributionsAndWarn(contributions)
    M <- this$PredictionMatrix(d=d, X.out=X.out, contributions=contributions)
    result <- M %*% d$xformedDpts
    if (untransform.result) {
      result <- d$Untransform(result)
    }
    return (result)
  })

#' Matrix connecting noisy data to true function
#'
#' A matrix relating function values at the output points to function values at
#' the input points.
#'
#' @S3method PredictionMatrix Model
#' @export PredictionMatrix PredictionMatrix.Model
#' @name PredictionMatrix.Model
#'
#' @param this The Model object.
#' @param d  The Dataset to train the Model on.
#' @param X.out  matrix; the X-points where we want to predict the function.
#' @param contributions  character vector; id's of the contributing Covariances
#'      we want to predict (other contributions are considered noise).
#' @param ... Not used.
#'
#' @return The matrix \code{M} which yields the predictions at \code{X.out}
#'    when multiplied by the datapoint vector.
#'
#' @seealso \code{\link{Model}}
setMethodS3("PredictionMatrix", "Model", conflict="quiet",
  function(this, d, X.out=d$X, contributions=this$getSignalIds(), ...) {
    # A matrix relating function values at the output points to function values
    # at the input points.
    #
    #
    K.in.out <- matrix(0, nrow=NumPoints(X.out), ncol=d$n)
    for (c.name in contributions) {
      covar <- this$.contributions[[c.name]]
      covar.K <- covar$KInOut(d=d, X.out=X.out)
      K.in.out <- K.in.out + covar.K
      rm(covar.K)
      gc()
    }
    M <- K.in.out %*% this$KInv(d=d)
    rm(K.in.out)
    gc()
    return (M)
  })

#' Best estimate, including uncertainty
#'
#' Computes bounds on the uncertainty (in terms of the standard deviation) in
#' the prediction at a given point, along with the prediction.
#'
#' @S3method PosteriorInterval Model
#' @export PosteriorInterval PosteriorInterval.Model
#' @name PosteriorInterval.Model
#'
#' @param this The Model object.
#' @param d  The Dataset to train the Model on.
#' @param X.out  (matrix) The X-points where we want to predict the function.
#' @param num.sd  The number of standard deviations from the mean our interval
#'     should include (defaults to 1).
#' @param contributions  (character vector) Id's of the contributing Covariances
#'      we want to predict (other contributions are considered noise);
#'      default is every contribution not named 'noise'.
#' @param ... Not used.
#'
#' @return A numeric vector with optimal predictions at every point in X.out.
#'
#' @seealso \code{\link{PredictionMatrix.Model}}
#' @seealso \code{\link{PosteriorMean.Model}}
#' @seealso \code{\link{Model}}
setMethodS3("PosteriorInterval", "Model", conflict="quiet",
  function(this, d, X.out=d$X, num.sd=1, contributions=this$getSignalIds(),
    ...) {
    M <- this$PredictionMatrix(d=d, X.out=X.out, contributions=contributions)
    std.dev <- this$PosteriorStandardDeviation(d=d, X.out=X.out,
      contributions=contributions, ...)
    prediction <- M %*% d$xformedDpts
    return (data.frame(X=X.out,
        mean=d$Untransform(prediction),
        lower=d$Untransform(prediction - num.sd * std.dev),
        upper=d$Untransform(prediction + num.sd * std.dev)))
  })

#' Pointwise uncertainty
#'
#' Computes the posterior \dQuote{sigma} at a given point.
#'
#' @S3method PosteriorStandardDeviation Model
#' @export PosteriorStandardDeviation PosteriorStandardDeviation.Model
#' @name PosteriorStandardDeviation.Model
#'
#' @param this The Model object.
#' @param d  The Dataset to train the Model on.
#' @param X.out  (matrix) the X-points where we want to predict the function.
#' @param contributions  (character vector) Id's of the contributing Covariances
#'      we want to predict (other contributions are considered noise);
#'      default is every contribution not named 'noise'.
#' @param ... Not used.
#'
#' @return A numeric vector with the posterior standard deviation at every point
#'    in X.out.
#'
#' @seealso \code{\link{PredictionMatrix.Model}}
#' @seealso \code{\link{PosteriorMean.Model}}
#' @seealso \code{\link{Model}}
setMethodS3("PosteriorStandardDeviation", "Model", conflict="quiet",
  function(this, d, X.out=d$X, contributions=this$getSignalIds(), ...) {
    contributions <- this$CheckContributionsAndWarn(contributions)
    # Calculate the posterior predictive mean.
    N.out <- NumPoints(X.out)
    K.in.out <- matrix(0, nrow=N.out, ncol=d$n)
    variance <- rep(0, N.out)
    for (c.name in contributions) {
      covar <- this$.contributions[[c.name]]
      K.in.out <- K.in.out + covar$KInOut(d=d, X.out=X.out)
      variance <- variance + covar$Variance(X=X.out)
      gc()
    }
    M <- this$PredictionMatrix(d=d, X.out=X.out, contributions=contributions)
    # The 'rowSums' bit is a fancy way to calculate only the diagonal elements
    # of the matrix product (K.in.out %*% this$KInv(d=d) %*% t(K.in.out)).
    # We force the variance to be non-negative; this is necessary because
    # roundoff errors in the matrix multiplication can cause it to go negative
    # by negligible amounts.
    variance <- pmax(0, variance - rowSums(M * K.in.out))
    std.dev <- sqrt(variance)
    return (std.dev)
  })

#' Pretty-printing for Model objects
#'
#' Prints out the id of the Model, followed by a list of its contributing
#' Covariances and their IDs.  Then, prints each parameter, including the name,
#' lower bound, current value, and upper bound.
#'
#' @method print Model
#'
#' @param this The Model object to print.
#' @param indent Aids in formatting: the number of spaces to print before every
#'    line.
#' @param ... Not used.
#'
#' @export
#' @seealso \code{\link{Model}}
print.Model <- function(this, indent=0, ...) {
  tab <- Spaces(num=indent)
  cat(sprintf("%s%s, id='%s'\n", tab, class(this)[1], this$id))
  cat(sprintf("%s%sCONTRIBUTING COVARIANCES:\n", tab, Spaces(2)))
  for (covar in this$.contributions) {
    cat(sprintf("%s%sid=%-20s (%s)\n", tab, Spaces(4), Wrap(covar$id, "'"),
        class(covar)[1]))
  }
  PrintParams(lower=this$lower, upper=this$upper, params=this$params,
    indent=indent)
  return (invisible(this))
}

#' Uncertainty about the noise level
#'
#' The noise level is not necessarily known precisely a priori.  If it is
#' uncertain, this function lets you set a range of possible values.
#'
#' Keep in mind that the noise is a \dQuote{scale-type} parameter, so we are
#' uncertain about its \emph{order of magnitude}.  Thus, it's a good idea to
#' pass boundaries like, say, \code{c(1e-7, 1e-3)}, rather than boundaries which
#' have the same order of magnitude.  (Unless of course you do know the order of
#' magnitude ahead of time!)
#'
#' @S3method SetNoiseBounds Model
#' @export SetNoiseBounds SetNoiseBounds.Model
#' @name SetNoiseBounds.Model
#'
#' @param this The Model object.
#' @param sigma.vals  A numeric vector, such that range(sigma.vals) sets the
#'    range of values for the noise level.
#' @param ... Not used.
#'
#' @seealso \code{\link{Model}}
setMethodS3("SetNoiseBounds", "Model", conflict="quiet",
  function(this, sigma.vals, ...) {
    this$AddCovariance(CovarianceNoise(id="noise", sigma.bounds=sigma.vals),
      on.duplicate.id="replace")
  })

#' Train a Model on a Dataset
#'
#' Optimize this Model's parameters so they describe the given data.
#'
#' @S3method Train Model
#' @export Train Train.Model
#' @name Train.Model
#'
#' @param this The Model object.
#' @param d  (Dataset) The data which our parameters should describe.
#' @param ... Not used.
#'
#' @details Presently, we use R's \code{\link{optim}} \dQuote{under the hood}.
#'    Specifically, we use the L-BFGS-B approach because we found it gave good
#'    results, and because it naturally incorporates boundaries.
#'
#' Future versions may have an option to calculate the Hessian analytically,
#' which would allow trust region approaches to be used.  This should speed
#' convergence.
#'
#' @seealso \code{\link{Model}}
setMethodS3("Train", "Model", conflict="quiet",
  function(this, d, force.retrain=FALSE, ...) {
    if (force.retrain || !d$Same(d=this$.last.trained.d,
        compare=c("X", "xformedDpts", "noiseVar"))) {
      old.params <- this$params
      lower <- this$getLower(for.training=TRUE)
      upper <- this$getUpper(for.training=TRUE)
      params <- this$getParams(for.training=TRUE)
      if (length(params) > 1) {
        this$.opt <- optim(method="L-BFGS-B", control=list(fnscale=-1),
          par=params, lower=lower, upper=upper,
          fn=LogML, gr=GradLogML,
          # Extra parameters needed by 'fn' and 'gr':
          model=this, d=d)
        if (this$.opt$convergence > 50) {  # L-BFGS-B warning or error
          warning(paste(sep='',
              "L-BFGS-B optimization ran into trouble; params NOT changed:\n    ",
              this$.opt$message))
          cat("By the way, the params were as follows:\n        ", this$.opt$par, "\n")
          this$params <- old.params
        } else {
          this$.last.trained.d <- clone(d)
        }
      } else if (length(params) == 1) {
        opt <- optimize(f=LogML1D, maximum=TRUE,
          lower=lower, upper=upper,
          model=this, d=d, name=names(params))
        best <- opt$maximum
        names(best) <- names(params)
        this$params <- best
      }
    }
    return (invisible(this))
  })

#-------------------------------------------------------------------------------
# (Model) PRIVATE METHODS:

setMethodS3("CheckContributionsAndWarn", "Model", private=TRUE, conflict="quiet",
  function(this, contributions, ...) {
    # Check the list of names in 'contributions' to see which we actually have,
    # returning a validated list of names, and warning if any were requested
    # which do not exist, or if no valid contributions remain.
    #
    # Args:
    #   contributions:  A character vector with the names of the contributions
    #      which we should check.
    #
    # Returns:
    #   A character vector of the unique contribution names which actually
    #   exist.

    # Check which requested contributions are present in the model.
    contributions <- unique(contributions)
    existing <- which(contributions %in% this$contributionIds)
    bad.names <- contributions[-existing]
    if (length(bad.names) > 0) {
      culprits <- paste(sep='', collapse=' ', '"', bad.names, '"')
      warning(paste("The following contribution IDs do not exist:\n",
          culprits, "\n"))
    }
    contributions <- contributions[existing]

    # Make sure we have at least one contribution!
    if (length(contributions) < 0) {
      warning("No nonzero contributions actually supplied!\n")
    }

    return (contributions)
  })

setMethodS3("ComputeL", "Model", private=TRUE, conflict="quiet",
  function(this, d, X.out=d$X, contributions, ...) {
    # Compute the Cholesky decomposition of the model's current covariance
    # matrix for datapoints 'd'.
    #
    # Args:
    #   d:  The Dataset for evaluating this model.
    #   X.out:  matrix; the X-points where we want to predict the function.
    #   contributions:  character vector; id's of the contributing Covariances
    #      we want to predict (other contributions are considered noise);
    #      default is every contribution not named 'noise'.
    #
    # Returns:
    #   The Cholesky decomposition of the output points' covariance matrix for
    #   this model and dataset.
    ingredients <- list(X=d$X, X.out=X.out, noiseVar=d$noiseVar,
      params=this$params)
    if (this$.L$NeedToRecalculate(ingredients=ingredients)) {
      this$ComputeKChol(d=d)
      K.out.out <- d$noiseVar * diag(NumPoints(X.out))
      K.in.out <- matrix(0, nrow=NumPoints(X.out), ncol=NumPoints(d$X))
      for (c.id in contributions) {
        covar <- this$.contributions[[c.id]]
        K.out.out <- K.out.out + covar$KOutOut(X.out=X.out)
        K.in.out <- K.in.out + covar$KInOut(d=d, X.out=X.out)
      }
      K.posterior <- K.out.out - (K.in.out %*% chol2inv(this$.K.chol$M) %*%
        t(K.in.out))
      this$.L$StoreMatrix(M=t(chol(K.posterior)), ingredients=ingredients)
    }
    return (invisible(this))
  })

setMethodS3("ComputeKChol", "Model", private=TRUE, conflict="quiet",
  function(this, d, ...) {
    # Compute the Cholesky decomposition of the model's current covariance
    # matrix for datapoints 'd'.
    #
    # Args:
    #   d:  The Dataset for evaluating this model.
    #
    # Returns:
    #   The Cholesky decomposition of the total covariance matrix for this
    #   model.
    ingredients <- list(X=d$X, noiseVar=d$noiseVar, params=this$params)
    if (this$.K.chol$NeedToRecalculate(ingredients=ingredients)) {
      K.tot <- this$KTotal(d=d)
      K.chol <- DebugIfError(chol.default, K.tot)
      this$.K.chol$StoreMatrix(M=K.chol, ingredients=ingredients)
    }
    return (invisible(this))
  })

setMethodS3("KDeriv", "Model", private=TRUE, conflict="quiet",
  function(this, d, param, ...) {
    # Computes the element-by-element derivative of the total K-matrix for
    # Dataset 'd', with respect to the parameter named 'param'.
    #
    # Args:
    #   d:  The Dataset for evaluating this model.
    #   param:  The name of the parameter with respect to which we're
    #      differentiating .
    #
    # Returns:
    #   The element-by-element derivative of the total K-matrix for Dataset
    #   'd', with respect to the parameter named 'param'.
    K.deriv <- matrix(0, nrow=d$n, ncol=d$n)
    in.logspace <- (length(grep(pattern=LogspaceTag(), x=param)) > 0)
    for (covar in this$.contributions) {
      if (in.logspace) {  # Decode the param name, if necessary
        param <- sub(pattern=LogspaceTag(), replacement="", x=param)
      }
      d.covar_d.param <- covar$KInInDeriv(d=d, param=param)
      if (in.logspace) {  # d/d(log(x)) = x*d/d(x)
        d.covar_d.param <- d.covar_d.param * this$params[param]
      }
      K.deriv <- K.deriv + d.covar_d.param
    }
    return (K.deriv)
  })

setMethodS3("K", "Model", private=TRUE, conflict="quiet",
  function(this, X, X.out=X, contributions, ...) {
    # Calculate a covariance matrix from X to X.out, including only the named
    # contributions.
    #
    # Args:
    #   X:  2-column numeric matrix; the input points (where we have data).
    #   X.out:  2-column numeric matrix; The output points (where we wish to
    #      make predictions).
    #   contributions:  character vector; a list of names of covariances to
    #      include.
    #
    # Returns:
    #   The covariance matrix from X to X.out.
    K <- matrix(0, nrow=NumPoints(X.out), ncol=NumPoints(X))
    # Add covariance matrix for each additive contribution:
    for (covar in this$.contributions) {
      if (covar$id %in% contributions) {
        K <- K + covar$K.specific(X=X, X.out=X.out)
      }
    }
    return (K)
  })

setMethodS3("KTotal", "Model", private=TRUE, conflict="quiet",
  function(this, d, ...) {
    # Compute the total covariance matrix for this model with respect to the
    # data in Dataset 'd'.
    #
    # Args:
    #   d:  The Dataset for evaluating this model.
    #
    # Returns:
    #   The total covariance matrix for this model.
    K.tot <- matrix(0, nrow=d$n, ncol=d$n)
    for (covar in this$.contributions) {
      K.tot <- K.tot + covar$KInIn(d=d)
      gc()
    }
    # Cap it off with the noise associated with the data:
    diag(K.tot) <- diag(K.tot) + d$noiseVar
    return (K.tot)
  })

setMethodS3("KInv", "Model", private=TRUE, conflict="quiet",
  function(this, d, ...) {
    # The inverse total covariance matrix for this model, with respect to the
    # given dataset.
    #
    # Args:
    #   d:  The Dataset for evaluating this model.
    #
    # Returns:
    #   The inverse total covariance matrix for this model.
    this$ComputeKChol(d=d)
    return (chol2inv(this$.K.chol$M))
  })

setMethodS3("LogDetK", "Model", private=TRUE, conflict="quiet",
  function(this, d, ...) {
    # Compute the logarithm of the determinant of the model's total covariance
    # matrix for the points in Dataset 'd'.
    #
    # Args:
    #   d:  Dataset at whose points we evaluate the covariance matrix.
    #
    # Returns:
    #   The logarithm of the determinant of the model's covariance matrix.
    this$ComputeKChol(d=d)
    return (2 * sum(log(diag(this$.K.chol$M))))
  })

