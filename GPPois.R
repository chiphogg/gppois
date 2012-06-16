# GPPois.R: Gaussian Process-based inference for Poisson-noised data.
# Copyright (c) 2011 Charles R. Hogg III
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Author: Charles R. Hogg III (2011) <charles.r.hogg@gmail.com>

# FILE DESCRIPTION:
# Provides classes and functions to analyze data using Gaussian Processes.  Key
# examples include:
#  - Dataset: the data to be analyzed
#  - Covariance: (virtual) superclass for the different possible kinds of
#      covariance structures (e.g. Squared-Exponential, periodic, Matern, ...)
#  - Model: an interface to the total covariance structure (different
#      Covariance terms contribute *additively*)
#
# These objects *need* to be mutable, so I am using the ''R.oo'' library to
# enable pass-by-reference.
#
# A lot of computations result in really big matrices.  These should be easy to
# save and restore to disk.
#
# Also, I am trying to adhere to Google's R Style Guide as much as possible:
# http://google-styleguide.googlecode.com/svn/trunk/google-r-style.html

# Object-orientation (with pass-by-reference!!)
library("R.oo")

################################################################################
# CLASS:                             Model
#
# A collection of covariance structures which can be trained on datasets and
# make predictions.
#
# Helper Functions:
#   LogML:  The logarithm of the (M)arginal (L)ikelihood for this model, given
#      a Dataset.
#   GradLogML:  The gradient (w.r.t. the parameter values) of the logarithm of
#      the (M)arginal (L)ikelihood for this model, given a Dataset.
#
# Virtual Fields (R = read only):
# (R) contributionIds:  The id strings of the Covariance objects which
#        contribute to this Model.
# ( ) id:  A character string used to identify this model.
# ( ) lower:  A named vector of lower bounds on parameters governing this
#     Model.
# ( ) params:  A named vector of parameters governing this Model.
# (R) signalIds:  The ids of contributions which are considered to be "signal",
#        i.e., not noise.
# ( ) upper:  A named vector of upper bounds on parameters governing this
#     Model.
# (R) varyingParamNames:  Names of parameters with nonzero range between upper
#     and lower bounds.
#
# Methods:
#   AddCovariance:  Add a new Covariance object to this Model.
#   clone:  Deep clone of this Model object (clones Covariance's as well).
#   Forget:  Clear all large matrices to save memory.
#   Freeze:  Make all parameter values constant.
#   L:  Lower-triangular Cholesky decomposition of the covariance matrix
#      (useful for generating random draws).
#   NamedCovariance:  Retrieve a clone of the first contributing Covariance
#      object with the given id.
#   PlotBubblingSurfaces2D:  Plot smoothly varying random surfaces to visualize
#      the uncertainty.
#   PlotCovariance2D:  Plot the covariance of points in 2D with respect to
#      another point; also plot the derivatives w.r.t. each parameter.
#   PosteriorInterval:  bounds on the uncertainty (in terms of the standard
#      deviation) in the prediction at a given point.
#   PosteriorMean:  The optimal prediction of the underlying function at a set
#      of points.
#   PosteriorStandardDeviation:  pointwise sigma (for *transformed* values) at
#      each prediction point.
#   PredictionMatrix:  A matrix relating function values at the output points
#      to function values at the input points.
#   print:  Prettied-up summary of this Model object.
#   SetNoiseBounds:  Add a 'noise' contribution with the given bounds.
#   Summary:  A data.frame summarizing predictions for multiple model subsets.
#   Train:  Optimize this Model's parameters to describe the given data.
#
# This class should work just fine, as long as
#   a) we are training on all the datapoints together (i.e., not breaking them
#      up into subregions to divide-and-conquer), and
#   b) this$params returns a vector which is amenable to simple optimization
#      routines (i.e., none of the parameters require special treatment).
# If either of these conditions fail, a specialized subclass should be created
# instead.  Note that *both* these conditions fail for *both* scenarios
# considered in our Journal of Applied Crystallography paper, despite the fact
# that I wrote this software to perform the analysis for that paper.  That's
# OK; having these classes still makes it very much easier to build the
# specialized functions I need.

setConstructorS3("Model",
  function(id="") {
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

#-------------------------------------------------------------------------------
# (Model) PUBLIC "HELPER" FUNCTIONS:

LogML <- function(par=model$getParams(for.training=TRUE), model, d,
  update.params=TRUE) {
  # Sets the parameter values to 'par' for Model 'model', and returns the log
  # of the (M)arginal (L)ikelihood for describing Dataset 'd'.
  #
  # Args:
  #   par:  The parameter values to test.
  #   model:  The Model object we're optimizing.
  #   d:  The Dataset we're training on.
  #   update.params:  logical; if TRUE, we should change model's params to the
  #      values in par.
  #
  # Returns:
  #   The log of the marginal likelihood (also has a side-effect of setting
  #   model$params <- par, unless update.params is FALSE (although this would
  #   create an extra copy of model, and thus probably be less efficient).
  if (!update.params) {
    model <- clone(model)
  }
  model$setParams(p=DecodeForTraining(par))
  Y <- d$xformedDpts
  # The following calculation is based on Equation 5.8 in Rasmussen and
  # Williams, "Gaussian Processes for Machine Learning".
  term.data.fit   <- -0.5 * t(Y) %*% model$KInv(d) %*% Y
  term.complexity <- -0.5 * model$LogDetK(d)
  term.num.dpts   <- -0.5 * d$n * log(2 * pi)
  return (term.data.fit + term.complexity + term.num.dpts)
}

LogML1D <- function(value, name, model, d, update.params=TRUE) {
  # Calculate the log of the marginal likelihood for Model 'model' on Dataset
  # 'd', if model parameter given by 'name' is replaced by 'value'.
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

GradLogML <- function(par=model$getParams(for.training=TRUE), model, d,
  update.params=TRUE) {
  # Sets the parameter values to 'par' for Model 'model', and returns the
  # gradient (w.r.t. par) of the log of the (M)arginal (L)ikelihood for
  # describing Dataset 'd'.
  #
  # Args:
  #   par:  The parameter values to test.
  #   model:  The Model object we're optimizing.
  #   d:  The Dataset we're training on.
  #   update.params:  logical; if TRUE, we should change model's params to the
  #      given value.
  #
  # Returns:
  #   The gradient of the log of the marginal likelihood (also has a
  #   side-effect of setting model$params <- par).
  if (!update.params) {
    model <- clone(model)
  }
  model$setParams(p=DecodeForTraining(par))
  Y <- d$xformedDpts
  # The following calculations are based on Equation 5.9 in Rasmussen and
  # Williams, "Gaussian Processes for Machine Learning".
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

#-------------------------------------------------------------------------------
# (Model) PUBLIC VIRTUAL FIELDS:

# Virtual Field: contributionIds
# A list of id's for all contributions in this Model.
#
# Returns:
#   The id's for this Model's contributing Covariance objects.
setMethodS3("getContributionIds", "Model", conflict="quiet",
  function(this, ...) {
    id.list <- c()
    for (covar in this$.contributions) {
      id.list <- c(id.list, covar$id)
    }
    return (id.list)
  })

# Virtual Field: id
# A character string identifying this Model object.
#
# Args:
#   id: the string to change the id to.
#
# Returns:
#   The id of this Model object.
setMethodS3("getId", "Model", conflict="quiet",
  function(this, ...) {
    return (this$.id)
  })
setMethodS3("setId", "Model", conflict="quiet",
  function(this, id, ...) {
    this$.id <- id
    return (this)
  })

# Virtual Field: params
# A named vector of parameters governing this Model.
#
# Args:
#   p: A (named) vector of new parameter values (we ONLY use ones which are
#      named, and whose names match up with names of parameters.)
#   for.training:  If TRUE, we ignore "constant" parameters (i.e., where
#      lower=upper) and return the *log* of any "scale" parameters (such as ell
#      or sigma.f for the SE model).
#
# Returns:
#   A vector with values for ell and sigma.f
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

# Virtual Field: lower
# Lower bounds for the parameter values.
#
# Args:
#   L: A (named) vector of new lower bounds (we ONLY use ones which are named,
#      and whose names match up with names of parameters.)
#   for.training:  If TRUE, we ignore "constant" parameters (i.e., where
#      lower=upper) and return the *log* of any "scale" parameters (such as ell
#      or sigma.f for the SE model).
#
# Returns:
#   The lower bounds for the parameters for this model.
setMethodS3("getLower", "Model", conflict="quiet",
  function(this, for.training=FALSE, ...) {
    L <- c()
    for (covar in this$.contributions) {
      L <- c(L, covar$getLower(for.training=for.training))
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

# Virtual Field: upper
# Upper bounds for the parameter values.
#
# Args:
#   U: A (named) vector of new upper bounds (we ONLY use ones which are named,
#      and whose names match up with names of parameters.)
#   for.training:  If TRUE, we ignore "constant" parameters (i.e., where
#      lower=upper) and return the *log* of any "scale" parameters (such as ell
#      or sigma.f for the SE model).
#
# Returns:
#   The upper bounds for the parameters for this model.
setMethodS3("getUpper", "Model", conflict="quiet",
  function(this, for.training=FALSE, ...) {
    U <- c()
    for (covar in this$.contributions) {
      U <- c(U, covar$getUpper(for.training=for.training))
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

# Virtual Field: signalIds (read-only)
# A list of id's for all non-noise contributions in this Model.
#
# Returns:
#   The id's for this Model's non-noise contributing Covariance objects.
setMethodS3("getSignalIds", "Model", conflict="quiet",
  function(this, ...) {
    # The ids of contributions which are considered to be "signal", i.e., not
    # noise.
    return (this$getContributionIds()[
        which(this$getContributionIds() != 'noise')])
  })

# Virtual Field: varyingParamNames
# Names of the parameters which are not constant
#
# Args:
#   U: A (named) vector of new upper bounds (we ONLY use ones which are named,
#      and whose names match up with names of parameters.)
#
# Returns:
#   The upper bounds for the parameters for this model.
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


#-------------------------------------------------------------------------------
# (Model) PUBLIC METHODS:

setMethodS3("AddCovariance", "Model", conflict="quiet",
  function(this, covariance, on.duplicate.id="rename", ...) {
    # Add another Covariance structure to this model (it contributes
    # additively).
    #
    # Args:
    #   covariance:  A Covariance object to be cloned and added to this model.
    #   on.duplicate.id:  Character string, one of ("rename", "replace"),
    #      directing how to handle a new contribution whose ID is the same as
    #      an existing one's.
    #
    # Returns:
    #   Used for its side-effect.

    # We CLONE it so we OWN it.  (Don't want anyone else to fiddle with the
    # Covariance object, EXCEPT this Model object.)
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

clone.Model <- function(this, ...) {
  # Deep clone this Model object (clones Covariance's too).
  #
  # Returns:
  #   A deep clone of the model object.
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

setMethodS3("Forget", "Model", conflict="quiet",
  function(this, ...) {
    this$.K.chol <- LazyMatrix()
    this$.L <- LazyMatrix()
    this$.last.trained.d <- NA
    gc()
  })

setMethodS3("Freeze", "Model", conflict="quiet",
  function(this, p.names=names(this$params), ...) {
    good.i <- which(p.names %in% names(this$params))
    for (p.name in p.names[good.i]) {
      this$lower <- this$params[p.name]
      this$upper <- this$params[p.name]
    }
  })

setMethodS3("L", "Model", conflict="quiet",
  function(this, d, X.out=d$X, contributions=this$getSignalIds(), ...) {
    # Compute the lower-triangular Cholesky decomposition of the covariance
    # matrix (useful for generating random draws).
    # 
    # Args:
    #   d:  The Dataset to evaluate the Model on.
    #   X.out:  matrix; the X-points where we want to predict the function.
    #   contributions:  character vector; id's of the contributing Covariances
    #      we want to predict (other contributions are considered noise);
    #      default is every contribution not named 'noise'.
    #
    # Returns:
    #   The lower-triangular Cholesky decomposition of the covariance matrix.
    this$ComputeL(d=d, X.out=X.out, contributions=contributions)
    return (this$.L$M)
  })

setMethodS3("NamedCovariance", "Model", conflict="quiet",
  function(this, id, ...) {
    # Retrieve a clone of the first contributing Covariance objects with the
    #   given id.
    #
    # Args:
    #   id:  The id (character) of the covariance to retrieve and clone.
    #
    # Returns:
    #   A clone of the Covariance object with the given id (or else NA if no
    #   such object exists).
    for (covar in this$.contributions) {
      if (covar$id == id) {
        Cov <- clone(covar)
        return (Cov)
      }
    }
    return (NA)
  })

setMethodS3("PlotBubblingSurfaces2D", "Model", conflict="quiet",
  function(this, d, X.out=d$X, contributions=this$getSignalIds(),
    n.surfaces=10, n.times=50, file.name=NULL, ...) {
    # Plot smoothly varying random surfaces to visualize the uncertainty.
    #
    # Args:
    #   d:  The Dataset to evaluate the Model on.
    #   X.out:  matrix; the X-points where we want to predict the function.
    #   contributions:  character vector; id's of the contributing Covariances
    #      we want to predict (other contributions are considered noise).
    #   n.surfaces:  The number of independent draws to take for each bubbler.
    #   n.times:  The final number of interpolated points.
    #   file.name:  Character, the basename of the file to plot to.
    #
    # Returns:
    #   Used for its side-effect (plotting a series of timesteps to PNG files).
    B <- BubblingRandomMatrix(n.pts=NumPoints(X.out), N=n.surfaces,
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

setMethodS3("PlotCovariance2D", "Model", conflict="quiet",
  function(this, d, file.name=NULL, file.type='png', i=1, ...) {
    # Plot the covariance of points in 2D with respect to another point; also
    # plot the derivatives w.r.t. each parameter.
    #
    # Args:
    #   d:  The Dataset to evaluate the Model on.
    #   file.name:  Character, the name of the file to plot to (file extension
    #      should be automatically appended if necessary)
    #   file.type:  Character, one of ("pdf", "png")
    #   i:  numeric; index of datapoint to use for calculating the covariance
    #
    # Returns:
    #   Used for its side-effects of plotting.
    #
    # Notes:
    #   Lay them out in a grid of two rows.  First column: top is covariance
    #   matrix, bottom is for text (parameter values).  Subsequent columns show
    #   derivative matrices.

    # Set the layout:
    n.rows <- 2
    n.cols <- 1 + ceiling(length(this$params) / n.rows)
    Layout <- grid.layout(nrow = n.rows, ncol = n.cols, 
      widths  = unit(rep(1, n.cols), rep("null", n.cols)), 
      heights = unit(rep(1, n.rows), rep("null", n.rows)))

    # Open the file:
    my.file <- SetupFileInfo(name=file.name, type=file.type)
    if (my.file$type == 'png') {
      if (!require("Cairo")) {
        png(filename=my.file$name, ...)
      } else {
        CairoPNG(filename=my.file$name, ...)
      }
    } else if (my.file$type == 'pdf') {
      cairo_pdf(filename=my.file$name, ...)
    } else {
      stop("Sorry, we can only do pdf or png output for now.")
    }
    LayoutNewGridPage(Layout=Layout)

    # Construct a data.frame for all the matrices
    point.out <- matrix(d$X[i, ], nrow=1)
    K.data <<- data.frame(X.1=d$X[, 1], X.2=d$X[, 2],
      Cov=this$KTotal(d=d)[i, ],
      M=t(this$PredictionMatrix(d=d, X.out=point.out,
          contributions=this$contributionIds[
            which(this$contributionIds != 'noise')]))
      )

    # Generic ggplot setup
    ps.out <- 1.3  # outer point size
    ps.in <- 1.0   # inner point size
    p.base <- (ggplot(data=K.data, aes(x=X.1, y=X.2))
      + scale_colour_gradientn(colours=c('white','red','blue'))
      + geom_point(colour='black', size=ps.out)
      + geom_vline(xintercept=K.data[i, 1])
      + geom_hline(yintercept=K.data[i, 2])
      )

    # Print the covariance matrix in the top left
    p.cov <- (p.base
      + geom_point(aes(colour=Cov), size=ps.in)
      )
    print(p.cov, vp=Subplot(1, 1))
    # Weight matrix in bottom left
    p.weight <- (p.base
      + geom_point(aes(colour=M), size=ps.in)
      + scale_colour_gradient2(colours=c('red', 'white', 'blue'))
      )
    print(p.weight, vp=Subplot(2, 1))

    dev.off()
    return (invisible(this))
  })

setMethodS3("PosteriorMean", "Model", conflict="quiet",
  function(this, d, X.out=d$X, contributions=this$getSignalIds(),
    untransform.result=TRUE, ...) {
    # Computes this Model's optimal prediction of the underlying function's
    # value at every point in 'X.out'.
    #
    # Args:
    #   d:  The Dataset to train the Model on.
    #   X.out:  matrix; the X-points where we want to predict the function.
    #   contributions:  character vector; id's of the contributing Covariances
    #      we want to predict (other contributions are considered noise);
    #      default is every contribution not named 'noise'.
    #   untransform.result:  logical; if TRUE, we transform back to the space
    #      of dpts (as opposed to the space of xformedDpts where training takes
    #      place).
    #
    # Returns:
    #   A numeric vector with optimal predictions at every point in X.out.
    contributions <- this$CheckContributionsAndWarn(contributions)
    M <- this$PredictionMatrix(d=d, X.out=X.out, contributions=contributions)
    result <- M %*% d$xformedDpts
    if (untransform.result) {
      result <- d$Untransform(result)
    }
    return (result)
  })

setMethodS3("PredictionMatrix", "Model", conflict="quiet",
  function(this, d, X.out=d$X, contributions=this$getSignalIds(), ...) {
    # A matrix relating function values at the output points to function values
    # at the input points.
    #
    # Args:
    #   d:  The Dataset to train the Model on.
    #   X.out:  matrix; the X-points where we want to predict the function.
    #   contributions:  character vector; id's of the contributing Covariances
    #      we want to predict (other contributions are considered noise).
    #
    # Returns:
    #   The matrix 'M' such that (Y.pred = M %*% d$xformedDpts).
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

setMethodS3("PosteriorInterval", "Model", conflict="quiet",
  function(this, d, X.out=d$X, num.sd=1, contributions=this$getSignalIds(),
    ...) {
    # Computes bounds on the uncertainty (in terms of the standard deviation)
    # in the prediction at a given point, along with the prediction.
    #
    # Args:
    #   d:  The Dataset to train the Model on.
    #   X.out:  matrix; the X-points where we want to predict the function.
    #   num.sd:  The number of standard deviations from the mean our interval
    #      should include.
    #   contributions:  character vector; id's of the contributing Covariances
    #      we want to predict (other contributions are considered noise);
    #      default is every contribution not named 'noise'.
    #
    # Returns:
    #   A numeric vector with optimal predictions at every point in X.out.
    M <- this$PredictionMatrix(d=d, X.out=X.out, contributions=contributions)
    std.dev <- this$PosteriorStandardDeviation(d=d, X.out=X.out,
      contributions=contributions, ...)
    prediction <- M %*% d$xformedDpts
    return (data.frame(X=X.out,
        mean=d$Untransform(prediction),
        lower=d$Untransform(prediction - num.sd * std.dev),
        upper=d$Untransform(prediction + num.sd * std.dev)))
  })

setMethodS3("PosteriorStandardDeviation", "Model", conflict="quiet",
  function(this, d, X.out=d$X, contributions=this$getSignalIds(), ...) {
    # Computes the posterior "sigma" at a given point.
    #
    # Args:
    #   d:  The Dataset to train the Model on.
    #   X.out:  matrix; the X-points where we want to predict the function.
    #   contributions:  character vector; id's of the contributing Covariances
    #      we want to predict (other contributions are considered noise);
    #      default is every contribution not named 'noise'.
    #
    # Returns:
    #   A numeric vector with the posterior standard deviation at every point
    #   in X.out.
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

print.Model <- function(this, indent=0, ...) {
  # Pretty-prints information for this Model object.
  #
  # Returns:
  #   Used for its side-effect.
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

setMethodS3("SetNoiseBounds", "Model", conflict="quiet",
  function(this, sigma.vals, ...) {
    # Set the range of values for the noise level in this model.
    #
    # Args:
    #   sigma.vals:  A numeric vector, whose range() sets the range of values
    #   for the noise level.
    #
    # Returns:
    #   Used for its side-effect.
    this$AddCovariance(CovarianceNoise(id="noise", sigma.bounds=sigma.vals),
      on.duplicate.id="replace")
  })

setMethodS3("Train", "Model", conflict="quiet",
  function(this, d, force.retrain=FALSE, ...) {
    # Optimize this Model's parameters so they describe the given data.
    #
    # Args:
    #   d:  Dataset; the data which our parameters should describe.
    #
    # Returns:
    #   Used for its side-effect (i.e., changing the values of this Model's
    #   parameters to reflect the given Dataset).
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

################################################################################
# SECTION:                    Standalone Functions
################################################################################

