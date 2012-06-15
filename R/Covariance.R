############################################################################/**
# @RdocClass Covariance
#
# @title "Superclass for all Covariance function objects"
#
# \description{
#   @classhierarchy
#
#   @get "title".
#
#   NB: Nobody will actually make an object of type \code{Covariance}; instead,
#   they will use one of the derived classes (such as
#   \code{link{CovarianceSE}}).
# }
#
# @synopsis
#
# \arguments{
#   \item{id}{(character) A string to identify this covariance object.}
#   \item{...}{Not used.}
# }
#
# \details{
#    Regarding parameter names: the SUBCLASS has the responsibility to provide
#    "plain-named" versions of all the virtual fields (paramNamesPlain,
#    paramsPlain, lowerPlain, etc.).  The SUPERCLASS will automatically handle
#    the "decorated" versions (paramNames, params, lower, etc.).
#
#    All Covariance subclasses should remember the param-values and Dataset
#    they last used to compute their K-matrix.  If they get asked to compute it
#    again, they will simply return the previously-computed result if these
#    values have not changed.
# }
#
# \section{Fields and Methods}{
#  @allmethods
#
# }
#
# @author
#*/###########################################################################
setConstructorS3("Covariance",
  function(id="", ...) {
    extend(Object(), "Covariance",
      .id = id
      )
  })

############################################################################/**
# @RdocMethod getId
#
# @title "The ID for this Covariance"
#
# \description{
#   @get "title":
#   a character string identifying this Covariance object.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#   The id of this Covariance object.
# }
#
# \seealso{
#   @seemethod "setId"
#   @seeclass
# }
#
# @alias id.Covariance
#
# @author
#*/###########################################################################
setMethodS3("getId", "Covariance", conflict="quiet",
  function(this, ...) {
    return (this$.id)
  })

############################################################################/**
# @RdocMethod setId
#
# @title "The ID for this Covariance"
#
# \description{
#   @get "title":
#   a character string identifying this Covariance object.
# }
#
# @synopsis
#
# \arguments{
#   \item{id}{A character string identifying this Covariance object.}
# }
#
# \seealso{
#   @seemethod "getId"
#   @seeclass
# }
#
# @author
#*/###########################################################################
setMethodS3("setId", "Covariance", conflict="quiet",
  function(this, id, ...) {
    this$.id <- id
    return (this)
  })

############################################################################/**
# @RdocMethod getParamNames
#
# @title "Decorated names of the covariance parameters"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#   Parameter names in the form "id.basename", where basename is the
#   "undecorated" parameter name.
# }
#
# \seealso{
#   @seemethod getParams
#   @seeclass
# }
#
# @author
#*/###########################################################################
setMethodS3("getParamNames", "Covariance", conflict="quiet",
  function(this, ...) {
    if (!is.character(getId(this)) || nchar(getId(this)) < 1) {
      return (getParamNamesPlain(this))
    }
    return (this$PrependId(getParamNamesPlain(this)))
  })

############################################################################/**
# @RdocMethod getParams
#
# @title "Parameter values for this Covariance"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{for.training}{(logical) If TRUE, we return values more suitable for
#      optimization: "scale"-type parameters are given as logarithms (which
#      amounts to using the Jeffreys prior), and "constant" params (i.e., where
#      the lower bound equals the upper bound) are omitted entirely.}
#   \item{...}{Not used.}
# }
#
# \value{
#   (named numeric vector) Parameter values for this Covariance.
# }
#
# \seealso{
#   @seemethod setParams
#   @seemethod getParamsPlain
#   @seeclass
# }
#
# @alias params.Covariance
#
# @author
#*/###########################################################################
setMethodS3("getParams", "Covariance", conflict="quiet",
  function(this, for.training=FALSE, ...) {
    p <- getParamsPlain(this)
    names(p) <- getParamNames(this)
    if (for.training) {
      p <- this$EncodeForTraining(p)
    }
    return (p)
  })

############################################################################/**
# @RdocMethod setParams
#
# @title "Parameter values for this Covariance"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{p}{(named numeric vector) New parameter values. (N.B: we \emph{only}
#      use ones which are named, and whose names match up with names of
#      parameters.)}
#   \item{for.training}{(logical) If TRUE, we return values more suitable for
#      optimization: "scale"-type parameters are given as logarithms (which
#      amounts to using the Jeffreys prior), and "constant" params (i.e., where
#      the lower bound equals the upper bound) are omitted entirely.}
#   \item{...}{Not used.}
# }
#
# \seealso{
#   @seemethod setParams
#   @seemethod getParamsPlain
#   @seeclass
# }
#
# @alias params.Covariance
#
# @author
#*/###########################################################################
setMethodS3("setParams", "Covariance", conflict="quiet",
  function(this, p, for.training=FALSE, ...) {
    if (for.training) {
      p <- this$DecodeForTraining(p)
    }
    p.plain <- this$UndecorateNames(p=p)
    this$setParamsPlain(p=p.plain)
    return (invisible(this))
  })

# Virtual Field: paramsPlain
# Gives a vector of parameter values, whose names are NOT decorated by the id
# of this Covariance object.
#
# Args:
#   p: A (named) vector of new parameter values (we ONLY use ones which are
#      named, and whose names match up with names of parameters.)
#
# Returns:
#   The parameters for this covariance function, but with names undecorated by
#   its id.
#
# NOTE: accessor method is responsibility of subclasses.

############################################################################/**
# @RdocMethod setParamsPlain
#
# @title "Set parameters using undecorated names"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{p}{A (named) vector of new parameter values (we \emph{only} use ones
#      which are named, and whose names match up with names of parameters.)}
#   \item{...}{Not used.}
# }
#
# \seealso{
#   @seemethod getParamsPlain
#   @seeclass
# }
#
# @author
#*/###########################################################################
setMethodS3("setParamsPlain", "Covariance", conflict="quiet",
  function(this, p, ...) {
    # First, we need a vector where every parameter name has a value: either
    # the new value from 'p', or if none was given, the current value.
    p.good.names <- this$getParamsPlain()
    to.change <- names(p)[which(names(p) %in% names(p.good.names))]
    p.good.names[to.change] <- p[to.change]
    # paramsPlainImplementation requires a vector with a value for every
    # parameter; we took care of this above.
    this$paramsPlainImplementation(p=p.good.names)
    p.clamped <- this$ClampedParamVals(p=p.good.names)
    p.names <- names(p.clamped)
    clamped <- which(p.good.names[p.names] != p.clamped[p.names])
    this$paramsPlainImplementation(p=p.clamped)
    if (any(clamped)) {
      culprits <- paste(sep='', collapse=' ', '"', p.names[clamped], '"')
      warning(paste("These parameters had to be clamped:\n", culprits, "\n"))
    }
    return (invisible(this))
  })

############################################################################/**
# @RdocMethod getLower
#
# @title "Lower bounds on each parameter vector"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{for.training}{(logical) If TRUE, we return values more suitable for
#      optimization: "scale"-type parameters are given as logarithms (which
#      amounts to using the Jeffreys prior), and "constant" params (i.e., where
#      the lower bound equals the upper bound) are omitted entirely.}
#   \item{...}{Not used.}
# }
#
# \value{
#   (named numeric vector) The lower bounds for the parameters for this
#   covariance function.
# }
#
# \seealso{
#   @seemethod setLower
#   @seemethod getUpper
#   @seeclass
# }
#
# @author
#*/###########################################################################
setMethodS3("getLower", "Covariance", conflict="quiet",
  function(this, for.training=FALSE, ...) {
    L <- this$getLowerPlain()
    names(L) <- this$getParamNames()
    if (for.training) {
      L <- this$EncodeForTraining(L)
    }
    return (L)
})

############################################################################/**
# @RdocMethod setLower
#
# @title "Lower bounds on each parameter vector"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{L}{A (named) vector of new lower bounds (we \emph{only} use ones
#      which are named, and whose names match up with names of parameters.)}
#   \item{for.training}{(logical) If TRUE, we return values more suitable for
#      optimization: "scale"-type parameters are given as logarithms (which
#      amounts to using the Jeffreys prior), and "constant" params (i.e., where
#      the lower bound equals the upper bound) are omitted entirely.}
#   \item{...}{Not used.}
# }
#
# \seealso{
#   @seemethod getLower
#   @seemethod setUpper
#   @seeclass
# }
#
# @author
#*/###########################################################################
setMethodS3("setLower", "Covariance", conflict="quiet",
  function(this, L, for.training=FALSE, ...) {
    if (for.training) {
      L <- DecodeForTraining(L)
    }
    L.plain <- this$UndecorateNames(L)
    this$setLowerPlain(L=L.plain)
    return (invisible(this))
})

############################################################################/**
# @RdocMethod getUpper
#
# @title "Upper bounds for the parameter values"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{for.training}{(logical) If TRUE, we return values more suitable for
#      optimization: "scale"-type parameters are given as logarithms (which
#      amounts to using the Jeffreys prior), and "constant" params (i.e., where
#      the lower bound equals the upper bound) are omitted entirely.}
#   \item{...}{Not used.}
# }
#
# \value{
#   (named numeric vector) The upper bounds for the parameters for this
#   covariance function.
# }
#
# \seealso{
#   @seemethod getLower
#   @seemethod setUpper
#   @seeclass
# }
#
# @author
#*/###########################################################################
setMethodS3("getUpper", "Covariance", conflict="quiet",
  function(this, for.training=FALSE, ...) {
    U <- this$getUpperPlain()
    names(U) <- this$getParamNames()
    if (for.training) {
      U <- this$EncodeForTraining(U)
    }
    return (U)
  })

############################################################################/**
# @RdocMethod setUpper
#
# @title "Upper bounds for the parameter values"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{U}{A (named) vector of new upper bounds (we \emph{only} use ones
#      which are named, and whose names match up with names of parameters.)}
#   \item{for.training}{(logical) If TRUE, we return values more suitable for
#      optimization: "scale"-type parameters are given as logarithms (which
#      amounts to using the Jeffreys prior), and "constant" params (i.e., where
#      the lower bound equals the upper bound) are omitted entirely.}
#   \item{...}{Not used.}
# }
#
# \seealso{
#   @seemethod setLower
#   @seemethod getUpper
#   @seeclass
# }
#
# @author
#*/###########################################################################
setMethodS3("setUpper", "Covariance", conflict="quiet",
  function(this, U, for.training=FALSE, ...) {
    if (for.training) {
      U <- this$DecodeForTraining(U)
    }
    U.plain <- this$UndecorateNames(U)
    this$setUpperPlain(U=U.plain)
    return (invisible(this))
  })

###########################################################################/**
# @RdocMethod clone
#
# @title "Deep-clones a Dataset"
#
# \description{
#   @get "title".  (Clones LazyMatrix objects too.)
# }
#
# @synopsis
#
# \arguments{
#   \item{this}{The Covariance object to clone.}
#   \item{...}{Not used.}
# }
#
# @author
#
# \value{A deep clone of the Covariance object.}
#
# \seealso{
#   @see "R.oo::clone.Object"
#   @seeclass
# }
#
#*/###########################################################################
clone.Covariance <- function(this, ...) {
  Cov <- clone.Object(this, ...)
  return (Cov)
}

############################################################################/**
# @RdocMethod EncodeForTraining
#
# @title "Makes parameters suitable for optimization"
#
# \description{
#   @get "title".
#   Specifically, we want "scale"-type parameters (i.e., where we're uncertain
#   about the magnitude) to be uncertain in logspace rather than regular-space.
#   This is equivalent to the Jeffreys prior.
# }
#
# @synopsis
#
# \arguments{
#   \item{this}{The Covariance object whose parameters to encode.}
#   \item{values}{A named vector of parameter values (with ID already
#      prepended).}
#   \item{...}{Not used.}
# }
#
# \value{
#   \code{values}, but with scale-type parameters changed to their log (and
#   appropriately renamed)
# }
#
# \seealso{
#   @seefunction DecodeForTraining
#   @seeclass
# }
#
# @author
#*/###########################################################################

setMethodS3("EncodeForTraining", "Covariance", conflict="quiet",
  function(this, values, ...) {
    # Takes a named (i.e., with the 'id' tag already prepended) vector of
    #    parameters (or lower or upper bounds); transforms "scale"-type
    #    parameters into logspace and change their names.
    #
    # Args:
    #   values:  A named vector of parameter values (with ID already
    #      prepended).
    #
    # Returns:
    #   'values', but with scale-type parameters changed to their log (and
    #   appropriately renamed)
    logspace.names <- this$PrependId(this$logspaceNames)
    i.log <- which(names(values) %in% logspace.names)
    values[i.log] <- log(values[i.log])
    names(values)[i.log] <- paste(sep="", names(values)[i.log], LogspaceTag())
    return (values)
  })

