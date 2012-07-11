############################################################################/**
# @RdocClass Covariance
#
# @title "Superclass for all Covariance function objects"
#
# \description{
#   \code{Covariance} is the superclass for more specific types (SE, Matern,
#   etc.).  Nobody will actually make an object of type \code{Covariance};
#   instead, they will use one of the derived classes (such as
#   \code{link{CovarianceSE}}).
#
#   @get "title".
#
#   @classhierarchy
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

#' The ID for this Covariance
#'
#' The ID is a character string identifying this Covariance object.
#'
#' @name getId.Covariance
#' @aliases Covariance$id getId.Covariance setId.Covariance
#' @S3method getId Covariance
#' @export getId getId.Covariance
#'
#' @param id A character string identifying this Covariance object.
#' @param ... Not used.
#'
#' @usage Covariance$id
#'
#' @return The id of this Covariance object.
#'
#' @seealso \code{\link{Covariance}}
setMethodS3("getId", "Covariance", conflict="quiet",
  function(this, ...) {
    return (this$.id)
  })
setMethodS3("setId", "Covariance", conflict="quiet",
  function(this, id, ...) {
    this$.id <- id
    return (this)
  })

#' ID-decorated names of Covariance parameters
#'
#' This function returns the names of the Covariance parameters, decorated with
#' the Covariance ID as a prefix.  This decoration helps distinguish the
#' parameters if a model has multiple copies of a given type of covariance.  For
#' example, if we used both coarse and fine lengthscales, we might have
#' parameters named "coarse.ell" and "fine.ell" instead of two parameters
#' ambiguously named "ell".
#'
#' @name getParamNames.Covariance
#' @aliases Covariance$paramNames getParamNames.Covariance
#' @S3method getParamNames Covariance
#' @export getParamNames getParamNames.Covariance
#'
#' @param ... Not used.
#'
#' @usage Covariance$paramNames
#'
#' @return Parameter names in the form "id.basename", where basename is the
#'    "undecorated" parameter name.
#'
#' @seealso \code{\link{getParams.Covariance}}
#' @seealso \code{\link{getParams.Model}}
#' @seealso \code{\link{Covariance}}
setMethodS3("getParamNames", "Covariance", conflict="quiet",
  function(this, ...) {
    if (!is.character(getId(this)) || nchar(getId(this)) < 1) {
      return (getParamNamesPlain(this))
    }
    return (gppois:::PrependId(this=this, strings=getParamNamesPlain(this)))
  })

#' Parameter values for this Covariance
#'
#' The current values of the parameters which govern this Covariance.
#'
#' \bold{Optimization mode}\cr
#' By default, we return all parameters.  However, this would if the caller
#' is an \emph{optimization} routine, there are at least two important
#' drawbacks.  First, if any parameters are fixed, it's wasteful (and
#' potentially hazardous) to pass these to the optimization routine.  Second,
#' and more seriously, it assumes a flat prior on \emph{all} parameters, even
#' "scale"-type parameters.  This causes the optimization to fail outright
#' even in many simple cases.
#'
#' To circumvent these problems, we provide an \emph{optimization mode},
#' which returns only non-constant parameters, and puts "scale"-type
#' parameters in logspace.  (The latter corresponds to the Jeffreys prior:
#' flat in \emph{log-space}, rather than real-space.  It represents
#' uncertainty about the \emph{order of magnitude} of the parameter.)
#'
#'
#' \bold{Handling crossed boundaries}\cr
#' The lower bound, upper bound, and value of every parameter must
#' \emph{always} be properly ordered:\cr
#' lower <= param <= upper\cr
#' Sometimes a proposed move might violate that condition.  The way it's
#' handled depends on whether it's a parameter that was moved, or a
#' boundary:\cr
#' \enumerate{
#'   \item Moved \emph{parameter}\cr
#'      The boundaries are hard.  The value will be clamped at the closest
#'      allowable value, and a warning will be given.  (Note that this implies a
#'      parameter can be set "constant" by setting its upper bound equal to its
#'      lower bound.)
#'   \item Moved \emph{boundary}\cr
#'       A moved boundary can "push" the parameter value, or even the other
#'       boundary.  e.g., if you have (lower=3, param=4, upper=5) and set upper
#'       to 2, the final state will be (lower=2, param=2, upper=2).
#' }
#'
#' @name getParams.Covariance
#' @aliases Covariance$params getParams.Covariance setParams.Covariance
#' @S3method setParams Covariance
#' @export setParams setParams.Covariance
#' @S3method getParams Covariance
#' @export getParams getParams.Covariance
#'
#' @param p (named numeric vector) New parameter values. (N.B: we \emph{only}
#'      use ones which are named, and whose names match up with names of
#'      parameters.)
#' @param for.training (logical) If TRUE, we return values more suitable for
#'      optimization: "scale"-type parameters are given as logarithms (which
#'      amounts to using the Jeffreys prior), and "constant" params (i.e., where
#'      the lower bound equals the upper bound) are omitted entirely.
#' @param ... Not used.
#'
#' @usage Covariance$params
#' @usage Covariance$params <- c(name1=value1, name2=value2, ...)
#'
#' @return (named numeric vector) Parameter values for this Covariance.
#'
#' @seealso \code{\link{getParamsPlain.Covariance}}
#' @seealso \code{\link{Covariance}}
setMethodS3("getParams", "Covariance", conflict="quiet",
  function(this, for.training=FALSE, ...) {
    p <- getParamsPlain(this)
    names(p) <- getParamNames(this)
    if (for.training) {
      p <- gppois:::EncodeForTraining(this, p)
    }
    return (p)
  })
setMethodS3("setParams", "Covariance", conflict="quiet",
  function(this, p, for.training=FALSE, ...) {
    if (for.training) {
      p <- gppois:::DecodeForTraining(p)
    }
    p.plain <- gppois:::UndecorateNames(this, p=p)
    this$setParamsPlain(p=p.plain)
    return (invisible(this))
  })


#' Set parameters using undecorated names
#'
#' Set parameters of this covariance using "undecorated" or "plain" names.
#'
#' \bold{Undecorated Names}\cr
#' Covariance objects decorate their parameter names with their ID: e.g., a
#' Covariance with id "SE" will have a parameter named "SE.ell".  This is
#' important when you have multiple Covariance objects; decorating with the ID
#' keeps the names unique.  However, it does make it less convenient to refer to
#' them.  This function allows to access them with the more intuitive,
#' \emph{undecorated} names.
#'
#' I do not provide anything like the \code{for.training} option from
#' \code{\link{getParams}}, because it doesn't make any sense.  The training
#' functions always have to use the \emph{decorated} names, or else risk name
#' collisions.
#'
#' \bold{Note for developers}\cr
#' For developers who plan to write new Covariance classes: the corresponding
#' \emph{accessor} method is the responsibility of subclasses.  i.e., you must
#' define a method "getParamsPlain" for your class.  (I suggest looking at the
#' corresponding methods for the \code{\link{CovarianceSE}} class to get
#' started.)
#'
#' @name setParamsPlain.Covariance
#' @aliases Covariance$paramsPlain setParamsPlain.Covariance getParamsPlain
#' @aliases getParamsPlain.Covariance paramsPlain
#' @S3method setParamsPlain Covariance
#' @export setParamsPlain setParamsPlain.Covariance
#'
#' @param p A (named) vector of new parameter values (we \emph{only} use ones
#'      which are named, and whose names match up with names of parameters.)
#' @param ... Not used.
#'
#' @usage Covariance$paramsPlain <- c(name1=value1, name2=value2, ...)
#'
#' @seealso \code{\link{getParamsPlain}}
#' @seealso \code{\link{getLowerPlain}}
#' @seealso \code{\link{setLowerPlain}}
#' @seealso \code{\link{Covariance}}
setMethodS3("setParamsPlain", "Covariance", conflict="quiet",
  function(this, p, ...) {
    # First, we need a vector where every parameter name has a value: either
    # the new value from 'p', or if none was given, the current value.
    p.good.names <- this$getParamsPlain()
    to.change <- names(p)[which(names(p) %in% names(p.good.names))]
    p.good.names[to.change] <- p[to.change]
    # paramsPlainImplementation requires a vector with a value for every
    # parameter; we took care of this above.
    gppois:::paramsPlainImplementation(this, p=p.good.names)
    p.clamped <- gppois:::ClampedParamVals(this, p=p.good.names)
    p.names <- names(p.clamped)
    clamped <- which(p.good.names[p.names] != p.clamped[p.names])
    gppois:::paramsPlainImplementation(this, p=p.clamped)
    if (any(clamped)) {
      culprits <- paste(sep='', collapse=' ', '"', p.names[clamped], '"')
      warning(paste("These parameters had to be clamped:\n", culprits, "\n"))
    }
    return (invisible(this))
  })

#' Lower bounds for parameters using full names
#'
#' Gets or sets the lower bounds for one or more parameters, using the full name
#' (ID + basename) of each parameter.  See section on Undecorated Names in help
#' page for \code{\link{paramsPlain}} for further commentary on name decoration.
#'
#' @name setLower.Covariance
#' @aliases Covariance$lower setLower.Covariance getLower getLower.Covariance
#' @S3method setLower Covariance
#' @export setLower setLower.Covariance
#' @S3method getLower Covariance
#' @export getLower getLower.Covariance
#'
#' @param L A (named) vector of new lower bounds on parameter values (we
#'    \emph{only} use ones which are named, and whose names match up with names
#'    of parameters.)
#' @param for.training (logical) If TRUE, we return values more suitable for
#'      optimization: "scale"-type parameters are given as logarithms (which
#'      amounts to using the Jeffreys prior), and "constant" params (i.e., where
#'      the lower bound equals the upper bound) are omitted entirely.
#' @param ... Not used.
#'
#' @usage Covariance$lower <- c(name1=value1, name2=value2, ...)
#'
#' @seealso \code{\link{getParamsPlain}}
#' @seealso \code{\link{getLowerPlain}}
#' @seealso \code{\link{setLowerPlain}}
#' @seealso \code{\link{Covariance}}
setMethodS3("getLower", "Covariance", conflict="quiet",
  function(this, for.training=FALSE, ...) {
    L <- this$getLowerPlain()
    names(L) <- this$getParamNames()
    if (for.training) {
      L <- gppois:::EncodeForTraining(this, L)
    }
    return (L)
})
setMethodS3("setLower", "Covariance", conflict="quiet",
  function(this, L, for.training=FALSE, ...) {
    if (for.training) {
      L <- gppois:::DecodeForTraining(L)
    }
    L.plain <- gppois:::UndecorateNames(this, L)
    this$setLowerPlain(L=L.plain)
    return (invisible(this))
})

#' Upper bounds for parameters using full names
#'
#' Gets or sets the upper bounds for one or more parameters, using the full name
#' (ID + basename) of each parameter.  See section on Undecorated Names in help
#' page for \code{\link{paramsPlain}} for further commentary on name decoration.
#'
#' @name setUpper.Covariance
#' @aliases Covariance$upper setUpper.Covariance getUpper getUpper.Covariance
#' @S3method setUpper Covariance
#' @export setUpper setUpper.Covariance
#' @S3method getUpper Covariance
#' @export getUpper getUpper.Covariance
#'
#' @param L A (named) vector of new upper bounds on parameter values (we
#'    \emph{only} use ones which are named, and whose names match up with names
#'    of parameters.)
#' @param for.training (logical) If TRUE, we return values more suitable for
#'      optimization: "scale"-type parameters are given as logarithms (which
#'      amounts to using the Jeffreys prior), and "constant" params (i.e., where
#'      the upper bound equals the upper bound) are omitted entirely.
#' @param ... Not used.
#'
#' @usage Covariance$upper <- c(name1=value1, name2=value2, ...)
#'
#' @seealso \code{\link{getParamsPlain}}
#' @seealso \code{\link{getUpperPlain}}
#' @seealso \code{\link{setUpperPlain}}
#' @seealso \code{\link{Covariance}}
setMethodS3("getUpper", "Covariance", conflict="quiet",
  function(this, for.training=FALSE, ...) {
    U <- this$getUpperPlain()
    names(U) <- this$getParamNames()
    if (for.training) {
      U <- gppois:::EncodeForTraining(this, U)
    }
    return (U)
  })
setMethodS3("setUpper", "Covariance", conflict="quiet",
  function(this, U, for.training=FALSE, ...) {
    if (for.training) {
      U <- gppois:::DecodeForTraining(U)
    }
    U.plain <- gppois:::UndecorateNames(this, U)
    this$setUpperPlain(U=U.plain)
    return (invisible(this))
  })

#' Deep-clone a Dataset
#'
#' Deep-clone a Dataset (i.e., clones LazyMatrix objects too).
#'
#' @method clone Covariance
#'
#' @param this The Covariance object to clone.
#' @param ... Not used.
#'
#' @return A deep clone of the Covariance object.
#'
#' @export
#' @seealso \code{\link{Covariance}}
clone.Covariance <- function(this, ...) {
  Cov <- clone.Object(this, ...)
  return (Cov)
}

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
    logspace.names <- gppois:::PrependId(this, this$logspaceNames)
    i.log <- which(names(values) %in% logspace.names)
    values[i.log] <- log(values[i.log])
    names(values)[i.log] <- paste(sep="", names(values)[i.log], LogspaceTag())
    return (values)
  })

#' Set a parameter to a constant value.
#'
#' Sets the lower and upper bounds equal for a given parameter, effectively
#' fixing that parameter as constant.
#'
#' @S3method FixConstParam Covariance
#' @export FixConstParam FixConstParam.Covariance
#' @name FixConstParam.Covariance
#'
#' @param p.name  The (undecorated) name of the parameter to change.
#' @param p.value  The new constant value for the named parameter.
#' @param ... Not used.
#'
#' @seealso \code{\link{getParams.Covariance}}
#' @seealso \code{\link{Covariance}}
setMethodS3("FixConstParam", "Covariance", conflict="quiet",
  function(this, p.name, p.value, ...) {
    p.old <- this$paramsPlain[p.name]
    p <- p.value
    names(p) <- p.name
    # Make the parameter into a constant:
    this$upperPlain <- p.old
    this$lowerPlain <- p.old
    # The order in which we move the values is chosen to preserve the ordering,
    # (upper>param>lower), and thereby avoid generating warnings:
    if (p.value > p.old) {  # Going up
      this$upperPlain <- p
      this$paramsPlain <- p
      this$lowerPlain <- p
    }
    if (p.value < p.old) {  # Going down
      this$lowerPlain <- p
      this$paramsPlain <- p
      this$upperPlain <- p
    }
    return (invisible(this))
  })

#' Covariance matrix
#'
#' Gives the covariance matrix K between an input X1 and output X2, where 
#' K[i, j] is the covariance function evaluated at datapoints X1[i] and X2[j]:
#'
#' \describe{
#'    \item{\code{KInIn}:}{X1 = d$X, X2 = d$X}
#'    \item{\code{KInOut}:}{X1 = X.out, X2 = d$X}
#'    \item{\code{KOutIn}:}{X1 = d$X, X2 = X.out}
#'    \item{\code{KOutOut}:}{X1 = X.out, X2 = X.out}
#' }
#'
#' @name K.Covariance
#' @aliases KInIn KInOut KOutIn KOutOut
#' @S3method KInIn Covariance
#' @export KInIn KInIn.Covariance
#' @S3method KInOut Covariance
#' @export KInOut KInOut.Covariance
#' @S3method KOutIn Covariance
#' @export KOutIn KOutIn.Covariance
#' @S3method KOutOut Covariance
#' @export KOutOut KOutOut.Covariance
#'
#' @param d  A Dataset object encapsulating the data to train on.
#' @param X.out  A matrix (with d$d columns) of X-locations where we want
#'      predictions.
#' @param ... Not used.
#'
#' @return A \code{ncol(X2) x ncol(X1)} matrix: the covariance from each input
#'    point to each output point.
#'
#' @seealso \code{\link{KInInDeriv}}
#' @seealso \code{\link{Covariance}}
setMethodS3("KInIn", "Covariance", conflict="quiet",
  function(this, d, ...) {
    return (this$K.specific(X=d$X))
  })
setMethodS3("KInOut", "Covariance", conflict="quiet",
  function(this, d, X.out=d$X, ...) {
    return (this$K.specific(X=d$X, X.out=X.out))
  })

setMethodS3("KOutIn", "Covariance", conflict="quiet",
  function(this, d, X.out=d$X, ...) {
    return (t(this$KInOut(d=d, X.out=X.out)))
  })
setMethodS3("KOutOut", "Covariance", conflict="quiet",
  function(this, X.out, ...) {
    return (this$K.specific(X=X.out))
  })

#' Element-wise derivatives of covariance matrix
#'
#' Computes element-wise derivatives of the "In-In" covariance matrix (i.e., the
#' covariance matrix \emph{from} the observed datapoints \code{d$X}, \emph{to}
#' those same datapoints \code{d$X}).  These derivatives are the only ones
#' important for computing derivatives of the posterior probability function.
#'
#' @S3method KInInDeriv Covariance
#' @export KInInDeriv KInInDeriv.Covariance
#' @name KInInDeriv.Covariance
#' @aliases KInInDeriv KInInDeriv.Covariance
#'
#' @param d  The Dataset whose X-values we are training on.
#' @param param  The name of the parameter with respect to which we're
#'      differentiating.
#' @param ... Not used.
#'
#' @note
#'   As this is the superclass, all we do here is to check whether this is
#'   our problem (i.e., whether this object actually has a parameter named
#'   \code{param}).  If not, just return a zero matrix.  Otherwise, we let the
#'   subclasses handle the actual computation.
#'
#' @return The element-by-element derivative of the total K-matrix for Dataset
#'    \code{d}, with respect to the parameter named \code{param}.
#'
#' @seealso \code{\link{K.Covariance}}
#' @seealso \code{\link{Covariance}}
setMethodS3("KInInDeriv", "Covariance", conflict="quiet",
  function(this, d, param, ...) {
    if (param %in% this$paramNames) {
      names(param) <- param
      param.plain <- names(gppois:::UndecorateNames(this, p=param))
      K.deriv <- this$KDerivImplementation(d=d, param=param.plain)
    } else {
      K.deriv <- matrix(0, nrow=d$n, ncol=d$n)
    }
    return (K.deriv)
  })

setMethodS3("PrependId", "Covariance", private=TRUE, conflict="quiet",
  function(this, strings, ...) {
    return (paste(sep='.', this$getId(), strings))
  })

#' Pretty-printing for Covariance objects
#'
#' Prints out the id of the Covariance, followed by a list of its parameters.
#' For each parameter, we print the name, lower bound, current value, and upper
#' bound.
#'
#' @method print Covariance
#'
#' @param this The Covariance object to print.
#' @param indent Aids in formatting: the number of spaces to print before every
#'    line.
#' @param ... Not used.
#'
#' @export
#' @seealso \code{\link{Covariance}}
print.Covariance <- function(this, indent=0, ...) {
  tab <- paste(collapse='', rep(' ', indent))
  cat(sprintf("%sCovariance, id=%-20s (%s)\n", tab, Wrap(this$id, "'"),
      class(this)[1]))
  PrintParams(lower=this$lowerPlain, params=this$paramsPlain,
    upper=this$upperPlain, indent=indent, ...)
  return (invisible(this))
}

LogspaceTag <- function() {
  # A tag we append to variable names to indicate they have been transformed
  # into logspace.
  #
  # Returns:
  #   A string: this tag.
  return (".__LOG__")
}

DecodeForTraining <- function(values, ...) {
  # Takes a named vector of parameters (or lower or upper bounds); transforms
  #    any tagged as "logspace" back to linear (changing names
  #    appropriately).
  #
  # Args:
  #   values:  A named vector of parameter values
  #
  # Returns:
  #   'values', but with log-transformed parameters changed back to linear
  #   scale (and appropriately renamed)
  log.pattern <- paste(sep='', LogspaceTag(), '$')
  i.log <- grep(pattern=log.pattern, x=names(values))
  values[i.log] <- exp(values[i.log])
  names(values)[i.log] <- sub(pattern=log.pattern, replacement='',
    x=names(values)[i.log])
  return (values)
}

setMethodS3("ClampedParamVals", "Covariance", private=TRUE, conflict="quiet",
  # Tells what the current parameter values *would* be, if clamped to lie
  # within the lower and upper bounds.  (Naturally, this is intended for
  # *internal* use, by every function which allows the parameter values to
  # change: they are expected to use *this* function so the parameter values
  # *always* stay within bounds.)
  #
  # Args:
  #   warn: Whether or not to print a warning in case params are outside
  #     bounds.
  #
  # Returns:
  #   A named vector of the new parameter values (with plain names).
  function(this, warn=TRUE, p=this$getParamsPlain(), ...) {
    lower <- this$getLowerPlain()
    upper <- this$getUpperPlain()
    p.new <- ClampNamed(x=p, lower=lower, upper=upper)
    if (warn && !isTRUE(all.equal(p.new[names(p)], p[names(p)]))) {
      culprits <- paste(sep='', collapse=' ', '"',
        names(p)[which(p != p.new)], '"')
      warning(paste("These parameters were outside the bounds:\n",
          culprits, "\n"))
    }
    return (p.new[names(p)])
  })

setMethodS3("ClampParams", "Covariance", private=TRUE, conflict="quiet",
  function(this, warn=TRUE, ...) {
    # Clamps the parameter values to fall between the upper and lower
    # boundaries; also, invalidates K if they change and optionally prints a
    # warning.
    #
    # Args:
    #   warn: If TRUE, we print a warning if anything actually changes.
    #
    # Returns:
    #   Object whose parameters got clamped, but used for its side-effects.

    # The following *should* be OK, as long as neither setParams nor
    # setParamsPlain calls ClampParams...
    this$setParamsPlain(p=gppois:::ClampedParamVals(this, warn=warn))
    return (invisible(this))
  })

setMethodS3("UndecorateNames", "Covariance", private=TRUE, conflict="quiet",
  function(this, p, ...) {
    # Returns entries in 'p' which match 'this$id', but with the id stripped
    # out of the names vector.
    #
    # Args:
    #   p: A (named) vector of values corresponding to the parameters of this
    #      Covariance (they could be the parameter values themselves, or
    #      perhaps upper or lower bounds).
    #
    # Returns:
    #   'p' stripped down so it only has the relevant entries with plain names.

    # Which entries in 'p' correspond to our parameters?
    relevant.indices <- which(names(p) %in% this$getParamNames())
    p.relevant <- p[relevant.indices]
    # Strip out our ID string
    id.prefix <- paste(sep='', this$getId(), '.')
    names(p.relevant) <- gsub(id.prefix, '', names(p.relevant))
    return (p.relevant)
  })

setMethodS3("PushUpperBounds", "Covariance", private=TRUE, conflict="quiet",
  function(this, U.min, warn.on.cross=TRUE, ...) {
    # Adjusts upper boundaries to "make way" for the new lower bounds in
    # 'U.min', giving warnings if any adjustments are actually necessary.
    #
    # Args:
    #   U.min:  New lower bounds for corresponding parameters of this Covariance,
    #      where names(U.min) are *already* undecorated.
    #   warn.on.cross:  If TRUE, we will warn if the bounds are actually moved
    #
    # Returns:
    #   Plain-named version of 'U.min' with only the relevant entries (though
    #   this method is mainly used for its side-effect).

    # Which ones will be affected?  (If none, just skip to the finish.)
    to.change <- names(U.min)[which(names(U.min) %in% this$paramNamesPlain)]
    if (length(to.change) < 1) {
      return (c())
    }
    # Remember old values so we can check later whether they've changed
    upper.old <- this$upperPlain
    # If the bounds cross, set them both to the new value
    bounds.crossed <- U.min[to.change] > upper.old[to.change]
    if (any(bounds.crossed)) {
      this$setUpperPlain(U=pmax(upper.old[to.change], U.min[to.change]))
      if (warn.on.cross) {
        culprits <- paste(sep='', collapse=' ', '"',
          to.change[which(bounds.crossed)], '"')
        warning(paste("Upper bounds forcibly moved for these parameters:\n",
            culprits, "\n"))
      }
    }
    return (U.min[to.change])
  })

setMethodS3("PushLowerBounds", "Covariance", private=TRUE, conflict="quiet",
  function(this, L.max, warn.on.cross=TRUE, ...) {
    # Adjusts lower boundaries to "make way" for the new upper bounds in 'L.max',
    # giving warnings if any adjustments are actually necessary.
    #
    # Args:
    #   L.max: New upper bounds for corresponding parameters of this Covariance,
    #      where names(L.max) are *already* undecorated.
    #   warn.on.cross:  If TRUE, we will warn if the bounds are actually moved
    #
    # Returns:
    #   Plain-named version of 'L.max' with only the relevant entries (though
    #   this method is mainly used for its side-effect).

    # Which ones will be affected?  (If none, just skip to the finish.)
    to.change <- names(L.max)[which(names(L.max) %in% this$paramNamesPlain)]
    if (length(to.change) < 1) {
      return (c())
    }
    # Remember old values so we can check later whether they've changed
    lower.old <- this$lowerPlain
    # If the bounds cross, set them both to the new value
    bounds.crossed <- L.max[to.change] < lower.old[to.change]
    if (any(bounds.crossed)) {
      this$setLowerPlain(L=pmin(lower.old[to.change], L.max[to.change]))
      if (warn.on.cross) {
        culprits <- paste(sep='', collapse=' ', '"',
          to.change[which(bounds.crossed)], '"')
        warning(paste("Lower bounds forcibly moved for these parameters:\n",
            culprits, "\n"))
      }
    }
    return (L.max[to.change])
  })

