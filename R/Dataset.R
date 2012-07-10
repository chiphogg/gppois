#############################################################################/**
# @RdocClass Dataset
#
# @title "A wrapper class for data being analyzed"
#
# \description{
#   This is the \emph{constructor} for a \code{Dataset} object:
#   @get "title".
#
#   Here is the class hierarchy: 
#   @classhierarchy
#
# }
#
# @synopsis
#
# \arguments{
#   \item{id}{(character) An id which identifies this Dataset.}
#   \item{data}{(\code{data.frame}) The raw data.}
#   \item{X.names}{(character vector) The names of the columns which represent
#      covariates.  (NB: the length of this vector implicitly defines the
#      dimensionality of the data: if you give 2 column names, you have a 2D
#      dataset.)}
#   \item{noise.var}{(named numeric vector) Desired noise levels for each
#      column (i.e., a threshhold for precision: we don't care about
#      differences smaller than the corresponding scale).}
#   \item{column}{(character) The name of the quantity to select for analysis.}
#   \item{poisson.names}{(character vector) The column names which represent
#      Poisson data.}
#   \item{tol.factor}{(numeric) The relative tolerance factor: by default, we
#      set the noise variance for each column to \code{tol.factor * sd(col)}.}
#   \item{data.offset}{(named numeric vector) An optional offset which is
#      subtracted from each column.  If an unnamed value is included, it is
#      treated as a default.  The "default default" is to subtract off the mean
#      for each column.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods
#
# }
#
# @author
#*/###########################################################################
setConstructorS3("Dataset",
  function(id="UNNAMED", data=data.frame(), X.names="X", noise.var=c(), column="",
    poisson.names=c(), tol.factor=1e-4, data.offset=c(), ...) {
    # Ignore "factor"-type variables
    data.trim <- data
    factor.i <- which(sapply(data, class) == 'factor')
    if (length(factor.i) > 0) {
      data.trim <- data[, -factor.i]
    }

    # Only choose the X columns which actually exist, and designate the rest as
    # candidate "data" columns.
    existing.X.names <- X.names[which(X.names %in% colnames(data.trim))]
    other.names <- colnames(data.trim)[
      -which(colnames(data.trim) %in% existing.X.names)]

    # Pick a default column for analysis.
    column <- column[1]
    if (length(other.names) > 0 && !any(column %in% other.names)) {
      column <- other.names[1]
    }

    # Compute the noise variance.
    n.var <- ComputeNoiseVariance(data=data.trim, specs=noise.var, tol=tol.factor,
      poisson.names=poisson.names)

    # Compute data offsets, according to the following logic:
    # 1. Any named entries get set to the given value.
    # 2. The first non-named entry (if any) is a default value.
    # 3. Unassigned entries get set to the default; if none exists, they are
    #    set to the mean of the data.
    if (length(data.offset) > 0 && is.null(names(data.offset))) {
      names(data.offset) <- rep('', length(data.offset))
    }
    unnamed <- which(nchar(names(data.offset)) == 0)
    default.offset <- ifelse(length(unnamed) > 0, data.offset[unnamed[1]], NA)
    offsets <- rep(default.offset, length(other.names))
    names(offsets) <- other.names
    # Add specifically-named values
    for (n in names(data.offset)[-unnamed]) {
      offsets[n] <- data.offset[n]
    }
    # Change NA to mean
    for (i in which(is.na(offsets))) {
      n <- names(offsets)[i]
      is.poisson <- n %in% poisson.names
      dpts.gauss <- dpts <- data.trim[, n]
      if (is.poisson) {
        dpts.gauss <- Anscombe(dpts)
      }
      offsets[n] <- mean(dpts.gauss)
    }

    extend(Object(), "Dataset",
      .data.frame      = data.trim,
      .data.full       = data,
      .data.offset     = offsets,
      .id              = id,
      # The computed noise variance for each column, based on the above specs:
      .noise.var       = n.var,
      # The user's original specifications about the noise variance (useful to
      # have, in case we delete rows and want to reconstruct noiseVar):
      .noise.var.specs = noise.var,
      # Columns which are Poisson distributed.
      .poisson         = colnames(data)[which(colnames(data) %in% poisson.names)],
      # The column we're actually using for predictions:
      .quantity        = column,
      .tol.factor      = tol.factor,
      .X.names         = existing.X.names
      )
  })

#-------------------------------------------------------------------------------
# (Dataset) "HELPER" FUNCTIONS:

ComputeNoiseVariance <- function(data, specs, tol, poisson.names=c()) {
  # Compute the noise variance for each column in 'data', using any custom
  # values in 'specs', or the relative tolerance 'tol' for entries lacking such
  # custom values.
  #
  # Args:
  #   data:  A data.frame holding the columns to specify.
  #   specs:  A named vector of values giving custom noise variances for
  #      correspondingly-named columns in 'data'.
  #   tol:  The default amount of noise we can tolerate in a column, relative
  #      to the sample standard deviation of the data in that column.
  #   poisson.names:  A list of columns which have the Poisson distribution;
  #      *these* columns should have a default value of 0.25, based on the
  #      Anscombe transform.
  #
  # Returns:
  #   A vector with one named entry per column of 'data', giving the Gaussian
  #   variance of the noise in this column.
  #
  # RULES FOR NOISE VARIANCE:
  # 1) Any column value manually specified takes top precedence
  # 2) Otherwise, if a column is Poisson, it gets set to 0.25
  # 3) Un-named values in 'specs':
  #  a) Any unnamed values get doled out as defaults
  #  b) If there IS NO default value, set unspecified to (tol * sd), to
  #     aid numerical stability
  #
  # NOTE that it is computationally advantageous to have nonzero noise in
  # almost every case.  Otherwise, K-matrices tend to be very close to
  # singular.  That's the idea behind specifying a relative tolerance: it makes
  # a relatively painless way to say you don't care about fluctuations smaller
  # than 'tol' (in units of your data's fluctuations, which we guess to be in
  # the same ballpark as the sample standard deviation).  Really though, users
  # should look at their data and manually decide what level fluctuations they
  # could live with for each column.

  # Setup the object
  noise.var <- rep(NA, ncol(data))
  names(noise.var) <- colnames(data)

  # Rule 1 implemented:
  given.names <- names(specs)[which(names(specs) %in% names(noise.var))]
  noise.var[given.names] <- specs[given.names]

  # Rule 2 implemented:
  is.Poisson <- names(noise.var) %in% poisson.names
  noise.var[which(is.na(noise.var) & is.Poisson)] <- 0.25

  # Rule 3 implemented:
  # 'leftover': names of columns which don't have a noise variance yet
  leftover <- names(noise.var)[which(is.na(noise.var))]
  if (length(leftover) > 0) {
    # 'defaults': initialize with the "Rule 3b" values
    defaults <- (tol * sd(data[, leftover])) ^ 2
    # If specs has any named entries, it's safe to look for unnamed ones
    # (i.e., Rule 3a)
    if (!is.null(names(specs)) ) {
      supplied.defaults <- specs[which(is.na(names(specs)) | names(specs)=="")]
      if (length(supplied.defaults) > 1) default <- supplied.defaults
    }
    noise.var[leftover] <- defaults
  }
  return (noise.var)
}

DatasetFromFile <- function(file, sep='\t', id=basename(file), ...) {
  # Read a file's contents into a data.frame; then, use that data.frame to
  # construct a Dataset object.
  #
  # Args:
  #   file:  A character string holding the name of the file to read.
  #   sep:  The character string which separates entries on a line.
  #
  # Returns:
  #   A Dataset object holding the data in the given file
  data.from.file <- read.table(header=TRUE, file=file, sep=sep)
  # Default ID is the basename of the file
  return (Dataset(data=data.from.file, id=id, ...))
}

#' The dimensionality of the Dataset
#'
#' The dimensionality of the independent variable in this Dataset. i.e., 1 for
#' curves, 2 for surfaces, etc.
#'
#' @S3method getD Dataset
#'
#' @name getD.Dataset
#' @aliases getD d getD.Dataset
#' @S3method getD Dataset
#'
#' @param ... Not used.
#'
#' @usage Dataset$d
#'
#' @return An integer representing the dimensionality of the independent
#'   variable.  (i.e., 1 for curves, 2 for surfaces, etc.)
#'
#' @seealso \code{\link{Dataset}}
setMethodS3("getD", "Dataset", conflict="quiet",
  function(this, ...) {
    return (length(this$.X.names))
  })

#' The offset for this dataset
#'
#' The offset of the currently-selected quantity in the dataset.
#'
#' Usually accesed using \code{this$d} syntax; read-only.
#'
#' @name getDataOffset.Dataset
#' @aliases dataOffset getDataOffset getDataOffset.Dataset
#' @S3method getDataOffset Dataset
#'
#' @param ... Not used.
#'
#' @usage Dataset$d
#'
#' @return The datapoint offset values for the currently-selected quantity.
#'
#' @seealso \code{\link{Dataset}}
setMethodS3("getDataOffset", "Dataset", conflict="quiet",
  function(this, ...) {
    return (this$.data.offset[this$getQuantity()])
  })


#' Datapoint values
#'
#' Values for the currently-selected quantity.
#'
#' Usually accesed using \code{this$dpts} syntax; read-only.
#'
#' @name getDpts.Dataset
#' @aliases dpts getDpts getDpts.Dataset
#' @S3method getDpts Dataset
#'
#' @param ... Not used.
#'
#' @usage Dataset$dpts
#'
#' @return The raw datapoint values for the currently-selected quantity.
#'
#' @seealso \code{\link{Dataset}}
#' @seealso \code{\link{xformedDpts}}
setMethodS3("getDpts", "Dataset", conflict="quiet",
  function(this, ...) {
    return (this$.data.frame[, this$getQuantity()])
  })

#' Transformed datapoint values
#'
#' Transformed values for the currently-selected quantity.  Transformations
#' include the Anscombe transform (for Poisson data), followed by subtracting
#' off the mean (if desired).
#'
#' Usually accesed using \code{this$xformedDpts} syntax; read-only.
#'
#' @name getXformedDpts.Dataset
#' @aliases xformedDpts getXformedDpts getXformedDpts.Dataset
#' @S3method getXformedDpts Dataset
#'
#' @param ... Not used.
#'
#' @usage Dataset$xformedDpts
#'
#' @return The datapoint values for the currently-selected quantity, transformed
#'    such that they can be assumed to have a Gaussian distribution.
#'
#' @seealso \code{\link{Dataset}}
#' @seealso \code{\link{dpts}}
setMethodS3("getXformedDpts", "Dataset", conflict="quiet",
  function(this, ...) {
    gaussian.dpts <- dpts <- this$getDpts()
    if (this$getIsPoisson()) {
      # Anscombe transform: turns Poisson data into Gaussian data with constant
      # variance of 1/4.
      gaussian.dpts <- Anscombe(dpts)
    }
    return (gaussian.dpts - this$getDataOffset())
  })

#' The id string for this Dataset
#'
#' Each Dataset has an id string to give a hint about what it is.  These are not
#' checked for uniqueness.
#'
#' Usually accesed using \code{this$id} syntax.
#'
#' @name Dataset$id
#' @aliases id getId setId getId.Dataset setId.Dataset
#' @S3method getId Dataset
#'
#' @param ... Not used.
#'
#' @usage Dataset$id
#' @usage Dataset$id <- "example_id"
#'
#' @return A string representing the ID for this Dataset.
#'
#' @seealso \code{\link{Dataset}}
setMethodS3("getId", "Dataset", conflict="quiet",
  function(this, ...) {
    return (this$.id)
  })
setMethodS3("setId", "Dataset", conflict="quiet",
  function(this, id, ...) {
    this$.id <- id
    return (this)
  })

#' Is the current column Poisson?
#'
#' Checks whether the currently-selected column consists of Poisson-noised data.
#'
#' Usually accesed using \code{this$isPoisson} syntax; read-only.
#'
#' @name getIsPoisson.Dataset
#' @aliases isPoisson getIsPoisson getIsPoisson.Dataset
#' @S3method getIsPoisson Dataset
#'
#' @param ... Not used.
#'
#' @usage Dataset$isPoisson
#'
#' @return (logical) TRUE if the raw data for the currently selected column is
#'    Poisson.
#'
#' @seealso \code{\link{Dataset}}
setMethodS3("getIsPoisson", "Dataset", conflict="quiet",
  function(this, ...) {
    return (any(this$.quantity %in% this$.poisson))
  })

#' The number of datapoints
#'
#' Gives the number of datapoints in this Dataset object.
#'
#' Usually accesed using \code{this$n} syntax; read-only.
#'
#' @name getN.Dataset
#' @aliases n getN getN.Dataset
#' @S3method getN Dataset
#'
#' @param ... Not used.
#'
#' @usage Dataset$n
#'
#' @return (numeric) The number of datapoints in this Dataset.
#'
#' @seealso \code{\link{Dataset}}
setMethodS3("getN", "Dataset", conflict="quiet",
  function(this, ...) {
    return (nrow(this$.data.frame))
  })

#' Noise variance for current column
#'
#' The noise variance for the quantity in this Dataset which is currently
#' selected.
#'
#' Usually accesed using \code{this$noiseVar} syntax; read-only.
#'
#' @name getNoiseVar.Dataset
#' @aliases noiseVar getNoiseVar getNoiseVar.Dataset
#' @S3method getNoiseVar Dataset
#'
#' @param ... Not used.
#'
#' @usage Dataset$noiseVar
#'
#' @return (numeric) The noise variance for the currently selected column.
#'
#' @seealso \code{\link{Dataset}}
#' @seealso \code{\link{quantity}}
setMethodS3("getNoiseVar", "Dataset", conflict="quiet",
  function(this, ...) {
    return (this$.noise.var[this$getQuantity()])
  })

#' Currently selected quantity (column)
#'
#' A Dataset consists of \code{d} covariates and one or more observed variables.
#' Only one of these variables can be analyzed at a time; \code{quantity} gives
#' the name of the variable currently selected.
#'
#' Usually accesed using \code{this$quantity} syntax.
#'
#' @name Dataset$quantity
#' @aliases quantity getQuantity setQuantity getQuantity.Dataset setQuantity.Dataset
#' @S3method getQuantity Dataset
#' @S3method setQuantity Dataset
#'
#' @param ... Not used.
#'
#' @usage Dataset$quantity
#' @usage Dataset$quantity <- new.quantity
#'
#' @return (character) The name of the quantity which is currently selected.
#'
#' @seealso \code{\link{Dataset}}
setMethodS3("getQuantity", "Dataset", conflict="quiet",
  function(this, ...) {
    return (this$.quantity)
  })
setMethodS3("setQuantity", "Dataset", conflict="quiet",
  function(this, col.name, ...) {
    if (col.name %in% colnames(this$.data.frame)) {
      this$.quantity <- col.name
    } else {
      warning(sprintf("You have selected a column ('%s') which does not exist! Column NOT changed.\n", col.name))
    }
    return (invisible(this))
  })

#' X-values
#'
#' Gives the X-values of the Dataset: i.e., the covariates where we have
#' observations.
#'
#' Usually accesed using \code{this$x} syntax; read-only.
#'
#' @name getX.Dataset
#' @aliases x getX getX.Dataset
#' @S3method getX Dataset
#'
#' @param ... Not used.
#'
#' @usage Dataset$x
#'
#' @return (\code{n} x \code{d}) matrix giving the X-values for this Dataset.
#'
#' @seealso \code{\link{Dataset}}
setMethodS3("getX", "Dataset", conflict="quiet",
  function(this, ...) {
    return (as.matrix(this$.data.frame[, this$.X.names]))
  })

#' Delete datapoints
#'
#' Removes the indicated datapoints from the Dataset.
#'
#' @S3method DeleteRows Dataset
#' @name DeleteRows.Dataset
#'
#' @param indices The indices of the datapoints to be removed.
#' @param ... Not used.
#'
#' @usage Dataset$DeleteRows(indices)
#'
#' @return \code{this}, but used for its side-effect of deleting rows.
#'
#' @seealso \code{\link{Dataset}}
setMethodS3("DeleteRows", "Dataset", conflict="quiet",
  function(this, indices, ...) {
    this$.data.frame <- this$.data.frame[-indices, ]
    this$.noise.var <- ComputeNoiseVariance(data=this$.data.frame,
      specs=this$.noise.var.specs, tol=this$.tol.factor,
      poisson.names=this$.poisson)
    return (invisible(this))
  })

#' Mean square residuals
#'
#' Mean squared residuals are taken with respect to the currently selected
#' column by default.
#'
#' @S3method MSR Dataset
#' @name MSR.Dataset
#'
#' @param test.data A sequence of datapoints, which live in the same space as
#'      xformedDpts (i.e., same offset, and Anscombe-transformed if Poisson).
#' @param reference.data The comparison function: defaults to our noisy
#'      datapoints, but could also be different (e.g., if we have access to the
#'      true function).
#' @param ... Not used.
#'
#' @return The mean squared difference between \code{test.data} and
#'    \code{reference.data}.
#'
#' @seealso \code{\link{Dataset}}
setMethodS3("MSR", "Dataset", conflict="quiet",
  function(this, test.data, reference.data=this$xformedDpts, ...) {
    if (length(test.data) != length(reference.data)) {
      stop("Trying to find MSR for sequences of different lengths")
    }
    return (mean((test.data - reference.data) ^ 2))
  })

#' Scatterplot for 2D dataset
#'
#' Use rgl to plot 2D data.  If there are too many datapoints, plotting and
#' interacting become very slow.  Thus, this function allows to limit the number
#' plotted using the \code{max.points} parameter.
#'
#' @S3method Plot2D Dataset
#' @name Plot2D.Dataset
#'
#' @param max.points The maximum number of points to show (these are randomly
#'    sampled from the available points).
#' @param dist.factor The ratio of the datapoint radius to the minimum
#'    datapoint separation.
#' @param ... Other graphical parameters for rgl functions.
#'
#' @return Used for its side-effect.
#'
#' @seealso \code{\link{Dataset}}
setMethodS3("Plot2D", "Dataset", conflict="quiet",
  function(this, max.points=1000, dist.factor=0.2, clear=TRUE, ...) {
    d <- clone(this)
    if (max.points < this$n) {
      i <- sample(x=1:this$n, size=max.points, replace=FALSE)
      d$DeleteRows(-i)
    }
    unit <- min(dist(d$X))
    if (require("rgl") == FALSE) {
      stop("The Plot2D() method requires the rgl library to be installed.")
    }
    if (clear) {
      rgl.clear()
    }
    rgl.spheres(x=d$X[, 1], z=d$X[, 2], y=d$dpts, radius=dist.factor * unit, ...)
    Y.scale <- ifelse(hasArg(Y.scale), list(...)$Y.scale,
      max(dist(d$X)) / diff(range(this$dpts)))
    ScaleY(Y.scale=Y.scale)
    return (invisible(this))
  })

#' Pretty-printing for Dataset objects
#'
#' Prints out the id of the Dataset, the number of datapoints, and the number of
#' dimensions.  Then, it lists the covariates (name and range), followed by the
#' various quantities which can be selected (name, type, range, and noise
#' variance).
#'
#' @method print Dataset
#'
#' @param this The Dataset object to print.
#' @param ... Not used.
#'
#' @export
#' @seealso \code{\link{Dataset}}
print.Dataset <- function(this, ...) {
  cat(sprintf("Dataset '%s': %d datapoints in %d dimensions.\n",
      this$id, this$n, this$d))
  # Ruler:
  #     1...5...10...15...20...25...30...35...40...45...50...55...60...65...70
  s <- "            NAME:            RANGE:"
  s.quantities <- paste(sep='', s,       "      NOISE VARIANCE:\n")
  s.covariates <- paste(sep='', s,       "\n")
  format.str <- "%17s%4s%14.8g%21s\n"
  raw.data <- this$.data.frame
  quantity.names <- colnames(raw.data)
  decoration <- "------------------"
  cat(sprintf("%16s%s%s%s\n", '', decoration, 'Covariates:', decoration))
  cat(s.covariates)
  X.indices <- which(quantity.names %in% this$.X.names)
  for (X.i in X.indices) {
    cat(sprintf(format.str, quantity.names[X.i], '',
        diff(range(raw.data[, X.i])), ''))
  }
  cat(sprintf("\n%16s%s%s%s\n%18s%s\n", '', decoration, 'Quantities:',
      decoration, '', "('*' = currently selected; 'P' = Poisson)"))
  cat(s.quantities)
  for (X.i in (1:ncol(raw.data))[-X.indices]) {
    sel <- ifelse(quantity.names[X.i] == this$quantity, '*', ' ')
    pois <- ifelse(quantity.names[X.i] %in% this$.poisson, 'P', ' ')
    flag <- sprintf('(%s%s)', sel, pois)
    cat(sprintf(format.str, quantity.names[X.i], flag,
        diff(range(raw.data[, X.i])),
        sprintf('%17.8g', this$.noise.var[X.i])))
  }
  return (invisible(this))
}

#' Remove datapoints within a given range
#'
#' Remove datapoints within a given range, for 1D data only.  This function
#' could be improved by building an engine that works on arbitrary dimensions,
#' maybe by inputting some kind of generalized metric, and using L2-norm for
#' spheres, L-oo norm for boxes, etc.  But for now: 1D only!
#'
#' @S3method RemoveRange Dataset
#' @name RemoveRange.Dataset
#'
#' @param X.min The left boundary of the X-range to remove.
#' @param X.max The right boundary of the X-range to remove.
#' @param ... Not used.
#'
#' @seealso \code{\link{Dataset}}
setMethodS3("RemoveRange", "Dataset", conflict="quiet",
  function(this, X.min=min(this$X), X.max=max(this$X), ...) {
    if (this$d > 1) {
      stop("RemoveRange() only works for 1-D data.")
    }
    if (X.min < X.max) {
      bad.indices <- which(this$X <= X.max & this$X >= X.min)
      this$.data.frame <- this$.data.frame[-bad.indices, ]
      this$.noise.var <- ComputeNoiseVariance(data=this$.data.frame,
        specs=this$.noise.var.specs, tol=this$.tol.factor,
        poisson.names=this$.poisson)
    }
    return (invisible(this))
  })

#' Check whether two Datasets are identical
#'
#' Check whether two Datasets have identical values for some subset of their
#' parameters.  Many different objects in \code{gppois} perform expensive
#' calculations on Dataset objects, then cache the results.  These objects need
#' to be able to check that nothing in a Dataset has changed... at least,
#' nothing that affects their huge calculation!  This lets these objects return
#' the cached result, saving considerable time.
#'
#' @S3method Same Dataset
#' @name Same.Dataset
#'
#' @param d The Dataset for comparison.
#' @param compare The attributes (e.g., X, noiseVar, etc.) to check for
#'      equality.
#' @param ... Not used.
#'
#' @return (logical) TRUE if all the requested attributes are the same.
#'
#' @seealso \code{\link{Dataset}}
setMethodS3("Same", "Dataset", conflict="quiet",
  function(this, d, compare, ...) {
    # First: guard against NULLs and the like!
    if (!("Dataset" %in% class(d) && "Dataset" %in% class(this))) {
      return (FALSE)
    }

    # A convenient helper function:
    field.passes.test <- function(field) {
      # Check whether a *single* 'field' is identical for Datasets 'this' and
      # 'd'.
      #
      # Args:
      #   field:  The field (e.g., X, noiseVar, etc.) to check.
      #
      # Returns:
      #   FALSE if this field is in the list that gets checked, *and* its value
      #   differs between 'this' and 'd'; otherwise, TRUE.
      if (!(field %in% compare)) {
        return (TRUE)
      }
      # Locate accessor function; check that both objects return equal values.
      substr(field, 1, 1) <- toupper(substr(field, 1, 1))
      funct <- get(paste(sep='', 'get', field))
      return (identical(funct(this), funct(d)))
    }

    # Check all the requested fields
    ok <- TRUE
    fields.to.check <- c(
      'X', 'noiseVar', 'dpts', 'xformedDpts', 'id', 'quantity')
    for (field in fields.to.check) {
      ok <- ok && field.passes.test(field)
    }

    # Warn user if any fields were requested that we don't check!
    not.checked <- compare[-which(compare %in% fields.to.check)]
    if (length(not.checked) > 0) {
      offenders <- paste(sep='', "'", not.checked, "'", collapse=', ')
      warning(sprintf(
          "The following fields are not checked by Dataset$Same:\n%s\n",
          offenders))
    }
    return (ok)
  })

#' Undo any transformations on the datapoints
#'
#' Datapoints are often transformed prior to being analyzed.  Examples include\cr
#'   * the Anscombe transform for Poisson-noised data, which makes it
#'     approximately Gaussian with constant variance\cr
#'   * subtracting off the mean\cr
#' This function provides a (statistical!) inverse for these transforms.  It is
#' intended to be applied to the results of the analysis (which take place in
#' the \emph{transformed} space).
#'
#' @S3method Untransform Dataset
#' @name Untransform.Dataset
#'
#' @param values A numeric vector whose values "live" in the
#'      \emph{transformed} scale.
#' @param ... Not used.
#'
#' @return Numbers corresponding to \code{values}, but in a non-transformed
#'    scale (i.e., instead of being comparable to \code{xformedDpts}, the
#'    returned quantities are comparable to \code{dpts}).
#'
#' @seealso \code{\link{Dataset}}
#' @seealso \code{\link{xformedDpts}}
#' @seealso \code{\link{dpts}}
setMethodS3("Untransform", "Dataset", conflict="quiet",
  function(this, values, ...) {
    values <- values + this$dataOffset
    if (this$getIsPoisson()) {
      values <- AnscombeInverse(values)
    }
    return (values)
  })

