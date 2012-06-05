#############################################################################/**
# @RdocClass Dataset
#
# @title "A wrapper class for data being analyzed"
#
# \description{
#   @classhierarchy
#
#   @get "title".
# }
#
# @synopsis
#
# \section{Fields and Methods}{
#   @allmethods
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
