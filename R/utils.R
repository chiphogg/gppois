#' Runs a function with just-in-time debugging
#'
#' This function tries running \code{FUN(...)} normally.  If it's fine, the user
#' should see no difference.  If not, it catches the exception, turns on
#' debugging, and re-runs to assist debugging.  Debugging is turned off for this
#' function at the end.
#'
#' @param FUN  The function to call (either a character string or just the name
#'    of the function).
#' @param ...  The argument list for \code{FUN.}
#'
#' @return The result of calling \code{FUN(...)}, or else the exception.
DebugIfError <- function(FUN, ...) {
  debug.exists <- !isTRUE(
    suppressMessages(suppressWarnings(require("debug"))) == FALSE)
  if (is.character(FUN)) {
    FUN <- get(FUN)
  }
  arg.list <- list(...)
  value <- tryCatch(
    expr=FUN(...),
    error=function(e) {
      save(arg.list, file="exact_arguments.rda")
      # Make sure we can load the debug package
      if (debug.exists) {
        cat(sprintf("Error encountered! '%s'\nTrying again in debug mode.\n",
            e$message))
        mtrace(FUN)
        value <- FUN(...)
        mtrace(FUN, FALSE)
      } else {
        value <- e
      }
      return (value)
    }, finally={
      if (debug.exists) {
        mtrace(FUN, FALSE)
      }
    })
  return (value)
}

#' Fast trace for matrix product
#'
#' Computes the trace of the matrix product \code{m1 m2}, without actually
#' evaluating that product.
#'
#' @param m1 One matrix.
#' @param m2 The other matrix.
#'
#' @export
#' @return Trace(\code{m1 m2})
SmartTrace <- function(m1, m2) {
  return (sum(as.vector(m1) * as.vector(t(m2))))
}

ItTakes <- function(my.task, how.i.do.it) {
  # Measure how long it takes to do a given task, and print it out in
  # human-readable format.
  #
  # Args:
  #   my.task:  A short label describing what I'm trying to do
  #   how.i.do.it:  A code block; the implementation of my.task.
  #
  # Returns:
  #   Used for its side-effect.
  p.env <- parent.frame()
  cat(sprintf("%60s: ", my.task))
  elapsed <- system.time(eval(substitute(how.i.do.it), env=p.env))
  cat(sprintf("%12g sec\n", elapsed["user.self"]))
}

#' Plots a numeric matrix
#'
#' Makes a quick plot of numeric matrix data.  Not very customizable, but handy
#' for a first glance.
#'
#' @param M   The numeric matrix to plot.
#' @param ...   Unused
#'
#' @note Adapted from:
#' \url{http://www.phaget4.org/R/image_matrix.html} 2011-01-06 15:35
PlotMatrixQuickAndDirty <- function(M, ...) {

  min <- min(M)
  max <- max(M)
  yLabels <- rownames(M)
  xLabels <- colnames(M)
  title <-c()
  # check for additional function arguments
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
      min <- Lst$zlim[1]
      max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ){
      yLabels <- c(Lst$yLabels)
    }
    if( !is.null(Lst$xLabels) ){
      xLabels <- c(Lst$xLabels)
    }
    if( !is.null(Lst$title) ){
      title <- Lst$title
    }
  }
  # check for null values
  if( is.null(xLabels) ){
    xLabels <- c(1:ncol(M))
  }
  if( is.null(yLabels) ){
    yLabels <- c(1:nrow(M))
  }

  layout(matrix(data=c(1, 2), nrow=1, ncol=2), widths=c(4, 1), heights=c(1, 1))

  # Red and green range from 0 to 1 while Blue ranges from 1 to 0
  ColorRamp <- rgb( seq(0, 1, length=256),  # Red
    seq(0, 1, length=256),  # Green
    seq(1, 0, length=256))  # Blue
  ColorLevels <- seq(min, max, length=length(ColorRamp))

  # Reverse Y axis
  reverse <- nrow(M) : 1
  yLabels <- yLabels[reverse]
  M <- M[reverse, ]

  # Data Map
  par(mar = c(3, 5, 2.5, 2))
  image(1:length(xLabels), 1:length(yLabels), t(M), col=ColorRamp, xlab="",
    ylab="", axes=FALSE, zlim=c(min, max))
  if( !is.null(title) ){
    title(main=title)
  }
  axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=0.7)
  axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,
    cex.axis=0.7)

  # Color Scale
  par(mar = c(3, 2.5, 2.5, 2))
  image(1, ColorLevels,
    matrix(data=ColorLevels, ncol=length(ColorLevels), nrow=1),
    col=ColorRamp,
    xlab="", ylab="",
    xaxt="n")

  layout(1)
}

#' Computes the width of each point in X.
#'
#' For a sorted sequence of one-dimensional points X, calculates the width of
#' each point's "territory" (i.e., the set of points closer to that point than
#' any other).  The sum of all widths is the range of X.
#'
#' @param X   A sorted numeric 1-D vector, representing points on a line.
#'
#' @export
#' @return The width of the "territory" for each point in X.
Widths <- function(X) {
  N <- length(X)
  return (diff(c(X[1], (X[-N] + X[-1]) * 0.5, X[N])))
}

#' Pairwise distances between points in X and X.out
#'
#' Computes the distance from every point in X to every point in X.out.  Both
#' arguments are assumed to be numeric matrices with as many columns as the
#' dimensionality of the space.  (i.e., N 2D points would be represented by an
#' (N x 2) matrix, etc.)  Vector arguments are assumed to be 1D points, and are
#' automatically converted to matrices.
#'
#' @param X   A numeric matrix of input points.
#' @param X.out   A matrix of output points, whose distance to every point in
#'    'X' is desired.
#'
#' @export
#' @return A matrix whose [i, j] component gives the Euclidean distance from
#'   X.out[i, ] to X[j, ].
DistanceMatrix <- function(X, X.out=X) {
  X <- as.matrix(X)
  X.out <- as.matrix(X.out)
  if (ncol(X) != ncol(X.out)) {
    stop("Trying to find distance between points of different dimensions!")
  }
  X.dim <- ncol(X)
  N <- NumPoints(X)
  N.out <- NumPoints(X.out)

  # Construct two matrices with (N * N.out) rows, so that X's matrix holds
  # N.out copies of each X-point, and X.out's matrix holds N copies of each
  # X.out-point, and every point in X is matched with every point in X.out
  # exactly once.
  M <- matrix(rep(X, each=N.out), ncol=X.dim)
  M.out <- t(matrix(rep(t(X.out), N), nrow=X.dim))
  dist.sq <- rowSums((M - M.out) ^ 2)
  return (matrix(sqrt(dist.sq), nrow=N.out, ncol=N))
}


#' Clamps \code{x} to within \code{bounds}
#'
#' Clamps every value in \code{x} (numeric) to lie within the range of
#' \code{bounds} (numeric).
#'
#' @param x Some numeric object whose values must be clamped.
#' @param bounds A numeric object whose range sets the range for clamping.
#'
#' @return \code{x}, but with values clamped according to \code{bounds}.
Clamp <- function(x, bounds) {
  bounds <- range(bounds)
  return (pmax(pmin(x, bounds[2]), bounds[1]))
}

ClampNamed <- function(x, lower, upper) {
  # Clamps every value in 'x' to lie between the values of 'lower' and of
  # 'upper' with the same name.
  #
  # Args:
  #   x: A named numeric vector to be clamped.
  #   lower: A numeric vector holding the lower bounds.
  #   upper: A numeric vector holding the upper bounds.
  #
  # Returns:
  #   'x', but with values clamped between 'lower' and 'upper'.
  n.x <- names(x)
  return (pmax(pmin(x[n.x], upper[n.x]), lower[n.x]))
}

InitializeBoundedQuantity <- function(quantity, bounds,
  default=1, ok.range=c(-Inf, Inf), logspace=FALSE) {
  # Process values for a quantity and its boundaries so they are mutually
  # consistent.
  #
  # Args:
  #   quantity:  The bounded quantity, a single numeric variable.
  #   bounds:  Numeric vector of arbitrary length, whose range() gives the
  #      lower and upper bounds.
  #   default:  The default value for quantity, in case bounds and quantity are
  #      both missing or inappropriate.
  #   ok.range:  The range of values which make sense for this quantity (e.g.,
  #      if it is positive-definite, ok.range should be c(0, Inf)).
  #   logspace:  If TRUE, we set the initial value to the *geometric* mean.
  #
  # Returns:
  #   A list x, where x$quantity is the processed value of quantity, and
  #   x$bounds is the processed value of bounds, such that:
  #   (x$bounds[1] <= x$quantity <= x$bounds[2]).
  #
  # Notes:
  #   1) If both quantity and bounds are given, clamp the quantity to the
  #      bounds.
  #   2) If only the quantity is given, set the bounds to the quantity itself
  #      (effectively making it a constant).
  #   3) If only the bounds are given, set the quantity to their mean.
  #   4) If neither is given, set bounds to range, quantity to default.

  have.quantity <- is.numeric(quantity) && length(quantity) == 1
  have.bounds <- is.numeric(bounds)
  if (have.bounds) {
    bounds <- Clamp(x=range(bounds), bounds=ok.range)
    if (!have.quantity) {
      if (logspace) {
        quantity <- exp(mean(log(bounds)))
      } else {
        quantity <- mean(bounds)
      }
    }
    quantity <- Clamp(x=quantity, bounds=bounds)
  } else {
    if (have.quantity) {
      quantity <- Clamp(x=quantity, bounds=ok.range)
      bounds <- range(quantity)
    } else {
      bounds <- ok.range
      quantity <- default
    }
  }
  return (list(quantity=quantity, bounds=bounds))
}

NumPoints <- function(X) {
  # Count the number of points represented by X, hopefully with a little
  # robustness as to its type (matrix, vector, etc.).
  #
  # Args:
  #   X:  A set of points: ideally, a matrix with the number of columns equal
  #      to the dimensionality, though for 1D points it could also turn out to
  #      be a vector.
  #
  # Returns:
  #   The number of points in X.
  N <- length(X)
  if (is.matrix(X)) {
    N <- nrow(X)
  }
  return (N)
}

SetupFileInfo <- function(name, type, name.default="DEFAULT") {
  # Given a filename and file type, either of which may be missing or
  # malformed, return a proper filename and filetype for a file to write to.
  #
  # Args:
  #   name:  Character; the name of the file to save.
  #   type:  Character; the file type (pdf, png, etc.)
  #
  # Returns:
  #   list;
  #     $name = character; the filename to save to.
  #     $type = character; the type of file to save.
  have.name <- is.character(name) && length(name) == 1
  have.type <- is.character(type) && length(type) == 1
  if (!have.name) {
    name <- name.default
  }
  ext.match <- gregexpr(name, pattern="\\.[^.]+", perl=TRUE)
  z <- length(ext.match[[1]])
  name.has.ext <- (ext.match[[1]][1] > 0)
  if (!have.type) {
    if (name.has.ext) {
      ext.start <- ext.match[[1]][z]
      ext.stop <- ext.start + attr(ext.match[[1]], "match.length")[z]
      type <- substr(x, ext.start + 1, ext.stop)
    } else {
      type <- ""
    }
  }
  if (!name.has.ext && (type != "")) {
    name <- paste(sep='', name, '.', type)
  }
  return (list(name=name, type=type))
}

#' Generate marginally normal smooth random timetraces
#'
#' This function is an "engine" for visualizing uncertainty using animations --
#' at least, as long as that uncertainty is normally distributed.  The output is
#' a matrix whose rows represent independent oscillators, and whose columns
#' represent timesteps.  Every column is i.i.d. normal.  Every row has a CDF
#' which approaches the normal distribution for large N, but which is
#' nevertheless smooth to all orders of time derivatives.  Increasing N improves
#' the statistical properties, while increasing n.times relative to N yields
#' smoother animations.
#'
#' @param n.pts  The number of independent oscillators.
#' @param N  One-half the number of independent draws to take for each
#'    oscillator.
#' @param n.times  The final number of interpolated points.
#'
#' @export
#' @return  A numeric matrix with n rows and n.times columns, showing the
#'    timetraces of the oscillators.
BubblingRandomMatrix <- function(n.pts, N=10, n.times=100, ...) {
  # First: "seed matrix", which contains enough independent random draws
  max.time <- n.draws <- 2 * N
  d.t <- max.time / n.times
  times <- seq(from=d.t, to=max.time, by=d.t)
  n.seed <- n.pts * n.draws
  seed <- matrix(rnorm(n=n.seed), ncol=n.draws)

  # Next: "Fourier matrix", which gives trigonometric basis functions
  F.sin <- t(sapply(X=1:N, FUN=function(i, N, times) sin(pi * i * times / N),
    N=N, times=times))
  F.cos <- t(sapply(X=1:N, FUN=function(i, N, times) cos(pi * i * times / N),
    N=N, times=times))
  F.tot <- rbind(F.sin, F.cos) / sqrt(N)

  # Put them together and we're done!
  return (seed %*% F.tot)
}

Spaces <- function(num=0) {
  # Constructs a string consisting of the requested number of spaces.  (Useful
  # for recursive tab alignment of printed output.)
  #
  # Args:
  #   num:  The number of spaces to return.
  #
  # Returns:
  #   A string with 'num' spaces.
  return (paste(collapse='', rep(' ', num)))
}

Wrap <- function(x, wrapper) {
  # Wrap each entry in 'x' with the string in 'wrapper'
  #
  # Args:
  #   x:  Character vector whose entries are to be wrapped.
  #   wrapper:  Character string to wrap it in.
  #
  # Returns:
  #   Strings in the format "<wrapper><x><wrapper>".
  return (paste(sep='', wrapper, x, wrapper))
}

PrintParams <- function(lower, params, upper, indent=0, ...) {
  # Print out a list of parameters, together with their upper and lower bounds.
  #
  # Args:
  #   lower: A numeric vector holding the lower bounds.
  #   params: A numeric vector holding the parameter values.
  #   upper: A numeric vector holding the upper bounds.
  # Ruler:
  #     1...5...10...15...20...25...30...35...40...45...50...55...60...65...70
  s <- "  PARAMETER NAME:         MINIMUM:   CURRENT VALUE:         MAXIMUM:\n"
  tab <- Spaces(num=indent)
  cat(tab, s, sep='')
  p.names <- names(params)
  for (name in p.names) {
    cat(sprintf("%s  %15s   %14.8g   %14.8g   %14.8g\n", tab, name,
        lower[name], params[name], upper[name]))
  }

}

#' Error function
#'
#' The standard function 'erf', based on the CDF of the standard normal
#' distribution.
#'
#' @param x The x-value for erf(x)
#'
#' @export
#' @return  erf(x)
erf <- function(x) {
  2 * pnorm(sqrt(2) * x) - 1
}

#' Pause during demos
#'
#' Prompts the user for input, to break up the demos
#'
#' @param msg The prompt message to display
#' @export
DemoPause <- function(msg="Press Enter to continue:") {
  cat(msg)
  readline()
  cat("\n")
}
