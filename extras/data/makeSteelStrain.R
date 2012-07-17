# This script takes the raw DIC data, and creates two data.frames, one
# corresponding to the training datapoints, and the other corresponding to the
# test datapoints in the gap.
#
# I'm not including the original datafile with the package, because it's huge
# (5.1 MB); including huge files would violate the guidelines for R packages.

norm.L2 <- function(x) apply(X=x, MARGIN=1, function(x) sqrt(sum(x ^ 2)))
norm.Loo <- function(x) apply(X=x, MARGIN=1, function(x) max(abs(x)))

fit.column.quadratic <- function(d.frame, i.ignore, col.name) {
  # Subtract off a quadratic fit for the named column,
  # ignoring the specified points.
  nogap.data <- d.frame[-i.ignore, ]
  col.num <- which(colnames(nogap.data) == col.name)
  colnames(nogap.data)[col.num] <- "VAR"
  q.fit <- with(nogap.data, lm(VAR ~ X + Y + I(X ^ 2) + I(Y ^ 2) + I(X * Y)))
  co <- q.fit$coefficients
  d.frame[, col.name] <- d.frame[, col.name] - with(d.frame,
    (co[1] + (co[2] * X) + (co[3] * Y)
      + (co[4] * X^2) + (co[5] * Y^2) + (co[6] * X * Y)))
  return (d.frame)
}

create_rda <- function(csv.filename="big_steel_dataset.csv") {
  raw.data <- read.table(file=csv.filename, header=TRUE, sep=",")
  # Remove outliers:
  i.bad <- which(raw.data$exx < mean(range(raw.data$exx)))
  raw.data <- raw.data[-i.bad, ]
  # Full set of points considered (we ignore points beyond the ram's edge,
  # since these are completely unhelpful in predicting the gap data)
  R.max <- 26
  i.far <- which(norm.L2(raw.data[, c("X", "Y")]) > R.max)
  raw.data <- raw.data[-i.far, ]
  # Don't use any gap data to train the quadratic (really shouldn't affect the
  # fit, but I'm being careful):
  max.gap <- 6  # mm
  i.gap <- which(norm.Loo(raw.data[, c("X", "Y")]) < (max.gap * 0.5))
  # Subtract off a quadratic fit for each data column
  for (col.name in c("exx", "eyy", "exy")) {
    raw.data <- fit.column.quadratic(d.frame=raw.data, i.ignore=i.gap,
      col.name=col.name)
  }
  # Now that we've got our quadratic fit, use a lot fewer points for the
  # Gaussian process.  The extra edge points don't help our predictions anyway,
  # but they reeeally slow it down!
  R.GP <- 13  # mm
  i.pred <- which(norm.L2(raw.data[, c("X", "Y")]) > R.GP)
  raw.data <- raw.data[-i.pred, ]
  # (We have to recalc the gap point indices after deleting points:)
  i.gap <- which(norm.Loo(raw.data[, c("X", "Y")]) < (max.gap * 0.5))
  ok.names <- c("X", "Y", "exx", "eyy", "exy")
  steelStrain    <- raw.data[-i.gap, ok.names]
  steelStrainGap <- raw.data[ i.gap, ok.names]
  rownames(steelStrain) <- NULL
  rownames(steelStrainGap) <- NULL
  save(steelStrain, steelStrainGap, file="../../data/steelStrain.rda")
}
