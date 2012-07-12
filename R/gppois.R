#!/usr/bin/R

#' gppois (package)
#'
#' Gaussian Processes for Poisson-noised Data
#'
#' @name gppois
#' @docType package
#' @import R.oo
NULL

#' Strain on a stretched steel plate
#'
#' A dataset containing digital image correlation (DIC) measurements of strain
#' from a stretched steel plate.  Measurements were taken by Adam Creuziger at
#' the National Institute of Standards and Technology.  Subsequently, Charles
#' Hogg subtracted off a best-fit quadratic, then removed two regions:
#' \itemize{
#'   \item a 6 mm square from the centre (this would be absent in a real
#'         measurement)
#'   \item all datapoints beyond 13 mm from the centre (far-away datapoints are
#'         unhelpful in predicting the centre, but slow down the Gaussian
#'         Process considerably).
#' }
#'
#' \itemize{
#'   \item X: X-coordinate on the plate
#'   \item Y: Y-coordinate on the plate
#'   \item exx:  xx-compenent of the strain
#' }
#'
#' @docType data
#' @keywords datasets
#' @name steelStrain
#' @usage data(steelStrain)
#' @format Two data.frames with 3 variables: \code{steelStrain} has 3460 rows of
#'     training data, and \code{steelStrainGap} has 245 rows of test data
#'     (representing the datapoints in the gap)
NULL

#' Flame speed vs. fuel concentration
#'
#' Flame speed measurement is a standard experiment, useful for calibrating
#' computer models of combustion.  A mixture of fuel and oxygen is fed
#' continuously into a long, straight tube.  The tube is ignited at one end, and
#' the speed of the resulting flame front is measured.  The flame is fastest for
#' some optimal fuel-to-air ratio, so the data look roughly like a "hump".
#'
#' These results were aggregated from a variety of papers in the literature.
#' Since the total number of datapoints is small -- less than 100 -- it
#' constitutes a very fast example system to illustrate Gaussian Process
#' regression.
#'
#' \itemize{
#'   \item fuelRatio: Ratio of fuel to oxygen (or oxygen to fuel?)
#'   \item source: The source in the literature of each measurement
#'   \item speed: The speed of the flame front
#' }
#'
#' @docType data
#' @keywords datasets
#' @name flameSpeed
#' @usage data(flameSpeed)
#' @format A data.frame with 3 variables and 78 rows
NULL
