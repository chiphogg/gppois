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
#' @format A data.frame with 3 variables and 3705 rows
NULL
