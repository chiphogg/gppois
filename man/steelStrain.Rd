\docType{data}
\name{steelStrain}
\alias{steelStrain}
\alias{steelStrainGap}
\title{Strain on a stretched steel plate}
\format{Two data.frames with 3 variables: \code{steelStrain} has 3460 rows of
    training data, and \code{steelStrainGap} has 245 rows of test data
    (representing the datapoints in the gap)}
\description{
  A dataset containing digital image correlation (DIC)
  measurements of strain from a stretched steel plate.
  Measurements were taken by Adam Creuziger at the National
  Institute of Standards and Technology.  Subsequently,
  Charles Hogg subtracted off a best-fit quadratic, then
  removed two regions: \itemize{ \item a 6 mm square from
  the centre (this would be absent in a real measurement)
  \item all datapoints beyond 13 mm from the centre
  (far-away datapoints are unhelpful in predicting the
  centre, but slow down the Gaussian Process considerably).
  }
}
\details{
  \itemize{ \item X: X-coordinate on the plate \item Y:
  Y-coordinate on the plate \item exx: xx-compenent of the
  strain \item exy: xy-compenent of the strain \item eyy:
  yy-compenent of the strain }
}
\keyword{datasets}

