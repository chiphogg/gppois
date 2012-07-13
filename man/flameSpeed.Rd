\docType{data}
\name{flameSpeed}
\alias{flameSpeed}
\title{Flame speed vs. fuel concentration}
\format{A data.frame with 3 variables and 78 rows}
\description{
  Flame speed measurement is a standard experiment, useful
  for calibrating computer models of combustion.  A mixture
  of fuel and oxygen is fed continuously into a long,
  straight tube.  The tube is ignited at one end, and the
  speed of the resulting flame front is measured.  The
  flame is fastest for some optimal fuel-to-air ratio, so
  the data look roughly like a "hump".
}
\details{
  These results were aggregated from a variety of papers in
  the literature. Since the total number of datapoints is
  small -- less than 100 -- it constitutes a very fast
  example system to illustrate Gaussian Process regression.

  \itemize{ \item fuelRatio: Ratio of fuel to oxygen (or
  oxygen to fuel?)  \item source: The source in the
  literature of each measurement \item speed: The speed of
  the flame front }
}
\keyword{datasets}

