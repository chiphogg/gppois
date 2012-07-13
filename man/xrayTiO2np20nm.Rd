\docType{data}
\name{xrayTiO2np20nm}
\alias{xrayTiO2np20nm}
\title{Powder x-ray diffraction from 20nm TiO2 NPs}
\format{A data.frame with 3 variables and 10769 rows}
\description{
  Powder X-ray diffraction data from titania (TiO2)
  nanoparticles, 20 nm in diameter.  Includes a high-noise
  version (low count times) and a lower-noise version as a
  proxy for the true function (which is of course unknonwn
  in experimental settings).  The latter function is based
  on a long-count-time run, corrected for various
  systematic effects by procedures described in the
  Supplemental Material for Hogg et al. (2012).
}
\details{
  \itemize{ \item Q: Q-values where we have measurements.
  (Q is known as the \dQuote{wavenumber}; it is the
  independent variable in scattering experiments.)  \item
  noisy: Noisy measurements of the scattered x-ray
  intensity, in units of counts.  \item low.noise: A
  long-count-time version of the measurement, processed to
  correct for source intensity fluctuations and the like. }
}
\references{
  Hogg, C., K. Mullen, and I. Levin (2012). A Bayesian
  approach for denoising one-dimensional data. Journal of
  Applied Crystallography, 45(3), pp. 471-481.
}
\keyword{datasets}

