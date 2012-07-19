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
  literature of each measurement (see references below)
  \item speed: The speed of the flame front }
}
\references{
  \emph{VE&L 94}: C. M. Vagelopoulos, F. N. Egolfopoulos,
  C. K. Law, Proc. Combust. Inst. 25 (1994), 1341-1347.

  \emph{Gu 00}: X. J. Gu, M. Z. Haq, M. Lawes, R. Woolley,
  Combust. Flame 121 (2000), 41-58.

  \emph{B&dG 08}: K. J. Bosschaart, L. P. H. de Goey,
  Combust. Flame 136 (2004), 261-269.

  \emph{EZ&L 90}: F. N. Egolfopoulos, D. L. Zhu, C. K. Law,
  Proc.  Combust. Inst.  23 (1991), 471-478.

  \emph{HA&F 98}: M. I. Hassan, K. T. Aung, G. M. Faeth,
  Combust. Flame 115 (1998), 539-550.

  \emph{V&E 98}: C. M. Vagelopoulos, F. N. Egolfopoulos,
  Proc. Combust.  Inst. 27 (1998), 513-519.
}
\keyword{datasets}

