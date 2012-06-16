############################################################################/**
# @RdocClass CovarianceSEAniso2D
#
# @title "2D Anisotropic SE Covariance"
#
# \description{
#   A 2D squared-exponential covariance whose eigenvalues of the covariance
#   matrix are not assumed identical.  In other words, features might vary more
#   rapidly in one direction than in the orthogonal direction.
#
#   @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#   \item{id}{(character) A string to identify this covariance object.}
#   \item{ell.1}{(numeric) One characteristic horizontal scale for features in
#      the functions.}
#   \item{ell.2}{(numeric) One characteristic horizontal scale for features in
#      the functions.}
#   \item{sigma.f}{(numeric) A characteristic vertical scale for features in
#      functions being modeled.}
#   \item{ell.bounds}{(numeric) The range of values which \code{ell} might
#      assume.}
#   \item{sigma.f.bounds}{(numeric) The range of values which \code{sigma.f}
#      might assume.}
#   \item{...}{Not used.}
# }
#
# \section{Covariance Parameters}{
#   This section lists the fit parameters corresponding to this type of
#   Covariance.  Any parameters marked as \dQuote{(Scale parameter)} will be
#   optimized in log-space, consistent with the Jeffreys prior.
#
#   \describe{
#     \item{ell}{(Scale parameter) The horizontal feature lengthscale.}
#     \item{sigma.f}{(Scale parameter) The vertical feature lengthscale.}
#   }
# }
#
# \section{Fields and Methods}{
#  @allmethods
#
# }
#
# @author
#*/###########################################################################
setConstructorS3("CovarianceSEAniso2D",
  function(..., id="Aniso2D", ell.1=NA, ell.2=NA, theta.1=NA, sigma.f=NA,
    ell.1.bounds=NA, ell.2.bounds=NA, theta.1.bounds=NA, sigma.f.bounds=NA) {
    pos.def.range <- c(0, Inf)
    ell.1.good <- InitializeBoundedQuantity(ok.range=pos.def.range,
      quantity=ell.1, bounds=ell.1.bounds, logspace=TRUE)
    ell.2.good <- InitializeBoundedQuantity(ok.range=pos.def.range,
      quantity=ell.2, bounds=ell.2.bounds, logspace=TRUE)
    theta.1.good <- InitializeBoundedQuantity(ok.range=pi * c(-1, 1),
      quantity=theta.1, bounds=theta.1.bounds)
    sigma.f.good <- InitializeBoundedQuantity(ok.range=pos.def.range,
      quantity=sigma.f, bounds=sigma.f.bounds, logspace=TRUE)
    # Construct the CovarianceSEAniso2D object:
    extend(Covariance(..., id=id), "CovarianceSEAniso2D",
      .ell.1          = ell.1.good$quantity,
      .ell.1.bounds   = ell.1.good$bounds,
      .ell.2          = ell.2.good$quantity,
      .ell.2.bounds   = ell.2.good$bounds,
      .theta.1        = theta.1.good$quantity,
      .theta.1.bounds = theta.1.good$bounds,
      .sigma.f        = sigma.f.good$quantity,
      .sigma.f.bounds = sigma.f.good$bounds)
  })

