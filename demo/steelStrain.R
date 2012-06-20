require("gppois")

data(steelStrain)

# First, setup a Dataset object from the steel strain data.
# If rgl is installed, plot the raw datapoints.
#
# NOTE: the rgl plot window is interactive!  Try rotating with the left mouse
# button.  Right and middle buttons change zoom and perspective.
# Get a feel for the data...
DemoPause()

d.strain <- Dataset(id="steel.strain", data=steelStrain, X.names=c("X", "Y"),
  column="exx", data.offset=0)
d.gap <- Dataset(id="gap.points", data=steelStrainGap, X.names=c("X", "Y"),
  column="exx", data.offset=0)
rgl.installed <- require("rgl")
if (rgl.installed) {
  d.strain$Plot2D(max.points=d.strain$n, Y.scale=500)
}
#
# Next, we'll setup a Model for this data.
DemoPause()

# Make the Model object
M.aniso <- Model(id="aniso")

# Setup and add the Covariance.  First, make some educated guesses about
# lengthscales.
ell.bounds <- c(0.1, 10)
sigma.f.relative <- c(0.1, 10)
sigma.n.bounds <- diff(range(d.strain$dpts)) * c(1e-7, 1e-3)
# Then, setup the Covariance object:
Cov.2d <- CovarianceSEAniso2D(id="signal", theta.1=0,
  ell.1.bounds=ell.bounds, ell.2.bounds=ell.bounds,
  sigma.f.bounds=sigma.f.relative * sd(d.strain$dpts))
# Then, add it to the Model
M.aniso$AddCovariance(Cov.2d)
M.aniso$SetNoiseBounds(sigma.n.bounds)
print(M.aniso)

# Now it's time to *Train* the Model.
# (In other words: find the Model parameters which best describe the data.)
# !! This might take a few minutes after you begin !!
DemoPause("Hit Enter to begin the training:")
ItTakes(my.task="Train the anisotropic model",
  how.i.do.it=M.aniso$Train(d=d.strain))
