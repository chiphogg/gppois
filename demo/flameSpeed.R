require("gppois")

# Load the data:
data(flameSpeed)
ratio.range <- range(flameSpeed$fuelRatio)
x.out <- seq(from=min(ratio.range), to=max(ratio.range), length.out=200)

# Plot the data
ggplot2.installed <- require("ggplot2")
if (ggplot2.installed) {
  p <- ggplot(data=flameSpeed, aes(x=fuelRatio, y=speed, colour=source))
  print(p + geom_point())
} else {
  with(flameSpeed, plot(fuelRatio, speed))
}

###############################################################################
# In gppois, we wrap our data in a Dataset object.  For help, type:
# > ?Dataset
DemoPause()
d.flame <- Dataset(data=flameSpeed, column="speed", X.names="fuelRatio")

###############################################################################
# Now, it's time to set up a Model for this system.
# We'll start with an empty Model object, then add some contributions.
DemoPause()
m <- Model(id="simple")
print(m)

###############################################################################
# From the plot, it looks like a simple SE covariance will work well.
# The SE uses two "feature lengthscale" params: horizontal and vertical.
# We'll give generous but reasonable boundaries, based on eyeballing the plot.
DemoPause()
c.SE <- CovarianceSE(id="signal",
  ell.bounds=c(0.2, 2), sigma.f.bounds=c(1, 50))
m$AddCovariance(c.SE)
print(m)

###############################################################################
# We also need to add a contribution which models the noise.
# One way to do this is to construct a CovarianceNoise object and add it to m.
# A simpler way is with the SetNoiseBounds function for the Model object.
# Again, we will guess a plausible range based on eyeballing the plot.
DemoPause()
m$SetNoiseBounds(sigma.vals=c(1, 10))
print(m)

###############################################################################
# Notice that the parameters have been assigned default values.
# These correspond to the *geometric* mean of the lower and upper bounds.
# There is no reason to think these values are any good!
# Let's plot the model's fit and see...
#m$PosteriorInterval(d=d.flame, X.out=x.out)
