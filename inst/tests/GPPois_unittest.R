# This function contains some unit tests I had written before I knew how to
# make R packages.  I am currently (2012-07-12) in the process of integrating
# them into the package, but much has changed and the task requires both time
# and care.

# This file should not be called.  It won't be, because it doesn't start with
# "test" or "helper", but just in case:
return()

################################################################################
# LazyMatrix tests

context("LazyMatrix")

test_that("I have written tests for LazyMatrix", {
    expect_that(FALSE, is_true())
})


################################################################################
# Dataset tests

# Strain example
adams.data <- DatasetFromFile(file="Creuziger_strain.dat", sep=',',
  X.names=c("X", "Y"), column="exx", tol.factor=1e-4)

test_that("Dataset objects can be created from datafiles", {
    expect_that(adams.data$id, is_equivalent_to("Creuziger_strain.dat"))
    expect_that(adams.data$d, is_equivalent_to(2))
    expect_that(adams.data$quantity, is_equivalent_to("exx"))
})

# Flame speed example
daves.data.frame <- read.table(sep='\t', header=TRUE, file="flame_speed.dat")
daves.melted <- melt(data=daves.data.frame, na.rm=TRUE, id.vars=1)
daves.data <- Dataset(data=daves.melted, column="value")
# Handy test object: same exact data, but shuffled, so that 
daves.melted.shuffled <- daves.melted[sample(daves.data$n),]
rownames(daves.melted.shuffled) <- rownames(daves.melted)
daves.data.shuffled <- Dataset(data=daves.melted.shuffled, column="value")
# Objects for plotting the flame data
plot.flame.base <- (
  ggplot(data=daves.melted, aes(x=X, y=value, colour=variable, fill=variable))
  + scale_colour_brewer()
  + scale_fill_brewer()
  )
plot.flame.smoothed <- (plot.flame.base
  + geom_point()
  + stat_smooth()
  )

test_that("Noise variance is handled correctly", {
    custom.tol <- 1e-2
    igors.data <- DatasetFromFile(file="TiO2_powder.dat", id="TiO2Powder", sep='\t',
      X.names="Q", poisson.names=c("noisy", "wavelet"), tol.factor=custom.tol,
      noise.var=c(wavelet=0.1, true=0.1234))
    # First, assuming we auto-selected the first non-X column, see that the
    # variance is 1/4 as expected
    expect_that(igors.data$quantity, is_equivalent_to("noisy"))
    expect_that(igors.data$isPoisson, is_equivalent_to(TRUE))
    expect_that(igors.data$noiseVar, is_equivalent_to(0.25))
    expect_that(sqrt(igors.data$dpts + 3/8),
      is_equivalent_to(igors.data$xformedDpts))
    # Now, check the auto-set noise variance :
    igors.data$quantity <- "waveletAnscombe"
    expect_that(igors.data$quantity, is_equivalent_to("waveletAnscombe"))
    noise.unit <- sd(igors.data$xformedDpts)
    expect_that(igors.data$noiseVar, is_equivalent_to((noise.unit * custom.tol)^2))
    # And check a manually-set noise variance:
    igors.data$quantity <- "true"
    expect_that(igors.data$noiseVar, is_equivalent_to(0.1234))
})

################################################################################
# Covariance tests

test_that("We can change upper and lower bounds for CovarianceSE ", {
    # First: does constructor clamping work?
    Cov <- CovarianceSE(id="bounds.test", ell=5, sigma.f=-4,
      ell.bounds=1:4, sigma.f.bounds=0:100)
    expect_that(Cov$paramsPlain["ell"], is_equivalent_to(4))
    expect_that(Cov$paramsPlain["sigma.f"], is_equivalent_to(0))
    # Pushing lower bound past upper bound should move everything
    expect_that(Cov$lowerPlain <- c(sigma.f=200), gives_warning())
    expect_that(Cov$lowerPlain["sigma.f"], is_equivalent_to(200))
    expect_that(Cov$upperPlain["sigma.f"], is_equivalent_to(200))
    expect_that(Cov$paramsPlain["sigma.f"], is_equivalent_to(200))
    # Test setter-methods for both plain and decorated param names
    Cov$upper <- c(bounds.test.sigma.f=1000)
    expect_that(Cov$upper["bounds.test.sigma.f"], is_equivalent_to(1000))
    expect_that(Cov$upperPlain["sigma.f"], is_equivalent_to(1000))
    expect_that(Cov$params <- c(bounds.test.sigma.f=2000), gives_warning())
    expect_that(Cov$paramsPlain["sigma.f"], is_equivalent_to(1000))
    # These are positive-definite quantities: make sure this is enforced
    expect_that(Cov$upperPlain <- c(ell=-1), gives_warning())
    expect_that(Cov$lowerPlain["ell"], is_equivalent_to(0))
})

test_that("Validity of K-matrix is properly tracked", {
    # Create a covariance matrix based on previously-calculated optimal values:
    # http://www-i.nist.gov/msel/ceramics/data/projects/logsheets/AdamStressStrain/Log.html
    Cov <- CovarianceSE(id="K.test", ell=1.26, sigma.f=sqrt(1.482e-6))
    adams.data$quantity <- "exx"
    adams.other.data <- clone(adams.data)
    # First, compute the K-matrix based on Adam's data.
    base.time <- system.time(K.base <- Cov$KInIn(adams.data))["user.self"]
    # I don't know exactly how much faster it should be... but I'd say for sure
    # if it takes less than half the time, that's a good sign we didn't
    # recompute the matrix.
    time.thresh <- 0.5 * base.time
    adams.data$quantity <- "eyy"
    other.object.time <- system.time(K <- Cov$KInIn(adams.other.data))["user.self"]
    expect_equivalent(K, K.base)
    expect_equivalent(other.object.time < time.thresh, TRUE)
    # Now let's change the covariates and make sure it recomputes everything.
    # Note that this is the SAME object, so if Covariance uses naive equality
    # checking (i.e. checking based on references), this test will fail.  (It
    # *should* keep a *clone* of the most recent dataset.)
    adams.other.data$DeleteRows(1)
    changed.data.time <- system.time(Cov$KInIn(adams.other.data))["user.self"]
    expect_equal(adams.other.data$n, adams.data$n - 1)
    expect_equivalent(changed.data.time > time.thresh, TRUE)
    # And when we recalculate it, of course it should be "fast".
    recalc.data.time <- system.time(K <- Cov$KInIn(adams.other.data))["user.self"]
    expect_equivalent(recalc.data.time < time.thresh, TRUE)
    # Changing one of the params should make it slow again.
    old.ell <- Cov$paramsPlain["ell"]

    # We shouldn't be able to change a param out of bounds:
    expect_that(Cov$paramsPlain <- c(ell=unname(old.ell) * 2), gives_warning())
    # The right way to do it:
    Cov$FixConstParam(p.name="ell", p.value=old.ell * 2)
    expect_that(Cov$paramsPlain["ell"], is_equivalent_to(old.ell * 2))

    newparam.data.time <- system.time(
      K.new.p <- Cov$KInIn(adams.other.data)
      )["user.self"]
    expect_equivalent(newparam.data.time > time.thresh, TRUE)
    expect_equivalent(all(K.new.p >= K), TRUE)

    ## Debugging output
    #fstr <- "%40s: %f sec\n"
    #cat(sprintf(paste(c("\n", rep(fstr, 5)), collapse=''),
    #    "Base", base.time,
    #    "Cloned Dataset", other.object.time,
    #    "Cloned Dataset, changed", changed.data.time,
    #    "Recalc with changed cloned object", recalc.data.time,
    #    "Changed parameters", newparam.data.time))

})

################################################################################
# CovarianceSE tests

context("CovarianceSE")

test_that("Derivatives of K are correctly calculated", {
    # Generate a Covariance object and calculate K and derivatives.
    Cov <- CovarianceSE(ell=1, sigma.f=1,
      ell.bounds=c(0.1, 10), sigma.f.bounds=c(0.1,10))
    model <- Model(id="deriv.test")
    model$AddCovariance(Cov)
    change <- 1e-5
    for (p.n in names(model$params)) {
      model.new <- clone(model)
      model.new$params <- model.new$params[p.n] + change
      K.deriv <- model$KDeriv(d=daves.data, param=p.n)
      expect_that(max(abs(range(model.new$KTotal(d=daves.data) - (
                model$KTotal(d=daves.data) + (K.deriv * change)))))
        < 1e-8, is_true())
    }
})

################################################################################
# CovarianceSEAniso2D tests

context("CovarianceSEAniso2D")

test_that("Derivatives of K are correctly calculated", {
    expect_that(FALSE, is_true())
}

################################################################################
# Model tests

context("Model")

test_that("I can create a Model object", {
    model.1.cov <- Model(id="Single Covariance")
    expect_that(model.1.cov$id, is_equivalent_to("Single Covariance"))
    expect_that(length(model.1.cov$params), is_equivalent_to(0))
    # Now, add a simple covariance to it
    model.1.cov$AddCovariance(
      CovarianceSE(ell.bounds=c(0.1, 10), sigma.f.bounds=1e-3 * c(0, 10)))
    expect_that(length(model.1.cov$params), is_equivalent_to(2))
    expect_that(model.1.cov$contributionIds, is_equivalent_to("SE"))
    # Test collision detection in contribution id's
    model.1.cov$AddCovariance(
      CovarianceSE(ell.bounds=c(0.1, 10), sigma.f.bounds=1e-3 * c(0, 10)))
    expect_that(model.1.cov$contributionIds, is_equivalent_to(c("SE", "SE.1")))
})

test_that("GradLogML really is gradient of LogML", {
    model <- Model(id='GradLogML.test')
    model$SetNoiseBounds(c(1, 10))
    model$AddCovariance(
      CovarianceSE(ell.bounds=c(0.2, 2), sigma.f.bounds=c(1, 100)))
    logML.old <- LogML(model=model, d=daves.data)
    grad <- GradLogML(model=model, d=daves.data)
    n <- length(model$params)
    rel.change <- 1e-6
    change <- rnorm(n=n, sd=rel.change)
    model.new <- clone(model)
    model.new$params <- model$params + change
    logML.new <- LogML(model=model.new, d=daves.data)
    actual <- as.numeric(logML.new - logML.old)
    expected <- as.numeric(sum(change * grad))
    expect_equal(actual, expected, tolerance=1e-4)
})

test_that("I can train a Model on a Dataset", {
    # These parameter bounds come from eyeballing the data; try the command:
    # > print(plot.flame.smoothed)
    flame.model <<- Model(id="Flame Speed")
    flame.model$SetNoiseBounds(c(1, 10))
    flame.model$AddCovariance(
      CovarianceSE(ell.bounds=c(0.2, 2), sigma.f.bounds=c(1, 100)))

    # Check: training the model should be faster the second time around
    # (because it has already been done):
    N <- 5
    clone.offset <- system.time({
        for (i in 1:N) {
          clone(flame.model)
        }
      })["user.self"]
    # How much time to calculate "cold"? (i.e., not in vicinity of optimum)
    first.time <- (system.time({
          for (i in 1:N) {
            clone(flame.model)$Train(d=daves.data)
          }
          flame.model$Train(d=daves.data)
        })["user.self"] - clone.offset) / (N + 1)
    # How long to calculate when we start close to the optimum?  (Multiply by
    # 0.5 to avoid randomness... we want to be REALLY sure we have a good
    # benchmark to beat.)
    bench.time <- 0.5 * (system.time({
          for (i in 1:N) {
            clone(flame.model)$Train(d=daves.data.shuffled)
          }
          clone(flame.model)$Train(d=daves.data.shuffled)
        })["user.self"] - clone.offset) / (N + 1)
    # How long when we redo?  (Should be SHORT: it should see it's already
    # optimized for this Dataset, and just skip it.)
    redo.time  <- (system.time({
          for (i in 1:N) {
            clone(flame.model)$Train(d=daves.data)
          }
          flame.model$Train(d=daves.data)
        })["user.self"] - clone.offset) / (N + 1)
    expect_true(unname(redo.time < bench.time))
    # How long when we FORCE the redo?  (Should be comparable to 2*bench.time,
    # or, what bench.time *would* have been w/o multiplying by 0.5.)
    force.redo.time  <- system.time(
      flame.model$Train(d=daves.data, force.retrain=TRUE))["user.self"]
    expect_true(unname(force.redo.time > bench.time / N))
    cat(sprintf("TIMES: first=%g, bench=%g, redo=%g, force.redo=%g\n",
        first.time, bench.time, redo.time, force.redo.time))
  })
