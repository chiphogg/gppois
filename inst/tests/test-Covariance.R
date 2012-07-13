context("CovarianceSE")

test_that("CovarianceSE is constructed correctly", {
    # First: create an object and test basic properties
    c.SE <- CovarianceSE(id="SE", ell=5, sigma.f.bounds=1:100)
    expect_equivalent(c.SE$id, "SE")
    expect_equivalent(c.SE$paramNames, c("SE.ell", "SE.sigma.f"))
    expect_equivalent(c.SE$paramNamesPlain, c("ell", "sigma.f"))
    expect_equivalent(c.SE$params, c(5, 10))
    expect_equivalent(c.SE$lowerPlain["ell"], c.SE$upperPlain["ell"])
    expect_equivalent(c.SE$lowerPlain["sigma.f"], 1)
    expect_equivalent(c.SE$upperPlain["sigma.f"], 100)
    # Now, change that object as if it were a reference
    change.lower <- function(c.SE, new.vals) {
      c.SE$lowerPlain <- new.vals
    }
    change.lower(c.SE=c.SE, new.vals=c(ell=3))
    expect_equivalent(c.SE$lowerPlain["ell"], 3)
})

