context("Utility functions")

test_that("SmartTrace really does give the trace", {
    m <- 1000
    n <- 2000
    m.n <- m * n
    A <- matrix(runif(n=m.n, min=0, max=100), nrow=m, ncol=n)
    B <- matrix(runif(n=m.n, min=0, max=100), nrow=n, ncol=m)
    dumb.time <- system.time(trace.AB <- sum(diag(A %*% B)))
    smrt.time <- system.time(smart.tr <- SmartTrace(A, B))
    expect_equal(trace.AB, smart.tr)
    expect_equivalent(smrt.time["user.self"] < dumb.time["user.self"], TRUE)
})

test_that("ClampNamed clamps and preserves names", {
    x <- c(a=1, b=2, c=3)
    lower <- c(c=1, b=3, a=5)
    upper <- c(8, a=8, b=3, 9, c=2)
    Clamped <- ClampNamed(x=x, lower=lower, upper=upper)
    expect_that(names(x), is_equivalent_to(names(Clamped)))
    n <- names(x)
    expect_that(Clamped[n], is_identical_to(c(a=5, b=3, c=2)[n]))
})

test_that("Widths() computes the correct interval sizes", {
    N <- 10
    x <- 1:N
    dx <- Widths(x)
    expect_equal(sum(dx), diff(range(x)))
    expect_equal(dx[1] * 2, dx[2])
    expect_equal(dx[N] * 2, dx[N - 1])
})
