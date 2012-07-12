context("Dataset")

data(steelStrain)
d.strain.base <- Dataset(id="steel.strain", data=steelStrain,
  X.names=c("X", "Y"), column="exx", data.offset=0)
data(flameSpeed)
d.flame.base <- Dataset(data=flameSpeed, column="speed")

test_that("Dataset objects are created correctly from data.frames", {
    d.strain <- clone(d.strain.base)
    expect_equivalent(d.strain$d, 2)
    expect_equivalent(d.strain$n, nrow(steelStrain))
    expect_equivalent(d.strain$X, as.matrix(ncol=2, steelStrain[, c("X", "Y")]))
    expect_equivalent(d.strain$id, "steel.strain")
    expect_equivalent(d.strain$quantity, "exx")
})

test_that("We can change the ID and selected quantity on the fly", {
    d.strain <- clone(d.strain.base)
    expect_equivalent(d.strain$id, "steel.strain")
    d.strain$id <- "new.id"
    expect_equivalent(d.strain$id, "new.id")
    expect_equivalent(d.strain$quantity, "exx")
    d.strain$quantity <- "eyy"
    expect_equivalent(d.strain$quantity, "eyy")
    expect_warning(d.strain$quantity <- "FOO", "does not exist")
    expect_equivalent(d.strain$quantity, "eyy")
})

test_that("Rows can be deleted", {
    d.strain <- clone(d.strain.base)
    n.rows <- d.strain$n
    d.strain$DeleteRows(1:10)
    expect_false(d.strain$n == n.rows)
})

test_that("We can check whether datasets have the same 'X'", {
    d.strain <- clone(d.strain.base)
    d.strain.clone <- clone(d.strain)
    d.strain.clone$quantity <- "exy"
    expect_true(d.strain$SameX(d.strain.clone))
    d.strain.clone$DeleteRows(1:10)
    expect_false(d.strain$SameX(d.strain.clone))
})

