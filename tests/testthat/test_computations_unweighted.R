context("unweighted computations")

# no ties
x <- matrix(sample.int(200), 20, 10)
colnames(x) <- letters[1:10]

# with ties
xt <- x
xt[1:5, 1:2] <- xt[6:10, 1:2]
colnames(xt) <- letters[1:10]

test_that("method pearson is correct", {
    expect_equal(wdm(x), cor(x))
    expect_equal(wdm(xt), cor(xt))
})

test_that("method kendall is correct", {
    expect_equal(wdm(x, method = "kendall"), cor(x, method = "kendall"))
    expect_equal(wdm(xt, method = "kendall"), cor(xt, method = "kendall"))
})

test_that("method spearman is correct", {
    expect_equal(wdm(x, method = "spearman"), cor(x, method = "spearman"))
    expect_equal(wdm(xt, method = "spearman"), cor(xt, method = "spearman"))
})

test_that("method hoeffding is correct", {
    expect_equal(wdm(x, method = "hoeffd"), Hmisc::hoeffd(x)$D)
    # no tie correction in wdm
})

test_that("method blomqvist is correct", {
    expect_equal(
        wdm(x[, 1], x[, 2], method = "blomqvist"),
        copula::betan(copula::pobs(x[, 1:2]))
    )
    expect_equal(
        wdm(xt[, 1], xt[, 2], method = "blomqvist"),
        copula::betan(copula::pobs(xt[, 1:2]))
    )
    # b/c of median, calculations are slightly different for odd sn
    expect_equal(
        wdm(x[-1, 1], x[-1, 2], method = "blomqvist"),
        copula::betan(copula::pobs(x[-1, 1:2]))
    )
    expect_equal(
        wdm(xt[-1, 1], xt[-1, 2], method = "blomqvist"),
        copula::betan(copula::pobs(xt[-1, 1:2]))
    )
})
