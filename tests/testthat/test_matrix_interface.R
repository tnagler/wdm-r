test_that("vectors and matrices can be combined", {
    x <- seq_len(10)
    y <- cbind(increasing = x, decreasing = rev(x))

    expect_equal(wdm(x, y), cor(x, y))
    expect_equal(wdm(y, x), cor(y, x))
})

test_that("one-column matrices return one-by-one matrices", {
    x <- matrix(seq_len(10), ncol = 1, dimnames = list(NULL, "x"))

    expect_equal(wdm(x), cor(x))
})

test_that("undefined matrix diagonals are not forced to one", {
    x <- cbind(constant = rep(1, 10), increasing = seq_len(10))

    expect_true(is.na(wdm(x, method = "pearson")["constant", "constant"]))
    expect_true(is.na(wdm(x, method = "kendall")["constant", "constant"]))
    expect_true(is.na(wdm(x, method = "spearman")["constant", "constant"]))
})
