test_that("logical flags are scalar and non-missing", {
    expect_error(wdm(1:5, 1:5, remove_missing = NA), "TRUE or FALSE")
    expect_error(
        wdm(1:5, 1:5, remove_missing = c(TRUE, FALSE)), "TRUE or FALSE"
    )
    expect_error(
        indep_test(1:5, 1:5, remove_missing = NA), "TRUE or FALSE"
    )
    expect_error(
        indep_test(1:5, 1:5, y_continuous = c(TRUE, FALSE)), "TRUE or FALSE"
    )
})

test_that("dimensions and lengths are validated in R", {
    expect_error(
        wdm(array(1:8, dim = c(2, 2, 2))), "vector or matrix"
    )
    expect_error(
        wdm(1:5, 1:5, weights = matrix(rep(1, 5))), "numeric vector"
    )
    expect_error(
        wdm(1:5, 1:5, weights = 1:4), "one value per observation"
    )
    expect_error(
        indep_test(1:5, 1:4), "same length"
    )
    expect_error(
        indep_test(1:5, 1:5, weights = 1:4), "one value per observation"
    )
})

test_that("weighted-rank arguments are validated in R", {
    expect_error(rank_wtd(matrix(1:4, ncol = 1)), "numeric vector")
    expect_error(rank_wtd(1:5, 1:4), "same length")
    expect_error(rank_wtd(1:5, c(1, 1, 1, 1, -1)), "nonnegative")
    expect_error(rank_wtd(1:5, rep(0, 5)), "positive sum")
    expect_error(rank_wtd(1:5, rep(1, 5), "unknown"))
})
