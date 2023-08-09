context("weighted ranks")

test_that("weight 1 corresponds to unweighted", {
    x <- rnorm(50)

    # x[c(2, 10, 20, 30, 50)] <- NA
    x <- c(x, rep(Inf, 10))
    w <- rep(1, 60)

    expect_equal(
        wdm:::rank_wtd(x, w, "average"),
        rank(x, ties = "average"),
    )

    expect_equal(
        wdm:::rank_wtd(x, w, "min"),
        rank(x, ties = "min"),
    )
})
