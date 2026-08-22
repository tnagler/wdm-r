test_that("Chatterjee's xi is directional", {
    x <- seq_len(20)
    y <- (x - 10)^2

    expect_false(isTRUE(all.equal(
        wdm(x, y, method = "chatterjee"),
        wdm(y, x, method = "chatterjee")
    )))

    estimates <- wdm(cbind(x, y), method = "chatterjee")
    expect_equal(estimates[1, 2], wdm(x, y, method = "chatterjee"))
    expect_equal(estimates[2, 1], wdm(y, x, method = "chatterjee"))
})

test_that("predictor-tie ordering is reproducible and configurable", {
    x <- rep(1:2, each = 30)
    y <- (37 * seq_along(x)) %% length(x)

    expect_identical(
        wdm(x, y, method = "chatterjee"),
        wdm(x, y, method = "chatterjee")
    )
    expect_gt(length(unique(vapply(seq_len(6), function(seed) {
        wdm(x, y, method = "chatterjee", seeds = seed)
    }, numeric(1)))), 1)
})

test_that("Chatterjee independence controls reach the backend", {
    x <- seq_len(20)
    y <- c(2:20, 1)

    expect_true(is.finite(indep_test(
        x, y, method = "chatterjee", alternative = "greater"
    )$p_value))
    expect_true(is.finite(indep_test(
        x, y, method = "chatterjee", y_continuous = FALSE
    )$p_value))
    expect_error(indep_test(
        x, y, method = "chatterjee", weights = seq_along(x),
        y_continuous = FALSE
    ), "discrete or tied response")
})

test_that("Chatterjee arguments are validated", {
    expect_error(
        wdm(1:5, 1:5, method = "chatterjee", seeds = 1.5),
        "finite integers"
    )
    expect_error(
        indep_test(1:5, 1:5, method = "chatterjee", y_continuous = NA),
        "TRUE or FALSE"
    )
})
