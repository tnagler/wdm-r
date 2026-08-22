#' Independence Tests for Weighted Dependence Measures
#'
#' Computes a dependence measure and its asymptotic independence test for two
#' numeric vectors.
#'
#' @param x,y numeric vectors of data values. `x` and `y` must have the same
#'   length.
#' @param method the dependence measure; see *Details* for possible values.
#' @param weights an optional vector of weights for the observations.
#' @param remove_missing if `TRUE`, all pairwise incomplete observations are
#'   removed; if `FALSE`, the function throws an error if there are incomplete
#'   observations.
#' @param alternative indicates the alternative hypothesis and must be one of
#'   `"two-sided"`, `"greater"` or `"less"`. You can specify just the initial
#'   letter. `"greater"` corresponds to positive association, `"less"` to
#'   negative association.
#' @param seeds an optional integer vector used to break predictor ties for
#'   Chatterjee's xi. The default uses a fixed, reproducible ordering.
#' @param y_continuous whether the response distribution is known to be
#'   continuous for Chatterjee inference. Set this to `FALSE` for a discrete
#'   response even if the sample contains no observed response ties.
#'
#' @details Available methods:
#' - `"pearson"`: Pearson correlation
#' - `"spearman"`: Spearman's \eqn{\rho}
#' - `"kendall"`: Kendall's \eqn{\tau}
#' - `"blomqvist"`: Blomqvist's \eqn{\beta}
#' - `"hoeffding"`: Hoeffding's \eqn{D}
#' - `"chatterjee"`: Chatterjee's \eqn{\xi}
#'
#' Partial matching of method names is enabled.
#' This implementation of Hoeffding's \eqn{D} does not support tied
#' observations; test results are invalid when ties are present. It supports
#' only the two-sided alternative. The natural one-sided alternative for
#' Chatterjee's \eqn{\xi} is `"greater"`.
#'
#' Chatterjee's \eqn{\xi} measures the dependence of `y` on `x`. Analytic
#' inference with unequal weights requires a continuous response and assumes
#' that weights are fixed or depend only on `x`. It is unavailable when the
#' response is discrete or tied and weights are unequal.
#'
#' @return A one-row data frame containing the estimate, transformed test
#'   statistic, p-value, effective sample size, method, and alternative.
#'
#' @export
#'
#' @examples
#' x <- rnorm(100)
#' y <- rpois(100, 1)
#' w <- runif(100)
#'
#' indep_test(x, y, method = "kendall")               # unweighted
#' indep_test(x, y, method = "kendall", weights = w)  # weighted
#'
indep_test <- function(x, y, method = "pearson", weights = NULL,
                       remove_missing = TRUE, alternative = "two-sided",
                       seeds = NULL, y_continuous = TRUE) {
    if (is.null(weights))
        weights <- numeric(0)
    check_indep_test_inputs(x, y, weights, remove_missing)
    method <- match.arg(method, allowed_methods)
    alternative <- match.arg(alternative, allowed_alternatives)
    seeds <- normalize_seeds(seeds)
    check_logical_scalar(y_continuous, "y_continuous")

    test <- indep_test_cpp(x, y, method, weights, remove_missing, alternative,
                           seeds, y_continuous)
    as.data.frame(test)
}

check_indep_test_inputs <- function(x, y, weights, remove_missing) {
    if (!(is.numeric(x) || is.logical(x)) || !is.null(dim(x)))
        stop("'x' must be a numeric vector")
    if (!(is.numeric(y) || is.logical(y)) || !is.null(dim(y)))
        stop("'y' must be a numeric vector")
    if (length(x) != length(y))
        stop("'x' and 'y' must have the same length")
    if (!is.numeric(weights) || !is.null(dim(weights)))
        stop("'weights' must be a numeric vector")
    if (length(weights) > 0 && length(weights) != length(x))
        stop("'weights' must have one value per observation")
    check_logical_scalar(remove_missing, "remove_missing")
}
