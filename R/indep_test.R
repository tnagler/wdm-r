#' Independence Tests for Weighted Dependence Measures
#'
#' Computes a (possibly weighted) dependence measure between `x` and `y` if
#' these are vectors. If `x` and `y` are matrices then the measure between the
#' columns of `x` and the columns of `y` are computed.
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
#' @export
#'
#' @examples
#' x <- rnorm(100)
#' y <- rpois(100, 1)  # all but Hoeffding's D can handle ties
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
    if (!is.logical(y_continuous) || length(y_continuous) != 1 ||
        is.na(y_continuous))
        stop("'y_continuous' must be TRUE or FALSE")

    test <- indep_test_cpp(x, y, method, weights, remove_missing, alternative,
                           seeds, y_continuous)
    as.data.frame(test)
}

allowed_methods <- c("pearson", "kendall", "spearman", "hoeffding",
                     "blomqvist", "chatterjee")
allowed_alternatives <- c("two-sided", "less", "greater")

check_indep_test_inputs <- function(x, y, weights, remove_missing) {
    if (!(is.numeric(x) | is.logical(x)) | (NCOL(x) != 1) )
        stop("'x' must be a numeric vector")
    if (!(is.numeric(y) | is.logical(y)) | (NCOL(y) != 1) )
        stop("'y' must be a numeric vector")
    if (!is.numeric(weights) | (NCOL(weights) != 1))
        stop("'weights' must be a numeric vector")
    stopifnot(is.atomic(x))
    stopifnot(is.atomic(y))
    stopifnot(is.atomic(weights))
    if (!is.logical(remove_missing))
        stop("remove_missing must be logical.")
}
