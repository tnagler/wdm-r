#' Independence Tests for Weighted Dependence Measures
#'
#' Computes a (possibly weighted) dependence measure between `x` and `y` if
#' these are vectors. If `x` and `y` are matrices then the measure between the
#' columns of `x` and the columns of `y` are computed.
#'
#' @param x, y numeric vectors of data values. `x` and `y` must have the same
#'   length.
#' @param method the dependence measure; see *Details* for possible values.
#' @param weights an optional vector of weights for the observations.
#' @param remove_missing if `TRUE`, all incomplete observations containing re
#'   removed; otherwise throws an error if there are any.
#'
#' @details Available methods: - `"pearson"`: Pearson correlation -
#' `"spearman"`: Spearman's \eqn{\rho} - `"kendall"`: Kendall's \eqn{\tau} -
#' `"blomqvist"`: Blomqvist's \eqn{\beta} - `"hoeffding"`: Hoeffding's \eqn{D}
#' Partial matching of method names is enabled.
#'
#' All methods besides `"hoeffding` work with discrete variables.
#'
#' @export
#'
#' @examples
indep_test <- function(x, y, method = "pearson", weights = NULL,
                       remove_missing = TRUE, alternative = "two-sided") {
    if (is.null(weights))
        weights <- numeric(0)
    check_indep_test_inputs(x, y, weights, remove_missing)
    method <- match.arg(method, allowed_methods)
    alternative <- match.arg(alternative, allowed_alternatives)

    test <- indep_test_cpp(x, y, method, weights, remove_missing, alternative)
    as.data.frame(test)
}

allowed_methods <- c("pearson", "kendall", "spearman", "hoeffding", "blomqvist")
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