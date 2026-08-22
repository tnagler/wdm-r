#' Computing weighted ranks
#'
#' For observations without ties, the weighted rank of \eqn{X_i} among
#' \eqn{X_1, \dots, X_n} with weights \eqn{w_1, \dots, w_n} is
#' \deqn{\frac{n}{\sum_{k = 1}^n w_k}
#'       \sum_{j = 1}^n w_j 1[X_j \le X_i].}
#' Thus, multiplying every weight by the same positive constant does not change
#' the ranks, and unit weights reproduce ordinary ranks. Tied values are
#' handled according to `ties_method`.
#'
#' @param x a numeric vector.
#' @param weights an optional vector of nonnegative weights with the same
#'   length as `x`.
#' @param ties_method how to treat ties; one of `"average"`, `"min"`,
#'   `"first"`, or `"random"`, as in [rank()].
#'
#' @return a vector of ranks.
#' @export
#'
#' @examples
#' x <- rnorm(100)
#' w <- rexp(100)
#' rank(x)
#' rank_wtd(x, w)
rank_wtd <- function(x, weights = numeric(), ties_method = "average") {
    ## preprocessing of arguments
    if (!is.numeric(x) || !is.null(dim(x)))
        stop("'x' must be a numeric vector")
    if (!is.numeric(weights) || !is.null(dim(weights)))
        stop("'weights' must be a numeric vector")
    if (length(weights) > 0 && length(weights) != length(x))
        stop("'weights' must have the same length as 'x'")
    if (length(weights) > 0 &&
        (any(!is.finite(weights)) || any(weights < 0) || sum(weights) <= 0))
        stop("'weights' must be finite, nonnegative, and have a positive sum")
    ties_method <- match.arg(
        ties_method, c("average", "min", "first", "random")
    )
    rank_wtd_cpp(x, weights, ties_method)
}
