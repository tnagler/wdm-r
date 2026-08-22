#' Computing weighted ranks
#'
#' The weighted rank of \eqn{X_i} among \eqn{X_1, \dots, X_n} with weights
#' \eqn{w_1, \dots, w_n} is defined as
#' \deqn{\frac 1 n \sum_{j = 1}^n w_i 1[X_j \le X_i].}
#'
#' @param x a numeric vector.
#' @param weights a vector of weights (same length as `x`).
#' @param ties_method Indicates how to treat ties; same as in R, see
#'  https://stat.ethz.ch/R-manual/R-devel/library/base/html/rank.html.
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
