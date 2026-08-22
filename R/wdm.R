#' Weighted Dependence Measures
#'
#' Computes a (possibly weighted) dependence measure between `x` and `y` if
#' these are vectors. If either argument is a matrix, the measures between all
#' corresponding columns are computed.
#'
#' @param x a numeric vector, matrix or data frame.
#' @param y `NULL` (default) or a vector, matrix or data frame with compatible
#'   dimensions to `x`. The default is equivalent to `y = x` (but more
#'   efficient).
#' @param method the dependence measure; see *Details* for possible values.
#' @param weights an optional vector of weights for the observations.
#' @param remove_missing if `TRUE`, all pairwise incomplete observations are
#'   removed; if `FALSE`, the function throws an error if there are incomplete
#'   observations.
#' @param seeds an optional integer vector used to break predictor ties for
#'   Chatterjee's xi. The default uses a fixed, reproducible ordering.
#'
#' @details Available methods:
#' - `"pearson"`: Pearson correlation
#' - `"spearman"`: Spearman's \eqn{\rho}
#' - `"kendall"`: Kendall's \eqn{\tau}
#' - `"blomqvist"`: Blomqvist's \eqn{\beta}
#' - `"hoeffding"`: Hoeffding's \eqn{D}
#' - `"chatterjee"`: Chatterjee's \eqn{\xi}
#' Partial matching of method names is enabled.
#'
#' Spearman's \eqn{\rho} and Kendall's \eqn{\tau} are corrected for ties if
#' there are any.
#' This implementation of Hoeffding's \eqn{D} does not support tied
#' observations; estimates are invalid when ties are present.
#' Chatterjee's \eqn{\xi} measures the dependence of `y` on `x` and is
#' generally asymmetric. Consequently, `wdm(x, method = "chatterjee")` need
#' not return a symmetric matrix.
#'
#' @return A numeric scalar when both inputs are vectors or one-column objects;
#'   otherwise, a matrix containing the dependence measure for every pair of
#'   columns.
#'
#' @export
#'
#' @examples
#' ##  dependence between two vectors
#' x <- rnorm(100)
#' y <- rpois(100, 1)
#' w <- runif(100)
#' wdm(x, y, method = "kendall")               # unweighted
#' wdm(x, y, method = "kendall", weights = w)  # weighted
#'
#' ##  dependence in a matrix
#' x <- matrix(rnorm(100 * 3), 100, 3)
#' wdm(x, method = "spearman")               # unweighted
#' wdm(x, method = "spearman", weights = w)  # weighted
#'
#' ##  dependence between columns of two matrices
#' y <- matrix(rnorm(100 * 2), 100, 2)
#' wdm(x, y, method = "hoeffding")               # unweighted
#' wdm(x, y, method = "hoeffding", weights = w)  # weighted
#'
wdm <- function(x, y = NULL, method = "pearson", weights = NULL,
                remove_missing = TRUE, seeds = NULL) {
    ## preprocessing of arguments
    if (is.null(weights))
        weights <- numeric(0)
    if (is.data.frame(y))
        y <- as.matrix(y)
    if (is.data.frame(x))
        x <- as.matrix(x)
    check_wdm_inputs(x, y, weights, remove_missing)
    method <- match.arg(method, allowed_methods)
    seeds <- normalize_seeds(seeds)

    ## computations
    if (is.null(y)) {
        out <- wdm_mat_cpp(x, method, weights, remove_missing, seeds)
        colnames(out) <- rownames(out) <- colnames(x)
    } else if (NCOL(x) == 1 && NCOL(y) == 1) {
        out <- wdm_cpp(x, y, method, weights, remove_missing, seeds)
    } else {
        out <- matrix(NA, NCOL(x), NCOL(y))
        for (i in seq_len(NCOL(x))) {
            for (j in seq_len(NCOL(y))) {
                out[i, j] <- wdm(
                    if (is.matrix(x)) x[, i] else x,
                    if (is.matrix(y)) y[, j] else y,
                    method, weights, remove_missing, seeds
                )
            }
        }
        rownames(out) <- colnames(x)
        colnames(out) <- colnames(y)
    }

    out[is.nan(out)] <- NA
    out
}

normalize_seeds <- function(seeds) {
    if (is.null(seeds))
        return(integer(0))
    if (!is.numeric(seeds) || !is.null(dim(seeds)) || anyNA(seeds) ||
        any(!is.finite(seeds)) ||
        any(seeds != trunc(seeds)) ||
        any(abs(seeds) > .Machine$integer.max))
        stop("'seeds' must contain finite integers")
    as.integer(seeds)
}

check_wdm_inputs <- function(x, y, weights, remove_missing) {
    if (!(is.numeric(x) || is.logical(x)))
        stop("'x' must be numeric")
    if (!is.null(dim(x)) && length(dim(x)) != 2)
        stop("'x' must be a vector or matrix")
    if (!is.matrix(x) && is.null(y))
        stop("supply both 'x' and 'y' or a matrix-like 'x'")
    if (!is.numeric(weights) || !is.null(dim(weights)))
        stop("'weights' must be a numeric vector")
    if (length(weights) > 0 && length(weights) != NROW(x))
        stop("'weights' must have one value per observation")
    if (!is.null(y)) {
        if (!(is.numeric(y) || is.logical(y)))
            stop("'y' must be numeric")
        if (!is.null(dim(y)) && length(dim(y)) != 2)
            stop("'y' must be a vector or matrix")
        if (NROW(x) != NROW(y))
            stop("'x' and 'y' must have the same number of rows")
    }
    check_logical_scalar(remove_missing, "remove_missing")
}

check_logical_scalar <- function(value, name) {
    if (!is.logical(value) || length(value) != 1 || is.na(value))
        stop(sprintf("'%s' must be TRUE or FALSE", name))
}
