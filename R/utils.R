#' Validate and Coerce Input Data
#'
#' Converts input to a numeric matrix with validation.
#'
#' @param x Input data (vector, matrix, data.frame, data.table).
#' @param name Name of the argument (for error messages).
#' @param min_n Minimum number of observations.
#' @param min_d Minimum number of columns.
#'
#' @return A numeric matrix.
#' @keywords internal
validate_input <- function(x, name = "x", min_n = 1L, min_d = 1L) {
  if (is.null(x)) {
    stop("`", name, "` must not be NULL.", call. = FALSE)
  }

  if (inherits(x, "data.table")) {
    x <- as.matrix(x[, vapply(.SD, is.numeric, logical(1)), with = FALSE])
  } else if (is.data.frame(x)) {
    num_cols <- vapply(x, is.numeric, logical(1))
    if (!all(num_cols)) {
      warning("Non-numeric columns in `", name, "` dropped.", call. = FALSE)
    }
    x <- as.matrix(x[, num_cols, drop = FALSE])
  } else {
    x <- as.matrix(x)
  }

  if (!is.numeric(x)) {
    stop("`", name, "` must be numeric.", call. = FALSE)
  }

  if (any(!is.finite(x))) {
    stop("`", name, "` contains non-finite values (NA, NaN, Inf).", call. = FALSE)
  }

  if (nrow(x) < min_n) {
    stop("`", name, "` must have at least ", min_n, " observations.", call. = FALSE)
  }

  if (ncol(x) < min_d) {
    stop("`", name, "` must have at least ", min_d, " columns.", call. = FALSE)
  }

  x
}

#' Compute Effective Sample Size
#'
#' Computes ESS from importance weights: ESS = (sum(w))^2 / sum(w^2).
#'
#' @param w Numeric vector of weights.
#'
#' @return Scalar effective sample size.
#'
#' @family diagnostics
#' @seealso [assess_overlap()], [plot_weights()].
#'
#' @examples
#' w <- runif(100, 0.5, 2)
#' effective_sample_size(w)
#'
#' @export
effective_sample_size <- function(w) {
  if (any(w < 0)) stop("Weights must be non-negative.", call. = FALSE)
  sw <- sum(w)
  if (sw == 0) return(0)
  sw^2 / sum(w^2)
}
