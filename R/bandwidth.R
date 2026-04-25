#' Select Kernel Bandwidth
#'
#' Computes a bandwidth (lengthscale) for RBF or Matern kernels using
#' the specified method.
#'
#' @param x Numeric matrix or vector.
#' @param method Character. `"median"` (default), `"scott"`, or a
#'   positive number for fixed bandwidth.
#'
#' @return A positive numeric scalar.
#'
#' @family kernels
#' @seealso [kernel_spec()], [kernel_matrix()].
#'
#' @details
#' - `"median"`: The median heuristic sets bandwidth to the square root
#'   of the median of pairwise squared distances. Robust default for most
#'   kernel tests (Gretton et al., 2012).
#' - `"scott"`: Scott's rule: `n^(-1/(d+4)) * sd_pooled`. Good for
#'   density estimation but may undersmooth for testing.
#'
#' @examples
#' x <- matrix(rnorm(200), 100, 2)
#' select_bandwidth(x, "median")
#'
#' @export
select_bandwidth <- function(x, method = "median") {
  x <- as.matrix(x)
  n <- nrow(x)
  d <- ncol(x)

  if (is.numeric(method)) {
    if (method <= 0) stop("`method` must be positive when numeric.", call. = FALSE)
    return(method)
  }

  method <- match.arg(method, c("median", "scott"))

  switch(method,
    median = median_bandwidth_cpp(x),
    scott = {
      sd_pooled <- sqrt(mean(apply(x, 2, var)))
      n^(-1 / (d + 4)) * sd_pooled
    }
  )
}
