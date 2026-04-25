#' Compute a Kernel Matrix
#'
#' Computes the kernel (Gram) matrix between two sets of observations.
#'
#' @param x Numeric matrix (n x d) or vector.
#' @param y Numeric matrix (m x d) or vector. If `NULL` (default),
#'   computes the kernel matrix of `x` with itself.
#' @param kernel A `kernel_spec` object. Default is RBF with median
#'   heuristic bandwidth. If the spec sets `approx = "nystrom"` and the
#'   call is symmetric (`y` is `NULL`), the matrix is returned via the
#'   Nystrom approximation; rectangular calls always use exact
#'   computation.
#'
#' @return An n x m numeric matrix.
#'
#' @family kernels
#' @seealso [kernel_spec()], [nystrom_factor()], [rff_features()].
#'
#' @examples
#' x <- matrix(rnorm(100), 50, 2)
#' K <- kernel_matrix(x)
#' dim(K)  # 50 x 50
#'
#' # Nystrom-approximated symmetric Gram matrix
#' k_nys <- kernel_spec("rbf", approx = "nystrom", approx_rank = 20, approx_seed = 1)
#' K_nys <- kernel_matrix(x, kernel = k_nys)
#' max(abs(K - K_nys))  # small for moderate rank
#'
#' # RFF-approximated kernel
#' k_rff <- kernel_spec("rbf", approx = "rff", approx_rank = 200, approx_seed = 1)
#' K_rff <- kernel_matrix(x, kernel = k_rff)
#' max(abs(K - K_rff))
#'
#' @export
kernel_matrix <- function(x, y = NULL, kernel = kernel_spec()) {
  x <- as.matrix(x)
  symmetric <- is.null(y)
  if (symmetric) y <- x else y <- as.matrix(y)

  if (ncol(x) != ncol(y)) {
    stop("`x` and `y` must have the same number of columns.", call. = FALSE)
  }

  # Resolve bandwidth if needed
  if (kernel$type %in% c("rbf", "matern") && identical(kernel$bandwidth, "median")) {
    kernel <- resolve_bandwidth(kernel, if (symmetric) x else rbind(x, y))
  }

  approx <- if (is.null(kernel$approx)) "none" else kernel$approx
  if (approx == "nystrom" && symmetric) {
    nf <- nystrom_factor(
      x,
      kernel = kernel,
      m = kernel$approx_rank,
      seed = kernel$approx_seed
    )
    return(nystrom_reconstruct(nf))
  }
  if (approx == "rff" && symmetric) {
    rf <- rff_features(
      x,
      kernel = kernel,
      D = kernel$approx_rank,
      seed = kernel$approx_seed
    )
    return(rff_reconstruct(rf))
  }

  exact_kernel_block(x, y, kernel)
}

#' Centre a Kernel Matrix
#'
#' Centres a kernel matrix in feature space (double centring).
#'
#' @param K Square numeric matrix.
#'
#' @return Centred kernel matrix.
#' @keywords internal
centre_kernel_matrix <- function(K) {
  n <- nrow(K)
  H <- diag(n) - 1 / n
  H %*% K %*% H
}
