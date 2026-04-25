#' Create a Kernel Specification
#'
#' Constructs a kernel specification object used throughout `kernR` for
#' computing kernel matrices. Supports RBF (Gaussian), Matern, linear,
#' and polynomial kernels, with optional low-rank approximation.
#'
#' @param type Character. Kernel type: `"rbf"` (default), `"matern"`,
#'   `"linear"`, or `"polynomial"`.
#' @param bandwidth Numeric or `"median"`. Lengthscale parameter for RBF
#'   and Matern kernels. If `"median"` (default), the median heuristic is
#'   used to select bandwidth automatically from the data.
#' @param nu Numeric. Smoothness parameter for the Matern kernel. Common
#'   choices: 0.5 (Laplace), 1.5, 2.5, Inf (RBF). Default is 2.5.
#' @param degree Integer. Degree for polynomial kernel. Default is 2.
#' @param offset Numeric. Offset for polynomial kernel. Default is 1.
#' @param approx Character. Approximation strategy:
#'   `"none"` (default, exact computation), `"nystrom"` (low-rank
#'   Nystrom approximation; symmetric Gram matrices only), or `"rff"`
#'   (Random Fourier Features; translation-invariant kernels only —
#'   currently `"rbf"` and `"matern"`, the latter via Student-t spectral
#'   sampling). The approximation is opt-in: callers must explicitly
#'   set this argument.
#' @param approx_rank Integer or `NULL`. For `"nystrom"`, number of
#'   landmark points (m). For `"rff"`, number of random frequencies (D);
#'   the resulting feature map has dimension `2 * D` (cos and sin
#'   blocks), so the approximated kernel has rank at most `2 * D`. If
#'   `NULL`, defaults to `max(50, ceiling(0.1 * n))`, capped at `n`.
#'   Larger values give a more accurate approximation at the cost of
#'   compute and memory.
#' @param approx_seed Integer or `NULL`. Seed for reproducible landmark
#'   sampling. Set this when you need bitwise-identical Nystrom output
#'   across runs.
#'
#' @return An object of class `"kernel_spec"`.
#'
#' @family kernels
#' @seealso [kernel_matrix()], [select_bandwidth()], [nystrom_factor()],
#'   [rff_features()].
#'
#' @details
#' The `approx` argument is honoured by [kernel_matrix()] for symmetric
#' Gram matrices (i.e. `y` is `NULL` or identical to `x`); rectangular
#' calls always compute the exact block. To compute and retain the
#' low-rank factor directly without materialising the n x n matrix,
#' use [nystrom_factor()] or [rff_features()] depending on the chosen
#' approximation.
#'
#' RFF is restricted to translation-invariant kernels: `"rbf"` and
#' `"matern"`. For Matern with finite smoothness `nu`, frequencies are
#' drawn from the multivariate Student-t spectral measure with `2 * nu`
#' degrees of freedom; for `nu = Inf` (or RBF), Gaussian frequencies.
#' Linear and polynomial kernels are not translation-invariant and will
#' error at compute time.
#'
#' @examples
#' # Default RBF kernel with median heuristic bandwidth
#' k <- kernel_spec()
#'
#' # RBF with fixed bandwidth
#' k <- kernel_spec("rbf", bandwidth = 1.0)
#'
#' # Matern kernel
#' k <- kernel_spec("matern", nu = 1.5)
#'
#' # Linear kernel (no bandwidth needed)
#' k <- kernel_spec("linear")
#'
#' # Nystrom-approximated RBF kernel with 100 landmarks
#' k <- kernel_spec("rbf", approx = "nystrom", approx_rank = 100, approx_seed = 1)
#'
#' # Random Fourier Features for an RBF kernel with 200 frequencies
#' k <- kernel_spec("rbf", approx = "rff", approx_rank = 200, approx_seed = 1)
#'
#' @export
kernel_spec <- function(type = c("rbf", "matern", "linear", "polynomial"),
                        bandwidth = "median",
                        nu = 2.5,
                        degree = 2L,
                        offset = 1.0,
                        approx = c("none", "nystrom", "rff"),
                        approx_rank = NULL,
                        approx_seed = NULL) {
  type <- match.arg(type)
  approx <- match.arg(approx)

  if (type %in% c("rbf", "matern")) {
    if (!identical(bandwidth, "median") && !is.numeric(bandwidth)) {
      stop("`bandwidth` must be numeric or \"median\".", call. = FALSE)
    }
    if (is.numeric(bandwidth) && bandwidth <= 0) {
      stop("`bandwidth` must be positive.", call. = FALSE)
    }
  }

  if (type == "matern" && (!is.numeric(nu) || nu <= 0)) {
    stop("`nu` must be a positive number.", call. = FALSE)
  }

  if (!is.null(approx_rank)) {
    if (!is.numeric(approx_rank) || length(approx_rank) != 1L || approx_rank < 1) {
      stop("`approx_rank` must be a positive integer or NULL.", call. = FALSE)
    }
    approx_rank <- as.integer(approx_rank)
  }
  if (!is.null(approx_seed)) {
    if (!is.numeric(approx_seed) || length(approx_seed) != 1L) {
      stop("`approx_seed` must be a single integer or NULL.", call. = FALSE)
    }
    approx_seed <- as.integer(approx_seed)
  }

  structure(
    list(
      type = type,
      bandwidth = bandwidth,
      nu = nu,
      degree = as.integer(degree),
      offset = offset,
      approx = approx,
      approx_rank = approx_rank,
      approx_seed = approx_seed
    ),
    class = "kernel_spec"
  )
}

#' @export
print.kernel_spec <- function(x, ...) {
  cat("Kernel specification:\n")
  cat("  Type:", x$type, "\n")
  if (x$type %in% c("rbf", "matern")) {
    bw <- if (identical(x$bandwidth, "median")) "median heuristic" else x$bandwidth
    cat("  Bandwidth:", bw, "\n")
  }
  if (x$type == "matern") cat("  nu:", x$nu, "\n")
  if (x$type == "polynomial") {
    cat("  Degree:", x$degree, "\n")
    cat("  Offset:", x$offset, "\n")
  }
  if (!is.null(x$approx) && x$approx != "none") {
    cat("  Approx:", x$approx, "\n")
    rank_str <- if (is.null(x$approx_rank)) "auto (max(50, 0.1n))" else as.character(x$approx_rank)
    cat("  Approx rank:", rank_str, "\n")
    if (!is.null(x$approx_seed)) cat("  Approx seed:", x$approx_seed, "\n")
  }
  invisible(x)
}

#' Resolve Kernel Bandwidth
#'
#' If `kernel$bandwidth` is `"median"`, compute the median heuristic from
#' the data. Otherwise return the fixed bandwidth.
#'
#' @param kernel A `kernel_spec` object.
#' @param x Numeric matrix (n x d).
#'
#' @return A `kernel_spec` with resolved numeric bandwidth.
#' @keywords internal
resolve_bandwidth <- function(kernel, x) {
  if (!identical(kernel$bandwidth, "median")) return(kernel)

  x <- as.matrix(x)
  bw <- median_bandwidth_cpp(x)
  kernel$bandwidth <- bw
  kernel
}
