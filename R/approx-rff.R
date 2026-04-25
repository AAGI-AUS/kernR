#' Default RFF Rank Heuristic
#'
#' Returns the default number of random frequencies used by
#' [rff_features()] when `D` is left `NULL` in [kernel_spec()].
#' Mirrors [default_nystrom_rank()]: `max(50, ceiling(0.1 * n))`,
#' capped at `n`.
#'
#' @param n Sample size.
#' @return Integer in `[1, n]`.
#' @keywords internal
default_rff_rank <- function(n) {
  d <- max(50L, as.integer(ceiling(0.1 * n)))
  min(d, as.integer(n))
}

#' Sample RFF Frequencies for a Translation-Invariant Kernel
#'
#' Draws `D` frequencies from the spectral measure associated with the
#' kernel. For RBF (and Matern with `nu = Inf`), this is
#' \eqn{\mathcal{N}(0, I_d / \sigma^2)}. For Matern with finite `nu`,
#' it is the multivariate Student-t with `2 * nu` degrees of freedom and
#' scale `1 / sigma^2` per dimension, drawn via the standard
#' Gaussian / chi-squared trick.
#'
#' @param kernel A `kernel_spec` object with resolved (numeric) bandwidth.
#' @param d Input dimension.
#' @param D Number of frequencies.
#' @param seed Optional integer seed.
#' @return A `D x d` matrix of frequencies.
#' @keywords internal
sample_rff_frequencies <- function(kernel, d, D, seed = NULL) {
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
      get(".Random.seed", envir = .GlobalEnv)
    } else {
      NULL
    }
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(as.integer(seed))
  }

  sigma <- kernel$bandwidth
  z <- matrix(stats::rnorm(D * d), D, d)

  is_rbf <- kernel$type == "rbf"
  matern_inf <- kernel$type == "matern" && (is.infinite(kernel$nu) || kernel$nu > 1e8)

  if (is_rbf || matern_inf) {
    return(z / sigma)
  }

  # Matern with finite nu: scale each row by sqrt(2*nu / u_i), u_i ~ chi^2_{2nu}
  nu <- kernel$nu
  u <- stats::rchisq(D, df = 2 * nu)
  scale <- sqrt(2 * nu / u) / sigma
  z * scale
}

#' Random Fourier Feature Map
#'
#' Computes the RFF feature map `Z` such that the rank-`(2 * D)`
#' approximation of the exact kernel matrix is \eqn{\tilde K = Z Z^\top}.
#' Storing `Z` (n x 2D) takes `O(n D)` memory and downstream computations
#' such as HSIC and MMD can in principle exploit the factorisation in
#' `O(n D^2)` time rather than `O(n^2)`.
#'
#' Bochner's theorem guarantees that for translation-invariant kernels
#' the spectral measure is non-negative; sampling `D` frequencies from it
#' and forming `(cos, sin)` features gives an unbiased Monte Carlo
#' estimator of the kernel.
#'
#' @param x Numeric matrix (n x d).
#' @param kernel A `kernel_spec` object. The kernel `type` must be
#'   `"rbf"` or `"matern"`. `bandwidth` is resolved against `x` if it is
#'   `"median"`.
#' @param D Integer number of frequencies. The resulting feature matrix
#'   has `2 * D` columns. If `NULL`, defaults to [default_rff_rank()].
#' @param seed Optional integer for reproducible frequency sampling.
#'
#' @return A list with class `"rff_features"`:
#'   \itemize{
#'     \item `Z`: n x (2 * D) feature matrix.
#'     \item `frequencies`: D x d matrix of sampled frequencies.
#'     \item `D`: number of frequencies.
#'     \item `kernel`: the kernel spec used (with bandwidth resolved).
#'   }
#'
#' @family approximations
#' @seealso [nystrom_factor()] for the Nystrom alternative, [kernel_spec()]
#'   for opt-in via the `approx` argument.
#'
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(200), 100, 2)
#' rf <- rff_features(x, kernel_spec("rbf", bandwidth = 1), D = 200, seed = 42)
#' dim(rf$Z)               # 100 x 400
#' K_hat <- rf$Z %*% t(rf$Z)
#' isSymmetric(K_hat)
#'
#' @export
rff_features <- function(x, kernel = kernel_spec("rbf"), D = NULL, seed = NULL) {
  x <- as.matrix(x)
  n <- nrow(x)
  d <- ncol(x)
  if (n == 0L) {
    stop("`x` has zero rows.", call. = FALSE)
  }

  if (!kernel$type %in% c("rbf", "matern")) {
    stop("RFF only supports translation-invariant kernels (\"rbf\", \"matern\"). ",
         "Got kernel type \"", kernel$type, "\".", call. = FALSE)
  }

  if (is.null(D)) {
    D <- default_rff_rank(n)
  }
  D <- as.integer(D)
  if (D < 1L) {
    stop("`D` must be at least 1.", call. = FALSE)
  }

  if (kernel$type %in% c("rbf", "matern") && identical(kernel$bandwidth, "median")) {
    kernel <- resolve_bandwidth(kernel, x)
  }

  omega <- sample_rff_frequencies(kernel, d, D, seed)

  # n x D matrix of inner products
  proj <- x %*% t(omega)
  scale <- 1 / sqrt(D)
  z_matrix <- cbind(cos(proj), sin(proj)) * scale

  structure(
    list(
      Z = z_matrix,
      frequencies = omega,
      D = D,
      kernel = kernel
    ),
    class = "rff_features"
  )
}

#' @export
print.rff_features <- function(x, ...) {
  cat("Random Fourier features:\n")
  cat("  n        :", nrow(x$Z), "\n")
  cat("  D        :", x$D, "\n")
  cat("  features :", ncol(x$Z), "(2 * D)\n")
  cat("  kernel   :", x$kernel$type, "\n")
  invisible(x)
}

#' Reconstruct an RFF-Approximated Kernel Matrix
#'
#' Materialises the `n x n` matrix \eqn{Z Z^\top} from an
#' [rff_features()] result. As with [nystrom_reconstruct()], most
#' downstream uses should keep the feature factor and exploit the
#' factorisation; this helper exists for parity testing and for callers
#' that need a dense matrix.
#'
#' @param features An `rff_features` object.
#' @return Symmetric n x n matrix (PSD by construction).
#' @keywords internal
rff_reconstruct <- function(features) {
  features$Z %*% t(features$Z)
}
