#' Default Nystrom Rank Heuristic
#'
#' Returns the default number of landmark points used by the Nystrom
#' approximation when `approx_rank` is left `NULL` in [kernel_spec()].
#' Defaults to `max(50, ceiling(0.1 * n))`, capped at `n`.
#'
#' @param n Sample size.
#' @return Integer rank in `[1, n]`.
#' @keywords internal
default_nystrom_rank <- function(n) {
  m <- max(50L, as.integer(ceiling(0.1 * n)))
  min(m, as.integer(n))
}

#' Sample Landmark Indices for Nystrom
#'
#' Draws `m` indices uniformly without replacement from `seq_len(n)`.
#'
#' @param n Sample size.
#' @param m Landmark count.
#' @param seed Optional integer seed.
#' @return Integer vector of length `m`.
#' @keywords internal
sample_landmarks <- function(n, m, seed = NULL) {
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
  sample.int(n, size = m, replace = FALSE)
}

#' Nystrom Low-Rank Factor for a Kernel
#'
#' Computes the Nystrom factor `U` such that the approximation
#' \eqn{\tilde K = U U^\top} is a rank-r (with \eqn{r \le m}) approximation
#' of the exact kernel matrix `K = kernel_matrix(x, x, kernel)`.
#'
#' This is the cheap representation: storing `U` (n x r) takes `O(nr)`
#' memory, and downstream HSIC/MMD computations can exploit the
#' factorisation in `O(n r^2)` time rather than `O(n^2)`. Reconstructing
#' `tilde_K = U %*% t(U)` collapses back to an `n x n` matrix and costs
#' `O(n^2 r)`.
#'
#' @param x Numeric matrix (n x d).
#' @param kernel A `kernel_spec` object. Bandwidth is resolved against
#'   `x` if it is `"median"`.
#' @param m Integer landmark count. If `NULL`, defaults to
#'   [default_nystrom_rank()].
#' @param seed Optional integer for reproducible landmark sampling.
#' @param tol Eigenvalue threshold (relative to `max(eigval(W))`).
#'   Eigenvalues at or below `tol * max(eigval(W))` are dropped, so the
#'   returned factor may have fewer than `m` columns.
#'
#' @return A list with class `"nystrom_factor"`:
#'   \itemize{
#'     \item `U`: n x r factor matrix.
#'     \item `landmarks`: integer vector of length m (sampled indices).
#'     \item `m`: requested landmark count.
#'     \item `rank`: effective rank `r` (columns of `U`).
#'     \item `kernel`: the kernel spec used (with bandwidth resolved).
#'   }
#'
#' @family approximations
#' @seealso [rff_features()] for the Random Fourier Features alternative,
#'   [kernel_spec()] for opt-in via the `approx` argument.
#'
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(200), 100, 2)
#' nf <- nystrom_factor(x, kernel_spec("rbf", bandwidth = 1), m = 30, seed = 42)
#' dim(nf$U)              # 100 x rank
#' K_hat <- nf$U %*% t(nf$U)
#' isSymmetric(K_hat)
#'
#' @export
nystrom_factor <- function(x, kernel = kernel_spec(), m = NULL,
                           seed = NULL, tol = 1e-10) {
  x <- as.matrix(x)
  n <- nrow(x)
  if (n == 0L) {
    stop("`x` has zero rows.", call. = FALSE)
  }

  if (is.null(m)) {
    m <- default_nystrom_rank(n)
  }
  m <- as.integer(m)
  if (m < 1L) {
    stop("`m` must be at least 1.", call. = FALSE)
  }
  if (m > n) {
    stop("`m` cannot exceed nrow(x); set `approx = \"none\"` instead.",
         call. = FALSE)
  }
  if (!is.numeric(tol) || tol < 0) {
    stop("`tol` must be a non-negative numeric.", call. = FALSE)
  }

  if (kernel$type %in% c("rbf", "matern") && identical(kernel$bandwidth, "median")) {
    kernel <- resolve_bandwidth(kernel, x)
  }

  idx <- sample_landmarks(n, m, seed)
  x_land <- x[idx, , drop = FALSE]

  c_block <- exact_kernel_block(x, x_land, kernel)
  w_block <- exact_kernel_block(x_land, x_land, kernel)

  u <- nystrom_factor_cpp(c_block, w_block, tol)

  structure(
    list(
      U = u,
      landmarks = idx,
      m = m,
      rank = ncol(u),
      kernel = kernel
    ),
    class = "nystrom_factor"
  )
}

#' @export
print.nystrom_factor <- function(x, ...) {
  cat("Nystrom factor:\n")
  cat("  n        :", nrow(x$U), "\n")
  cat("  m        :", x$m, "\n")
  cat("  rank     :", x$rank, "\n")
  cat("  kernel   :", x$kernel$type, "\n")
  invisible(x)
}

#' Reconstruct a Nystrom-Approximated Kernel Matrix
#'
#' Materialises the `n x n` matrix \eqn{U U^\top} from a
#' [nystrom_factor()] result. For most downstream uses you should keep
#' the factor and exploit the factorisation directly; this helper exists
#' for parity testing and for callers that need a dense matrix.
#'
#' @param factor A `nystrom_factor` object.
#' @return Symmetric n x n matrix (PSD by construction).
#' @keywords internal
nystrom_reconstruct <- function(factor) {
  factor$U %*% t(factor$U)
}

#' Internal: Exact Kernel Block (no approx dispatch)
#'
#' Computes a kernel matrix block via the same primitives as
#' [kernel_matrix()] but bypasses the approx dispatch. Used by Nystrom
#' to evaluate the C and W blocks without recursing into approximation
#' logic.
#'
#' @keywords internal
#' @noRd
exact_kernel_block <- function(x, y, kernel) {
  switch(
    kernel$type,
    rbf = rbf_kernel_matrix_cpp(x, y, kernel$bandwidth),
    matern = matern_kernel_matrix_cpp(x, y, kernel$bandwidth, kernel$nu),
    linear = linear_kernel_matrix_cpp(x, y),
    polynomial = polynomial_kernel_matrix_cpp(x, y, kernel$degree, kernel$offset),
    stop("Unknown kernel type: ", kernel$type, call. = FALSE)
  )
}
