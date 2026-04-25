#' MMD Two-Sample Test
#'
#' Tests whether two samples come from the same distribution using the
#' Maximum Mean Discrepancy (MMD). Uses a permutation test for inference.
#'
#' @param x Numeric vector, matrix, or data.frame. First sample.
#' @param y Numeric vector, matrix, or data.frame. Second sample.
#' @param kernel Kernel specification. Default is RBF with median heuristic.
#' @param n_permutations Integer. Number of permutations when
#'   `adaptive = FALSE`. Default is 500.
#' @param alpha Numeric. Significance level. Used both as the reporting
#'   threshold and, when `adaptive = TRUE`, as the early-stopping
#'   threshold. Default is 0.05.
#' @param adaptive Logical. If `TRUE`, use the sequential / adaptive
#'   permutation procedure of Besag & Clifford (1991) and Gandy (2009)
#'   via [adaptive_permutation_pvalue()]. Default is `FALSE` (fixed-`B`
#'   permutation, byte-identical to previous behaviour).
#' @param B_max Integer. Hard upper bound on permutations when
#'   `adaptive = TRUE`. Default is 9999.
#' @param batch_size Integer. Permutations drawn per adaptive batch.
#'   Default is 100.
#' @param seed Integer or `NULL`. Random seed for reproducibility.
#'
#' @family tests-base
#' @seealso [hsic_test()] for independence testing; [bd_hsic_test()],
#'   [dr_date_test()] for causal versions.
#'
#' @return An object of class `"kernel_test_result"` with components:
#'   \describe{
#'     \item{statistic}{The observed MMD^2 statistic (unbiased).}
#'     \item{p_value}{Permutation p-value.}
#'     \item{method}{`"MMD"`.}
#'     \item{n}{Total sample size (n_x + n_y).}
#'     \item{n_permutations}{Number of permutations used.}
#'     \item{null_distribution}{Vector of permuted MMD^2 values.}
#'     \item{kernel_x}{Kernel specification used.}
#'     \item{call}{The matched call.}
#'     \item{n_perms_used}{(Adaptive only) Number of permutations
#'       actually drawn.}
#'     \item{stop_reason}{(Adaptive only) Reason the sequential
#'       procedure terminated.}
#'   }
#'
#' @references
#' Gretton, A., Borgwardt, K. M., Rasch, M. J., Scholkopf, B., & Smola,
#' A. (2012). A kernel two-sample test. *JMLR*, 13, 723-773.
#'
#' @examples
#' set.seed(42)
#'
#' # Same distribution
#' x <- matrix(rnorm(200), 100, 2)
#' y <- matrix(rnorm(200), 100, 2)
#' result <- mmd_test(x, y)
#' print(result)
#'
#' # Different distributions
#' y_shifted <- matrix(rnorm(200, mean = 1), 100, 2)
#' result <- mmd_test(x, y_shifted)
#' print(result)
#'
#' @export
mmd_test <- function(x, y,
                     kernel = kernel_spec(),
                     n_permutations = 500L,
                     alpha = 0.05,
                     adaptive = FALSE,
                     B_max = 9999L,
                     batch_size = 100L,
                     seed = NULL) {
  cl <- match.call()

  x <- as.matrix(x)
  y <- as.matrix(y)
  nx <- nrow(x)
  ny <- nrow(y)

  if (ncol(x) != ncol(y)) {
    stop("`x` and `y` must have the same number of columns.", call. = FALSE)
  }
  if (nx < 5 || ny < 5) {
    stop("Each sample must have at least 5 observations.", call. = FALSE)
  }
  n_permutations <- as.integer(n_permutations)

  if (!is.null(seed)) set.seed(seed)

  # Resolve bandwidth on pooled data
  pooled <- rbind(x, y)
  kernel <- resolve_bandwidth(kernel, pooled)

  # Compute kernel matrices
  K_pool <- kernel_matrix(pooled, kernel = kernel)
  Kxx <- K_pool[1:nx, 1:nx, drop = FALSE]
  Kyy <- K_pool[(nx + 1):(nx + ny), (nx + 1):(nx + ny), drop = FALSE]
  Kxy <- K_pool[1:nx, (nx + 1):(nx + ny), drop = FALSE]

  # Observed statistic
  stat_obs <- mmd2_unbiased_cpp(Kxx, Kyy, Kxy)

  if (isTRUE(adaptive)) {
    perm_fun <- function(n_perms) {
      as.numeric(permutation_mmd_cpp(K_pool, nx, ny, as.integer(n_perms)))
    }
    adp <- adaptive_permutation_pvalue(
      observed = stat_obs,
      perm_fun = perm_fun,
      B_min = max(99L, as.integer(batch_size)),
      B_max = as.integer(B_max),
      batch_size = as.integer(batch_size),
      alpha = alpha,
      seed = NULL
    )
    return(structure(
      list(
        statistic = stat_obs,
        p_value = adp$p_value,
        method = "MMD",
        n = nx + ny,
        n_permutations = adp$n_perms_used,
        null_distribution = adp$null_distribution,
        ess = NA_real_,
        weights = NULL,
        kernel_x = kernel,
        kernel_y = NULL,
        n_perms_used = adp$n_perms_used,
        stop_reason = adp$stop_reason,
        call = cl
      ),
      class = "kernel_test_result"
    ))
  }

  # Permutation null distribution
  null_dist <- permutation_mmd_cpp(K_pool, nx, ny, n_permutations)

  # p-value
  p_value <- (1 + sum(null_dist >= stat_obs)) / (1 + n_permutations)

  structure(
    list(
      statistic = stat_obs,
      p_value = p_value,
      method = "MMD",
      n = nx + ny,
      n_permutations = n_permutations,
      null_distribution = as.numeric(null_dist),
      ess = NA_real_,
      weights = NULL,
      kernel_x = kernel,
      kernel_y = NULL,
      call = cl
    ),
    class = "kernel_test_result"
  )
}
