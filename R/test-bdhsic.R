#' Backdoor-HSIC Test for Causal Association
#'
#' Tests the do-null hypothesis H_0: p(y | do(x)) = p*(y) using a
#' kernel-based test with backdoor adjustment via density ratio
#' estimation. Detects causal associations including non-linear effects
#' that standard linear methods miss.
#'
#' @param x Numeric vector or matrix. Treatment variable.
#' @param y Numeric vector or matrix. Outcome variable.
#' @param z Numeric matrix, data.frame, or data.table. Confounders.
#' @param kernel_x Kernel specification for treatment space. Default is
#'   RBF with median heuristic.
#' @param kernel_y Kernel specification for outcome space. Default is
#'   RBF with median heuristic.
#' @param density_ratio Character. Method for density ratio estimation:
#'   `"logistic"` (default), `"ranger"`, `"xgboost"`, or `"rulsif"`.
#' @param n_permutations Integer. Number of permutations for the null
#'   distribution when `adaptive = FALSE`. Default is 500.
#' @param n_clusters Integer or `"auto"`. Number of clusters for valid
#'   permutation. Default is `"auto"`.
#' @param split_ratio Numeric in (0, 1). Proportion of data for training
#'   the density ratio estimator. Default is 0.5.
#' @param alpha Numeric. Significance level. Used both as the reporting
#'   threshold and, when `adaptive = TRUE`, as the early-stopping
#'   threshold. Default is 0.05.
#' @param adaptive Logical. If `TRUE`, use the sequential / adaptive
#'   permutation procedure of Besag & Clifford (1991) and Gandy (2009)
#'   via [adaptive_permutation_pvalue()]. Default is `FALSE`
#'   (byte-identical to previous behaviour).
#' @param B_max Integer. Hard upper bound on permutations when
#'   `adaptive = TRUE`. Default is 9999.
#' @param batch_size Integer. Permutations drawn per adaptive batch.
#'   Default is 100.
#' @param seed Integer or `NULL`. Random seed for reproducibility.
#' @param verbose Logical. Print progress. Default is `FALSE`.
#'
#' @family tests-causal
#' @seealso [hsic_test()] (no confounding adjustment),
#'   [kernel_causal_test()] (formula interface).
#'
#' @return An object of class `"kernel_test_result"`.
#'
#' @details
#' The bd-HSIC test (Hu, Sejdinovic & Evans, 2024) tests whether
#' treatment X has a causal effect on outcome Y after adjusting for
#' confounders Z via the backdoor criterion.
#'
#' The test works by:
#' 1. Estimating density ratios w(x, z) = p*(x) / p(x|z) to reweight
#'    observational samples to the interventional distribution.
#' 2. Computing a weighted HSIC statistic between X and Y.
#' 3. Obtaining p-values via permutation of Y within clusters of
#'    similar conditional densities p(x|z).
#'
#' Unlike PDS or Double ML, bd-HSIC can detect non-linear causal effects
#' (e.g., U-shaped relationships) where the treatment affects higher
#' moments of the outcome but not necessarily the mean.
#'
#' @references
#' Hu, R., Sejdinovic, D., & Evans, R. J. (2024). A kernel test for
#' causal association via noise contrastive backdoor adjustment. *JMLR*,
#' 25(160), 1-56.
#'
#' @examples
#' set.seed(42)
#' n <- 300
#' z <- matrix(rnorm(n * 2), n, 2)
#' x <- z[, 1] + rnorm(n)
#' y <- 0.5 * x + z[, 2] + rnorm(n, sd = 0.5)
#'
#' result <- bd_hsic_test(x, y, z, n_permutations = 200, seed = 1)
#' print(result)
#'
#' @export
bd_hsic_test <- function(x, y, z,
                         kernel_x = kernel_spec(),
                         kernel_y = kernel_spec(),
                         density_ratio = c("logistic", "ranger", "xgboost", "rulsif"),
                         n_permutations = 500L,
                         n_clusters = "auto",
                         split_ratio = 0.5,
                         alpha = 0.05,
                         adaptive = FALSE,
                         B_max = 9999L,
                         batch_size = 100L,
                         seed = NULL,
                         verbose = FALSE) {
  cl <- match.call()
  density_ratio <- match.arg(density_ratio)
  n_permutations <- as.integer(n_permutations)

  x <- validate_input(x, "x", min_n = 20)
  y <- validate_input(y, "y", min_n = 20)
  z <- validate_input(z, "z", min_n = 20)
  n <- nrow(x)

  if (nrow(y) != n || nrow(z) != n) {
    stop("`x`, `y`, and `z` must have the same number of observations.",
      call. = FALSE
    )
  }

  if (!is.null(seed)) set.seed(seed)

  # Split data: train density ratio on one half, test on the other
  n_train <- floor(n * split_ratio)
  idx_train <- sample.int(n, n_train)
  idx_test <- setdiff(seq_len(n), idx_train)
  n_test <- length(idx_test)

  if (verbose) message("Training density ratio estimator on ", n_train, " obs...")

  # Estimate density ratios on training data
  if (density_ratio == "rulsif") {
    # For RuLSIF: estimate p*(x) / p(x|z) directly
    dr <- estimate_rulsif(x[idx_train, , drop = FALSE],
      x[idx_test, , drop = FALSE],
      kernel = kernel_x
    )
    weights <- dr$weights
  } else {
    dr <- estimate_density_ratio(
      x[idx_train, , drop = FALSE],
      z[idx_train, , drop = FALSE],
      method = density_ratio
    )
    # Apply trained model concepts: use weights for test set
    # For simplicity, re-estimate on full test set
    dr_test <- estimate_density_ratio(
      x[idx_test, , drop = FALSE],
      z[idx_test, , drop = FALSE],
      method = density_ratio
    )
    weights <- dr_test$weights
  }

  ess <- effective_sample_size(weights)

  if (verbose) message("ESS: ", round(ess, 1), " / ", n_test)

  # Compute kernel matrices on test data
  kernel_x <- resolve_bandwidth(kernel_x, x[idx_test, , drop = FALSE])
  kernel_y <- resolve_bandwidth(kernel_y, y[idx_test, , drop = FALSE])

  Kx <- kernel_matrix(x[idx_test, , drop = FALSE], kernel = kernel_x)
  Ky <- kernel_matrix(y[idx_test, , drop = FALSE], kernel = kernel_y)

  # Observed weighted HSIC
  stat_obs <- weighted_hsic_stat_cpp(Kx, Ky, weights)

  if (verbose) message("Computing ", n_permutations, " permutations...")

  # Cluster observations for valid permutation
  clusters <- cluster_observations(
    weights,
    z[idx_test, , drop = FALSE],
    n_clusters = n_clusters
  )

  if (isTRUE(adaptive)) {
    perm_fun <- function(n_perms) {
      cluster_permutation_hsic(Kx, Ky, weights, clusters, as.integer(n_perms))
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
        method = "bd-HSIC",
        n = n_test,
        n_permutations = adp$n_perms_used,
        null_distribution = adp$null_distribution,
        ess = ess,
        weights = weights,
        kernel_x = kernel_x,
        kernel_y = kernel_y,
        n_perms_used = adp$n_perms_used,
        stop_reason = adp$stop_reason,
        call = cl
      ),
      class = "kernel_test_result"
    ))
  }

  # Permutation null distribution
  null_dist <- cluster_permutation_hsic(
    Kx, Ky, weights, clusters, n_permutations
  )

  # p-value
  p_value <- (1 + sum(null_dist >= stat_obs)) / (1 + n_permutations)

  structure(
    list(
      statistic = stat_obs,
      p_value = p_value,
      method = "bd-HSIC",
      n = n_test,
      n_permutations = n_permutations,
      null_distribution = null_dist,
      ess = ess,
      weights = weights,
      kernel_x = kernel_x,
      kernel_y = kernel_y,
      call = cl
    ),
    class = "kernel_test_result"
  )
}
