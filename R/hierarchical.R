#' Hierarchical Kernel Causal Test
#'
#' Extends bd-HSIC and DR-DATE/DR-DETT to hierarchical (nested/clustered)
#' data by decomposing the test statistic into within-cluster and
#' between-cluster components.
#'
#' @param y Numeric vector or matrix. Outcome.
#' @param treatment Treatment variable (binary for DR tests, any for bd-HSIC).
#' @param covariates Numeric matrix of confounders.
#' @param cluster_id Factor or integer vector identifying clusters.
#' @param method Character. `"dr-date"` (default), `"dr-dett"`, or `"bd-hsic"`.
#' @param kernel_y Kernel specification for outcomes.
#' @param n_permutations Integer. Number of permutations. Default is 500.
#' @param weight_method Character. How to weight within/between components:
#'   `"equal"` (default), `"icc"` (variance decomposition), or `"within_only"`.
#' @param seed Integer or `NULL`.
#' @param verbose Logical.
#' @param ... Additional arguments passed to the underlying test.
#'
#' @family tests-hierarchical
#' @seealso [bd_hsic_test()], [dr_date_test()], [dr_dett_test()] for the
#'   non-hierarchical analogues.
#'
#' @return An object of class `"kernel_test_result"` with additional
#'   `hierarchical` component containing within/between statistics.
#'
#' @details
#' For clustered data (e.g., patients within hospitals, plots within
#' farms), standard kernel tests may have inflated type I error because
#' observations within the same cluster are not independent.
#'
#' This function decomposes the test into:
#' - **Within-cluster**: Average of within-cluster test statistics
#'   (tests for treatment effects within each cluster).
#' - **Between-cluster**: Test on cluster-level mean embeddings
#'   (tests for treatment effects across clusters).
#'
#' The combined statistic is a weighted sum, with weights determined
#' by `weight_method`. Permutation is performed within clusters to
#' preserve the hierarchical structure.
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' n_clusters <- 20
#' n_per <- 30
#' n <- n_clusters * n_per
#' cluster_id <- rep(1:n_clusters, each = n_per)
#'
#' # Cluster-level random effects
#' cluster_effect <- rnorm(n_clusters, sd = 1)[cluster_id]
#' x <- matrix(rnorm(n * 2), n, 2)
#' t <- rbinom(n, 1, plogis(0.3 * x[, 1]))
#' y <- 0.5 * t + cluster_effect + x[, 1] + rnorm(n)
#'
#' result <- hierarchical_test(y, t, x, cluster_id,
#'   method = "dr-date",
#'   n_permutations = 100,
#'   seed = 1
#' )
#' print(result)
#' }
#'
#' @export
hierarchical_test <- function(y, treatment, covariates, cluster_id, # nolint: cyclocomp_linter.
                              method = c("dr-date", "dr-dett", "bd-hsic"),
                              kernel_y = kernel_spec(),
                              n_permutations = 500L,
                              weight_method = c("equal", "icc", "within_only"),
                              seed = NULL,
                              verbose = FALSE,
                              ...) {
  cl <- match.call()
  method <- match.arg(method)
  weight_method <- match.arg(weight_method)

  y <- validate_input(y, "y")
  covariates <- validate_input(covariates, "covariates")
  n <- nrow(y)
  cluster_id <- as.factor(cluster_id)
  clusters <- levels(cluster_id)
  n_clust <- length(clusters)

  if (n_clust < 3) {
    stop("At least 3 clusters are required.", call. = FALSE)
  }

  if (!is.null(seed)) set.seed(seed)

  kernel_y <- resolve_bandwidth(kernel_y, y)

  # Within-cluster statistics
  within_stats <- numeric(n_clust)
  within_n <- integer(n_clust)

  for (j in seq_along(clusters)) {
    idx <- which(cluster_id == clusters[j])
    nj <- length(idx)
    within_n[j] <- nj

    if (nj < 10) next  # Skip tiny clusters

    tryCatch(
      {
        if (method == "bd-hsic") {
          res_j <- bd_hsic_test(
            x = treatment[idx],
            y = y[idx, , drop = FALSE],
            z = covariates[idx, , drop = FALSE],
            kernel_y = kernel_y,
            n_permutations = 0L,
            ...
          )
        } else {
          test_fn <- if (method == "dr-date") dr_date_test else dr_dett_test
          res_j <- test_fn(
            y = y[idx, , drop = FALSE],
            treatment = treatment[idx],
            covariates = covariates[idx, , drop = FALSE],
            kernel_y = kernel_y,
            n_permutations = 0L,
            ...
          )
        }
        within_stats[j] <- res_j$statistic
      },
      error = function(e) {
        within_stats[j] <<- NA_real_
      }
    )
  }

  # Between-cluster statistic: cluster-level mean embeddings
  Ky_full <- kernel_matrix(y, kernel = kernel_y)

  # Compute cluster-level mean embedding kernel
  K_between <- matrix(0, n_clust, n_clust)
  for (j1 in seq_along(clusters)) {
    idx1 <- which(cluster_id == clusters[j1])
    for (j2 in j1:n_clust) {
      idx2 <- which(cluster_id == clusters[j2])
      K_between[j1, j2] <- mean(Ky_full[idx1, idx2])
      K_between[j2, j1] <- K_between[j1, j2]
    }
  }

  # Between-cluster treatment: majority treatment per cluster
  cluster_treatment <- tapply(as.integer(treatment), cluster_id, function(tt) {
    as.integer(mean(tt) > 0.5)
  })

  # Between-cluster MMD
  idx_t1 <- which(cluster_treatment == 1)
  idx_t0 <- which(cluster_treatment == 0)

  if (length(idx_t1) >= 2 && length(idx_t0) >= 2) {
    Kxx <- K_between[idx_t1, idx_t1, drop = FALSE]
    Kyy <- K_between[idx_t0, idx_t0, drop = FALSE]
    Kxy <- K_between[idx_t1, idx_t0, drop = FALSE]
    between_stat <- mmd2_unbiased_cpp(Kxx, Kyy, Kxy)
  } else {
    between_stat <- 0
  }

  # Combine within and between
  valid_within <- !is.na(within_stats)
  avg_within <- if (sum(valid_within) > 0) {
    weighted.mean(within_stats[valid_within], within_n[valid_within])
  } else {
    0
  }

  combined_stat <- switch(weight_method,
    equal = 0.5 * avg_within + 0.5 * between_stat,
    icc = {
      # ICC-like: fraction of variance between clusters
      var_between <- var(tapply(y[, 1], cluster_id, mean))
      var_total <- var(y[, 1])
      icc <- max(0, min(1, var_between / (var_total + 1e-10)))
      icc * between_stat + (1 - icc) * avg_within
    },
    within_only = avg_within
  )

  # Permutation: permute treatment within clusters
  null_dist <- numeric(n_permutations)
  for (p in seq_len(n_permutations)) {
    t_perm <- treatment
    for (j in seq_along(clusters)) {
      idx <- which(cluster_id == clusters[j])
      if (length(idx) > 1) {
        t_perm[idx] <- treatment[idx[sample.int(length(idx))]]
      }
    }

    # Recompute within stats
    perm_within <- numeric(n_clust)
    for (j in seq_along(clusters)) {
      idx <- which(cluster_id == clusters[j])
      nj <- length(idx)
      if (nj < 10) next

      tryCatch(
        {
          if (method == "bd-hsic") {
            res_j <- bd_hsic_test(
              x = t_perm[idx],
              y = y[idx, , drop = FALSE],
              z = covariates[idx, , drop = FALSE],
              kernel_y = kernel_y,
              n_permutations = 0L,
              ...
            )
          } else {
            test_fn <- if (method == "dr-date") dr_date_test else dr_dett_test
            res_j <- test_fn(
              y = y[idx, , drop = FALSE],
              treatment = t_perm[idx],
              covariates = covariates[idx, , drop = FALSE],
              kernel_y = kernel_y,
              n_permutations = 0L,
              ...
            )
          }
          perm_within[j] <- res_j$statistic
        },
        error = function(e) {
          perm_within[j] <<- NA_real_
        }
      )
    }

    perm_avg <- if (sum(!is.na(perm_within)) > 0) {
      weighted.mean(perm_within[!is.na(perm_within)], within_n[!is.na(perm_within)])
    } else {
      0
    }

    # Between
    perm_ct <- tapply(as.integer(t_perm), cluster_id, function(tt) {
      as.integer(mean(tt) > 0.5)
    })
    p_t1 <- which(perm_ct == 1)
    p_t0 <- which(perm_ct == 0)
    if (length(p_t1) >= 2 && length(p_t0) >= 2) {
      p_between <- mmd2_unbiased_cpp(
        K_between[p_t1, p_t1, drop = FALSE],
        K_between[p_t0, p_t0, drop = FALSE],
        K_between[p_t1, p_t0, drop = FALSE]
      )
    } else {
      p_between <- 0
    }

    null_dist[p] <- switch(weight_method,
      equal = 0.5 * perm_avg + 0.5 * p_between,
      icc = {
        icc <- max(0, min(1, var_between / (var_total + 1e-10)))
        icc * p_between + (1 - icc) * perm_avg
      },
      within_only = perm_avg
    )
  }

  p_value <- (1 + sum(null_dist >= combined_stat)) / (1 + n_permutations)

  structure(
    list(
      statistic = combined_stat,
      p_value = p_value,
      method = paste0("Hierarchical-", toupper(sub("-", "", method))),
      n = n,
      n_permutations = n_permutations,
      null_distribution = null_dist,
      ess = NA_real_,
      weights = NULL,
      kernel_x = NULL,
      kernel_y = kernel_y,
      call = cl,
      hierarchical = list(
        n_clusters = n_clust,
        within_stats = within_stats,
        between_stat = between_stat,
        weight_method = weight_method
      )
    ),
    class = "kernel_test_result"
  )
}
