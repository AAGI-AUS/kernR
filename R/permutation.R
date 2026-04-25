#' Cluster-Based Permutation for bd-HSIC
#'
#' Permutes Y indices within clusters of similar conditional densities
#' p(x|z), ensuring valid exchangeability under the null.
#'
#' @param Kx n x n kernel matrix for X.
#' @param Ky n x n kernel matrix for Y.
#' @param weights Density ratio weights.
#' @param clusters Integer vector of cluster assignments.
#' @param n_permutations Number of permutations.
#'
#' @return Vector of permuted weighted HSIC statistics.
#' @keywords internal
cluster_permutation_hsic <- function(Kx, Ky, weights, clusters,
                                     n_permutations) {
  n <- nrow(Kx)
  n_clusters <- max(clusters)
  results <- numeric(n_permutations)

  for (p in seq_len(n_permutations)) {
    # Permute Y indices within each cluster
    perm <- seq_len(n)
    for (k in seq_len(n_clusters)) {
      idx <- which(clusters == k)
      if (length(idx) > 1) {
        perm[idx] <- idx[sample.int(length(idx))]
      }
    }
    Ky_perm <- Ky[perm, perm]
    results[p] <- weighted_hsic_stat_cpp(Kx, Ky_perm, weights)
  }

  results
}

#' Bin-Based Permutation for DR Tests
#'
#' Permutes treatment labels within propensity score bins.
#'
#' @param treatment Binary treatment vector.
#' @param propensity_scores Propensity score vector.
#' @param n_bins Number of bins.
#'
#' @return Permuted treatment vector.
#' @keywords internal
bin_permute_treatment <- function(treatment, propensity_scores, n_bins = 10L) {
  bins <- as.integer(cut(propensity_scores,
    breaks = n_bins,
    labels = FALSE,
    include.lowest = TRUE
  ))

  perm_idx <- stratified_permute_cpp(bins, n_bins)
  # C++ returns 0-indexed
  treatment[perm_idx + 1L]
}

#' Adaptive Sequential Permutation P-Value
#'
#' Implements the sequential / adaptive permutation procedure of Besag &
#' Clifford (1991) and Gandy (2009): instead of running a fixed number of
#' permutations, draw permuted statistics in batches and stop early when
#' the resampled p-value's Wilson confidence interval clearly excludes
#' the significance level `alpha`. Type-I error is controlled by
#' construction.
#'
#' @param observed Numeric scalar. The observed test statistic.
#' @param perm_fun A function `function(n_perms)` returning a numeric
#'   vector of `n_perms` permuted statistics drawn under the null. The
#'   helper is test-agnostic: `perm_fun` encapsulates the sampling.
#' @param B_min Integer. Minimum number of permutations before stopping
#'   is considered. Default 99.
#' @param B_max Integer. Hard upper bound on the number of permutations.
#'   Default 9999.
#' @param batch_size Integer. Permutations drawn per batch. Default 100.
#' @param alpha Numeric in (0, 1). Significance level used to decide
#'   early stopping. Default 0.05.
#' @param wilson_level Numeric in (0, 1). Confidence level for the Wilson
#'   interval used to bound `p_hat`. Default 0.95.
#' @param seed Optional integer seed. If supplied, `set.seed(seed)` is
#'   called before any permutation is drawn, providing reproducibility
#'   across the entire adaptive procedure.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{p_value}{Final permutation p-value with the standard +1
#'       correction.}
#'     \item{n_perms_used}{Total permutations drawn.}
#'     \item{stop_reason}{One of `"upper_below_alpha"`,
#'       `"lower_above_alpha"`, or `"max_reached"`.}
#'     \item{null_distribution}{Numeric vector of length `n_perms_used`
#'       containing every permuted statistic drawn.}
#'   }
#'
#' @details
#' At each step the running p-value is computed as
#' `p_hat = (1 + n_extreme) / (1 + n_total_perms)` and a Wilson
#' confidence interval at level `wilson_level` is built around `p_hat`
#' using `n_total_perms` as the sample size. The procedure stops when
#' (a) the upper bound is below `alpha`, (b) the lower bound is above
#' `alpha`, or (c) `n_total_perms` reaches `B_max`. The first batch is
#' sized `max(batch_size, B_min)` so that no decision is made before
#' `B_min` permutations have been drawn.
#'
#' @references
#' Besag, J., & Clifford, P. (1991). Sequential Monte Carlo p-values.
#' *Biometrika*, 78(2), 301-304.
#'
#' Gandy, A. (2009). Sequential implementation of Monte Carlo tests with
#' uniformly bounded resampling risk. *Journal of the American
#' Statistical Association*, 104(488), 1504-1511.
#'
#' @keywords internal
adaptive_permutation_pvalue <- function(observed, perm_fun, # nolint: cyclocomp_linter.
                                        B_min = 99L, B_max = 9999L,
                                        batch_size = 100L,
                                        alpha = 0.05,
                                        wilson_level = 0.95,
                                        seed = NULL) {
  # Argument validation
  if (!is.function(perm_fun)) {
    stop("`perm_fun` must be a function of one argument (n_perms).",
      call. = FALSE
    )
  }
  if (!is.numeric(observed) || length(observed) != 1L || !is.finite(observed)) {
    stop("`observed` must be a finite numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L ||
        alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single number in (0, 1).", call. = FALSE)
  }
  if (!is.numeric(wilson_level) || length(wilson_level) != 1L ||
        wilson_level <= 0 || wilson_level >= 1) {
    stop("`wilson_level` must be a single number in (0, 1).", call. = FALSE)
  }
  B_min <- as.integer(B_min)
  B_max <- as.integer(B_max)
  batch_size <- as.integer(batch_size)
  if (length(B_min) != 1L || is.na(B_min) || B_min < 2L) {
    stop("`B_min` must be an integer >= 2.", call. = FALSE)
  }
  if (length(B_max) != 1L || is.na(B_max) || B_max < B_min) {
    stop("`B_max` must be an integer >= `B_min`.", call. = FALSE)
  }
  if (length(batch_size) != 1L || is.na(batch_size) || batch_size < 1L) {
    stop("`batch_size` must be a positive integer.", call. = FALSE)
  }

  if (!is.null(seed)) set.seed(seed)

  # Wilson interval critical value
  z <- stats::qnorm(1 - (1 - wilson_level) / 2)

  # Pre-allocate buffer up to B_max for efficiency
  null_dist <- numeric(B_max)
  n_total <- 0L
  n_extreme <- 0L
  stop_reason <- "max_reached"

  # First batch is at least B_min, then batch_size thereafter
  next_batch <- max(batch_size, B_min)

  repeat {
    # Cap by remaining budget
    take <- min(next_batch, B_max - n_total)
    if (take <= 0L) break

    new_perms <- perm_fun(take)
    if (!is.numeric(new_perms) || length(new_perms) != take) {
      stop("`perm_fun(n_perms)` must return a numeric vector of length ",
        "`n_perms`.",
        call. = FALSE
      )
    }
    new_perms <- as.numeric(new_perms)
    null_dist[(n_total + 1L):(n_total + take)] <- new_perms
    n_extreme <- n_extreme + sum(new_perms >= observed)
    n_total <- n_total + take

    p_hat <- (1 + n_extreme) / (1 + n_total)

    # Wilson interval around p_hat with sample size n_total
    denom <- 1 + z^2 / n_total
    centre <- (p_hat + z^2 / (2 * n_total)) / denom
    halfw <- (z * sqrt(p_hat * (1 - p_hat) / n_total +
                         z^2 / (4 * n_total^2))) / denom
    lower <- centre - halfw
    upper <- centre + halfw

    if (n_total >= B_min) {
      if (upper < alpha) {
        stop_reason <- "upper_below_alpha"
        break
      }
      if (lower > alpha) {
        stop_reason <- "lower_above_alpha"
        break
      }
    }

    if (n_total >= B_max) {
      stop_reason <- "max_reached"
      break
    }

    next_batch <- batch_size
  }

  null_dist <- null_dist[seq_len(n_total)]
  p_value <- (1 + n_extreme) / (1 + n_total)

  list(
    p_value = p_value,
    n_perms_used = n_total,
    stop_reason = stop_reason,
    null_distribution = null_dist
  )
}

#' Simple K-Means Clustering for Permutation Groups
#'
#' Clusters observations based on conditional density embeddings
#' using standard k-means on the density ratio weight space.
#'
#' @param weights Weight vector (density ratios or propensity scores).
#' @param z Confounder matrix.
#' @param n_clusters Number of clusters. If `"auto"`, selects by
#'   silhouette score (2 to 10 clusters).
#'
#' @return Integer vector of cluster assignments.
#' @keywords internal
cluster_observations <- function(weights, z, n_clusters = "auto") {
  z <- as.matrix(z)
  # Cluster on (weights, z) jointly
  features <- cbind(scale(weights), scale(z))

  if (identical(n_clusters, "auto")) {
    # Try 2 to min(10, n/5) clusters, pick by within-SS
    max_k <- min(10L, as.integer(nrow(features) / 5))
    max_k <- max(max_k, 2L)

    best_k <- 2L
    best_ratio <- Inf
    for (k in 2:max_k) {
      km <- stats::kmeans(features, centers = k, nstart = 5, iter.max = 50)
      ratio <- km$tot.withinss / km$totss
      if (ratio < best_ratio) {
        best_ratio <- ratio
        best_k <- k
      }
    }
    n_clusters <- best_k
  }

  km <- stats::kmeans(features, centers = as.integer(n_clusters),
    nstart = 10, iter.max = 100
  )
  km$cluster
}
