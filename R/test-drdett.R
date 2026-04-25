#' Doubly Robust Distributional Effect on the Treated Test (DR-DETT)
#'
#' Tests whether the distribution of the treated potential outcome Y(1)
#' differs from the control potential outcome Y(0) among the treated
#' subpopulation. Requires only one-sided overlap.
#'
#' @inheritParams dr_date_test
#'
#' @return An object of class `"kernel_test_result"`.
#'
#' @family tests-causal
#' @seealso [dr_date_test()] (average treatment effect distribution),
#'   [kernel_causal_test()] (formula interface).
#'
#' @details
#' DR-DETT is analogous to DR-DATE but focuses on the **effect on the
#' treated** (ETT) rather than the average treatment effect (ATE). It
#' requires only one-sided overlap: P(T=1|X) > epsilon (not bounded away
#' from 1), making it applicable in settings where positivity is
#' violated for the control group.
#'
#' @references
#' Fawkes, J., Hu, R., Evans, R. J., & Sejdinovic, D. (2024). Doubly
#' robust kernel statistics for testing distributional treatment effects.
#' *Transactions on Machine Learning Research*.
#'
#' @examples
#' set.seed(42)
#' n <- 300
#' x <- matrix(rnorm(n * 2), n, 2)
#' t <- rbinom(n, 1, plogis(0.5 * x[, 1]))
#' y <- t * rnorm(n, sd = 2) + (1 - t) * rnorm(n, sd = 1) + x[, 1]
#'
#' result <- dr_dett_test(y, t, x, n_permutations = 200, seed = 1)
#' print(result)
#'
#' @export
dr_dett_test <- function(y, treatment, covariates,
                         kernel_y = kernel_spec(),
                         propensity_model = c("logistic", "ranger", "xgboost"),
                         outcome_model = c("krr", "zero"),
                         n_permutations = 500L,
                         n_bins = 10L,
                         regularisation = "cv",
                         alpha = 0.05,
                         adaptive = FALSE,
                         B_max = 9999L,
                         batch_size = 100L,
                         seed = NULL,
                         verbose = FALSE) {
  cl <- match.call()
  propensity_model <- match.arg(propensity_model)
  outcome_model <- match.arg(outcome_model)
  n_permutations <- as.integer(n_permutations)

  y <- validate_input(y, "y", min_n = 30)
  covariates <- validate_input(covariates, "covariates", min_n = 30)
  treatment <- as.integer(treatment)
  n <- nrow(y)

  if (length(treatment) != n || nrow(covariates) != n) {
    stop("All inputs must have the same number of observations.", call. = FALSE)
  }
  if (!all(treatment %in% c(0L, 1L))) {
    stop("`treatment` must be binary (0/1).", call. = FALSE)
  }

  if (!is.null(seed)) set.seed(seed)

  # Estimate propensity scores
  ps <- estimate_propensity(treatment, covariates, method = propensity_model)
  e_hat <- ps$scores

  kernel_y <- resolve_bandwidth(kernel_y, y)
  Ky <- kernel_matrix(y, kernel = kernel_y)

  # Compute DR-DETT statistic
  stat_obs <- compute_dr_dett_stat(Ky, treatment, e_hat)

  if (isTRUE(adaptive)) {
    perm_fun <- function(n_perms) {
      out <- numeric(n_perms)
      for (p in seq_len(n_perms)) {
        t_perm <- bin_permute_treatment(treatment, e_hat, n_bins)
        out[p] <- compute_dr_dett_stat(Ky, t_perm, e_hat)
      }
      out
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
        method = "DR-DETT",
        n = n,
        n_permutations = adp$n_perms_used,
        null_distribution = adp$null_distribution,
        ess = NA_real_,
        weights = e_hat,
        kernel_x = NULL,
        kernel_y = kernel_y,
        n_perms_used = adp$n_perms_used,
        stop_reason = adp$stop_reason,
        call = cl
      ),
      class = "kernel_test_result"
    ))
  }

  # Permutation null
  null_dist <- numeric(n_permutations)
  for (p in seq_len(n_permutations)) {
    t_perm <- bin_permute_treatment(treatment, e_hat, n_bins)
    null_dist[p] <- compute_dr_dett_stat(Ky, t_perm, e_hat)
  }

  p_value <- (1 + sum(null_dist >= stat_obs)) / (1 + n_permutations)

  structure(
    list(
      statistic = stat_obs,
      p_value = p_value,
      method = "DR-DETT",
      n = n,
      n_permutations = n_permutations,
      null_distribution = null_dist,
      ess = NA_real_,
      weights = e_hat,
      kernel_x = NULL,
      kernel_y = kernel_y,
      call = cl
    ),
    class = "kernel_test_result"
  )
}

#' Compute DR-DETT Test Statistic
#'
#' @param Ky n x n outcome kernel matrix.
#' @param treatment Binary treatment vector.
#' @param e_hat Propensity scores.
#'
#' @return Scalar test statistic.
#' @keywords internal
compute_dr_dett_stat <- function(Ky, treatment, e_hat) {
  n <- length(treatment)
  t <- treatment
  n1 <- sum(t)

  if (n1 < 2) return(0)

  # DETT uses propensity odds weighting for the control counterfactual
  # w(x) = (1 - e(x)) / e(x) for control units
  # Treated: w_1 = t / n1 (uniform over treated)
  # Control counterfactual: w_0 = (1-t) * (1-e)/(e * sum((1-t)*(1-e)/e))

  w1 <- t / n1

  odds <- (1 - e_hat) / e_hat
  w0 <- (1 - t) * odds
  sw0 <- sum(w0)
  if (sw0 > 0) w0 <- w0 / sw0 else w0 <- rep(0, n)

  term_11 <- as.numeric(crossprod(w1, Ky %*% w1))
  term_00 <- as.numeric(crossprod(w0, Ky %*% w0))
  term_10 <- as.numeric(crossprod(w1, Ky %*% w0))

  mmd2 <- term_11 + term_00 - 2 * term_10
  max(mmd2, 0)
}
