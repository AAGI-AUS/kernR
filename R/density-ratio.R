#' Estimate Density Ratios via Logistic NCE
#'
#' Estimates the density ratio p*(x) / p(x|z) using Noise Contrastive
#' Estimation. Trains a classifier to distinguish joint samples (x, z)
#' from product-of-marginals samples (x_perm, z), where x_perm is drawn
#' from the reference distribution p*(x).
#'
#' @param x Numeric vector or matrix. Treatment variable.
#' @param z Numeric matrix or data.frame. Confounders.
#' @param method Character. Classification method: `"logistic"` (default),
#'   `"ranger"`, or `"xgboost"`.
#' @param n_noise Integer. Number of noise samples per real sample.
#'   Default is 1.
#' @param seed Integer or `NULL`. Random seed.
#'
#' @return A list of class `"density_ratio_fit"` with components:
#'   \describe{
#'     \item{weights}{Estimated density ratios for each observation.}
#'     \item{ess}{Effective sample size.}
#'     \item{method}{Method used.}
#'   }
#'
#' @family nuisance
#' @seealso [estimate_propensity()], [bd_hsic_test()].
#'
#' @examples
#' set.seed(42)
#' n <- 200
#' z <- matrix(rnorm(n * 2), n, 2)
#' x <- z[, 1] + rnorm(n)
#' dr <- estimate_density_ratio(x, z)
#' dr$ess
#'
#' @export
estimate_density_ratio <- function(x, z,
                                   method = c("logistic", "ranger", "xgboost"),
                                   n_noise = 1L,
                                   seed = NULL) {
  method <- match.arg(method)
  x <- as.matrix(x)
  z <- as.matrix(z)
  n <- nrow(x)

  if (nrow(z) != n) {
    stop("`x` and `z` must have the same number of rows.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)

  # Create joint samples (label = 1) and product-of-marginals (label = 0)
  # Product-of-marginals: permute x independently of z
  n_neg <- n * n_noise
  x_perm <- x[sample.int(n, n_neg, replace = TRUE), , drop = FALSE]
  z_rep <- z[rep(seq_len(n), n_noise), , drop = FALSE]

  features_pos <- cbind(x, z)
  features_neg <- cbind(x_perm, z_rep)
  features <- rbind(features_pos, features_neg)
  labels <- c(rep(1L, n), rep(0L, n_neg))

  # Train classifier
  probs <- switch(method,
    logistic = {
      df <- data.frame(y = labels, features)
      fit <- glm(y ~ ., data = df, family = binomial())
      predict(fit, type = "response")
    },
    ranger = {
      if (!requireNamespace("ranger", quietly = TRUE)) {
        stop("Package 'ranger' required for method = 'ranger'.", call. = FALSE)
      }
      df <- data.frame(y = factor(labels), features)
      fit <- ranger::ranger(y ~ ., data = df, probability = TRUE, num.trees = 500)
      fit$predictions[, "1"]
    },
    xgboost = {
      if (!requireNamespace("xgboost", quietly = TRUE)) {
        stop("Package 'xgboost' required for method = 'xgboost'.", call. = FALSE)
      }
      dtrain <- xgboost::xgb.DMatrix(data = features, label = labels)
      fit <- xgboost::xgb.train(
        params = list(objective = "binary:logistic", max_depth = 4, eta = 0.1),
        data = dtrain,
        nrounds = 100,
        verbose = 0
      )
      predict(fit, features)
    }
  )

  # Extract probabilities for the real (joint) samples
  probs_joint <- probs[1:n]

  # Density ratio: r(x, z) = p(joint) / p(marginal) ~ p / (1-p) * (n_noise/1)
  # Clamp probabilities to avoid division by zero
  probs_joint <- pmax(pmin(probs_joint, 1 - 1e-6), 1e-6)
  ratios <- (probs_joint / (1 - probs_joint)) * n_noise

  # Normalise to get importance weights
  weights <- ratios / sum(ratios) * n

  ess <- effective_sample_size(weights)

  structure(
    list(
      weights = as.numeric(weights),
      ratios = as.numeric(ratios),
      ess = ess,
      method = method,
      n = n
    ),
    class = "density_ratio_fit"
  )
}

#' Estimate Density Ratios via RuLSIF
#'
#' Relative unconstrained Least-Squares Importance Fitting.
#' Kernel-based closed-form density ratio estimation.
#'
#' @param x_num Numeric matrix. Numerator samples.
#' @param x_den Numeric matrix. Denominator samples.
#' @param kernel Kernel specification. Default is RBF with median heuristic.
#' @param lambda Regularisation parameter. Default is 0.1.
#' @param alpha Relative parameter (0 = LSIF, 0.5 = symmetric). Default is 0.
#'
#' @return Named list with `weights` and `ess`.
#'
#' @keywords internal
estimate_rulsif <- function(x_num, x_den,
                            kernel = kernel_spec(),
                            lambda = 0.1,
                            alpha = 0) {
  x_num <- as.matrix(x_num)
  x_den <- as.matrix(x_den)

  n_num <- nrow(x_num)
  n_den <- nrow(x_den)

  kernel <- resolve_bandwidth(kernel, rbind(x_num, x_den))

  # Kernel matrices
  K_num_num <- kernel_matrix(x_num, kernel = kernel)
  K_den_num <- kernel_matrix(x_den, x_num, kernel = kernel)

  # H = (1-alpha)/n_den * K_den_num^T K_den_num + alpha/n_num * K_num_num
  H <- (1 - alpha) / n_den * crossprod(K_den_num) +
    alpha / n_num * K_num_num +
    lambda * diag(n_num)

  # h = colMeans(K_num_num)
  h <- colMeans(K_num_num)

  # Solve H theta = h
  theta <- solve(H, h)
  theta <- pmax(theta, 0)  # Non-negativity

  # Compute weights for denominator samples
  weights <- as.numeric(K_den_num %*% theta)
  weights <- pmax(weights, 1e-8)
  weights <- weights / sum(weights) * n_den

  ess <- effective_sample_size(weights)

  list(weights = weights, ess = ess, method = "rulsif")
}

#' @export
print.density_ratio_fit <- function(x, ...) {
  cat("Density ratio estimation (", x$method, ")\n")
  cat("  N:   ", x$n, "\n")
  cat("  ESS: ", formatC(x$ess, digits = 1, format = "f"), "\n")
  cat("  Weight range: [",
    formatC(min(x$weights), digits = 3, format = "g"), ", ",
    formatC(max(x$weights), digits = 3, format = "g"), "]\n",
    sep = ""
  )
  invisible(x)
}
