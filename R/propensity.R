#' Estimate Propensity Scores
#'
#' Estimates P(T = 1 | X) using the specified model, with built-in
#' cross-fitting support.
#'
#' @param treatment Binary vector (0/1). Treatment indicator.
#' @param covariates Numeric matrix or data.frame. Confounders.
#' @param method Character. `"logistic"` (default), `"ranger"`, or
#'   `"xgboost"`.
#' @param cross_fit Logical. If `TRUE`, uses 5-fold cross-fitting to
#'   produce out-of-sample propensity estimates. Default is `TRUE`.
#' @param n_folds Integer. Number of cross-fitting folds. Default is 5.
#' @param trim Numeric. Trim extreme propensity scores to `[trim, 1-trim]`.
#'   Default is 0.01.
#'
#' @return A list of class `"propensity_fit"` with components:
#'   \describe{
#'     \item{scores}{Estimated propensity scores P(T=1|X).}
#'     \item{method}{Method used.}
#'     \item{trim}{Trimming threshold applied.}
#'     \item{n_trimmed}{Number of scores that were trimmed.}
#'   }
#'
#' @family nuisance
#' @seealso [estimate_density_ratio()], [dr_date_test()],
#'   [dr_dett_test()], [assess_overlap()].
#'
#' @examples
#' set.seed(42)
#' n <- 300
#' x <- matrix(rnorm(n * 3), n, 3)
#' logit_p <- 0.5 * x[, 1] - 0.3 * x[, 2]
#' t <- rbinom(n, 1, plogis(logit_p))
#' ps <- estimate_propensity(t, x)
#' summary(ps$scores)
#'
#' @export
estimate_propensity <- function(treatment,
                                covariates,
                                method = c("logistic", "ranger", "xgboost"),
                                cross_fit = TRUE,
                                n_folds = 5L,
                                trim = 0.01) {
  method <- match.arg(method)
  treatment <- as.integer(treatment)
  covariates <- as.matrix(covariates)
  n <- length(treatment)

  if (nrow(covariates) != n) {
    stop("`treatment` and `covariates` must have the same number of observations.",
      call. = FALSE
    )
  }
  if (!all(treatment %in% c(0L, 1L))) {
    stop("`treatment` must be binary (0/1).", call. = FALSE)
  }

  if (cross_fit && n >= n_folds * 5) {
    scores <- cross_fit_propensity(treatment, covariates, method, n_folds)
  } else {
    scores <- fit_propensity_single(treatment, covariates, method)
  }

  # Trim
  n_trimmed <- sum(scores < trim | scores > (1 - trim))
  scores <- pmax(pmin(scores, 1 - trim), trim)

  structure(
    list(
      scores = scores,
      method = method,
      trim = trim,
      n_trimmed = n_trimmed,
      n = n
    ),
    class = "propensity_fit"
  )
}

#' Cross-Fit Propensity Scores
#'
#' @param treatment Binary vector.
#' @param covariates Numeric matrix.
#' @param method Classification method.
#' @param n_folds Number of folds.
#'
#' @return Vector of out-of-sample propensity scores.
#' @keywords internal
cross_fit_propensity <- function(treatment, covariates, method, n_folds) {
  n <- length(treatment)
  folds <- sample(rep(1:n_folds, length.out = n))
  scores <- numeric(n)

  for (k in 1:n_folds) {
    train_idx <- which(folds != k)
    test_idx <- which(folds == k)

    pred <- fit_propensity_single(
      treatment[train_idx],
      covariates[train_idx, , drop = FALSE],
      method,
      newdata = covariates[test_idx, , drop = FALSE]
    )
    scores[test_idx] <- pred
  }

  scores
}

#' Fit Propensity Model (Single Fold)
#'
#' @param treatment Binary vector.
#' @param covariates Numeric matrix.
#' @param method Classification method.
#' @param newdata Optional matrix for prediction. If NULL, predicts on
#'   training data.
#'
#' @return Vector of predicted probabilities.
#' @keywords internal
fit_propensity_single <- function(treatment, covariates, method,
                                  newdata = NULL) {
  if (is.null(newdata)) newdata <- covariates

  switch(method,
    logistic = {
      df_train <- data.frame(y = treatment, covariates)
      fit <- glm(y ~ ., data = df_train, family = binomial())
      df_new <- data.frame(newdata)
      names(df_new) <- names(df_train)[-1]
      as.numeric(predict(fit, newdata = df_new, type = "response"))
    },
    ranger = {
      if (!requireNamespace("ranger", quietly = TRUE)) {
        stop("Package 'ranger' required for method = 'ranger'.", call. = FALSE)
      }
      df_train <- data.frame(y = factor(treatment), covariates)
      fit <- ranger::ranger(y ~ ., data = df_train, probability = TRUE, num.trees = 500)
      df_new <- data.frame(newdata)
      names(df_new) <- names(df_train)[-1]
      predict(fit, data = df_new)$predictions[, "1"]
    },
    xgboost = {
      if (!requireNamespace("xgboost", quietly = TRUE)) {
        stop("Package 'xgboost' required for method = 'xgboost'.", call. = FALSE)
      }
      dtrain <- xgboost::xgb.DMatrix(data = covariates, label = treatment)
      fit <- xgboost::xgb.train(
        params = list(objective = "binary:logistic", max_depth = 4, eta = 0.1),
        data = dtrain,
        nrounds = 100,
        verbose = 0
      )
      predict(fit, xgboost::xgb.DMatrix(data = newdata))
    }
  )
}

#' @export
print.propensity_fit <- function(x, ...) {
  cat("Propensity score estimation (", x$method, ")\n")
  cat("  N:       ", x$n, "\n")
  cat("  Trimmed: ", x$n_trimmed, " obs to [",
    x$trim, ", ", 1 - x$trim, "]\n",
    sep = ""
  )
  cat("  Range:   [",
    formatC(min(x$scores), digits = 3, format = "f"), ", ",
    formatC(max(x$scores), digits = 3, format = "f"), "]\n",
    sep = ""
  )
  cat("  Mean:    ", formatC(mean(x$scores), digits = 3, format = "f"), "\n")
  invisible(x)
}
