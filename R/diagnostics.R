#' Assess Propensity Score Overlap
#'
#' Diagnoses overlap (positivity) between treated and control groups
#' by summarising the propensity score distributions.
#'
#' @param propensity A `propensity_fit` object or a numeric vector of scores.
#' @param treatment Binary treatment vector. Required if `propensity` is a
#'   numeric vector.
#'
#' @return A list of class `"overlap_diagnostic"` with:
#'   \describe{
#'     \item{treated}{Summary statistics of propensity scores for treated.}
#'     \item{control}{Summary statistics for controls.}
#'     \item{overlap_warning}{Logical. TRUE if overlap is poor.}
#'   }
#'
#' @family diagnostics
#' @seealso [estimate_propensity()], [plot_weights()],
#'   [effective_sample_size()].
#'
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(200), 100, 2)
#' t <- rbinom(100, 1, plogis(0.5 * x[, 1]))
#' ps <- estimate_propensity(t, x, method = "logistic")
#' assess_overlap(ps, t)
#'
#' @export
assess_overlap <- function(propensity, treatment = NULL) {
  if (inherits(propensity, "propensity_fit")) {
    scores <- propensity$scores
  } else {
    scores <- as.numeric(propensity)
  }

  if (is.null(treatment)) {
    stop("`treatment` is required.", call. = FALSE)
  }
  treatment <- as.integer(treatment)

  s1 <- scores[treatment == 1]
  s0 <- scores[treatment == 0]

  treated_summary <- c(
    min = min(s1), q25 = unname(quantile(s1, 0.25)),
    median = median(s1), q75 = unname(quantile(s1, 0.75)),
    max = max(s1)
  )
  control_summary <- c(
    min = min(s0), q25 = unname(quantile(s0, 0.25)),
    median = median(s0), q75 = unname(quantile(s0, 0.75)),
    max = max(s0)
  )

  # Overlap warning: if ranges barely intersect
  overlap_min <- max(min(s1), min(s0))
  overlap_max <- min(max(s1), max(s0))
  overlap_frac <- max(0, overlap_max - overlap_min) /
    (max(max(s1), max(s0)) - min(min(s1), min(s0)) + 1e-10)

  structure(
    list(
      treated = treated_summary,
      control = control_summary,
      overlap_fraction = overlap_frac,
      overlap_warning = overlap_frac < 0.5
    ),
    class = "overlap_diagnostic"
  )
}

#' @export
print.overlap_diagnostic <- function(x, ...) {
  cat("Propensity Score Overlap Diagnostic\n\n")
  cat(
    "Treated:  ",
    paste(names(x$treated), "=", formatC(x$treated, digits = 3, format = "f"), collapse = ", "),
    "\n"
  )
  cat(
    "Control:  ",
    paste(names(x$control), "=", formatC(x$control, digits = 3, format = "f"), collapse = ", "),
    "\n"
  )
  cat("Overlap:  ", formatC(x$overlap_fraction * 100, digits = 1, format = "f"), "%\n")
  if (x$overlap_warning) {
    cat("WARNING: Poor overlap detected. Consider using DETT (requires only one-sided overlap).\n")
  }
  invisible(x)
}

#' Plot Weight Diagnostics
#'
#' Plots the distribution of importance weights with effective sample
#' size annotation.
#'
#' @param weights Numeric vector of importance weights.
#' @param main Title. Default is "Weight Distribution".
#'
#' @return Invisibly returns `weights`.
#'
#' @family diagnostics
#' @seealso [assess_overlap()], [effective_sample_size()],
#'   [estimate_density_ratio()].
#'
#' @examples
#' set.seed(1)
#' w <- exp(rnorm(200))
#' plot_weights(w)
#'
#' @export
plot_weights <- function(weights, main = "Weight Distribution") {
  ess <- effective_sample_size(weights)
  n <- length(weights)

  hist(weights,
    breaks = 30,
    main = main,
    xlab = "Weight",
    col = "grey80",
    border = "grey60",
    freq = FALSE
  )
  abline(v = 1, col = "#0072B2", lwd = 2, lty = 2)
  legend("topright",
    legend = c(
      paste("ESS =", formatC(ess, digits = 1, format = "f"),
        "/", n
      ),
      paste("ESS ratio =", formatC(ess / n * 100, digits = 1, format = "f"), "%")
    ),
    bty = "n"
  )

  invisible(weights)
}
