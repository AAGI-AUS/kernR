#' Unified Kernel Causal Test
#'
#' A convenience wrapper that dispatches to the appropriate test function
#' based on the `method` argument.
#'
#' @param formula A formula of the form `y ~ treatment | confounders`.
#' @param data A data.frame or data.table.
#' @param method Character. Test method: `"dr-date"` (default), `"dr-dett"`,
#'   or `"bd-hsic"`.
#' @param ... Additional arguments passed to the specific test function.
#'   This includes the sequential / adaptive permutation controls
#'   `adaptive`, `B_max`, and `batch_size` — see [dr_date_test()],
#'   [dr_dett_test()], and [bd_hsic_test()] for details.
#'
#' @return An object of class `"kernel_test_result"`.
#'
#' @family tests-causal
#' @seealso [bd_hsic_test()], [dr_date_test()], [dr_dett_test()] for the
#'   underlying tests; the formula syntax dispatches to one of these.
#'
#' @examples
#' set.seed(42)
#' n <- 200
#' dat <- data.frame(
#'   y = rnorm(n),
#'   treatment = rbinom(n, 1, 0.5),
#'   x1 = rnorm(n),
#'   x2 = rnorm(n)
#' )
#' dat$y <- dat$y + 0.5 * dat$treatment + 0.3 * dat$x1
#'
#' result <- kernel_causal_test(y ~ treatment | x1 + x2,
#'   data = dat, method = "dr-date",
#'   n_permutations = 100, seed = 1
#' )
#' print(result)
#'
#' @export
kernel_causal_test <- function(formula, data,
                               method = c("dr-date", "dr-dett", "bd-hsic"),
                               ...) {
  method <- match.arg(method)

  # Parse formula: y ~ treatment | confounders
  parsed <- parse_causal_formula(formula, data)

  switch(method,
    "dr-date" = dr_date_test(
      y = parsed$y,
      treatment = parsed$treatment,
      covariates = parsed$covariates,
      ...
    ),
    "dr-dett" = dr_dett_test(
      y = parsed$y,
      treatment = parsed$treatment,
      covariates = parsed$covariates,
      ...
    ),
    "bd-hsic" = bd_hsic_test(
      x = parsed$treatment,
      y = parsed$y,
      z = parsed$covariates,
      ...
    )
  )
}

#' Parse Causal Formula
#'
#' Parses `y ~ treatment | confounders` into component matrices.
#'
#' @param formula Formula.
#' @param data data.frame.
#'
#' @return List with `y`, `treatment`, `covariates`.
#' @keywords internal
parse_causal_formula <- function(formula, data) {
  # Convert formula to string and split on |
  fstr <- deparse(formula, width.cutoff = 500)

  if (!grepl("\\|", fstr)) {
    stop("Formula must include confounders separated by '|', ",
      "e.g., y ~ treatment | x1 + x2",
      call. = FALSE
    )
  }

  parts <- strsplit(fstr, "\\|")[[1]]
  lhs_rhs <- strsplit(trimws(parts[1]), "~")[[1]]

  y_name <- trimws(lhs_rhs[1])
  t_name <- trimws(lhs_rhs[2])
  conf_formula <- stats::as.formula(paste("~", trimws(parts[2])))

  y <- as.matrix(data[[y_name]])
  treatment <- data[[t_name]]
  covariates <- stats::model.matrix(conf_formula, data = data)[, -1, drop = FALSE]

  list(y = y, treatment = treatment, covariates = covariates)
}
