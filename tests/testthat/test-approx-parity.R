# Parity tests for the opt-in approximations.
#
# Strategy: for each (test, kernel, approximation) tuple we compare the
# approximated kernel's statistic against the exact statistic. Where the
# approximation rank is close to the sample size, the approximation
# error should be small in relative terms; where the rank is moderate
# (~0.1n) we expect "close but not identical" — captured with a
# loosened tolerance reflecting the documented rule of thumb.
#
# We intentionally compare *statistics* rather than permutation p-values:
# permutation p-values incorporate a layer of Monte Carlo noise on top
# of the kernel approximation, which would conflate the two sources of
# variation and turn pass/fail into a flaky question of seed choice.

# ------------------------------------------------------------------ #
# Helper: relative Frobenius error.
relfro <- function(A, B) sqrt(sum((A - B)^2) / sum(A^2))

# ------------------------------------------------------------------ #
# Kernel matrix parity (the foundation that the tests below rely on).
test_that("Nystrom matches exact RBF kernel within tolerance at near-full rank", {
  set.seed(1)
  x <- matrix(rnorm(160), 80, 2)
  k_exact <- kernel_spec("rbf", bandwidth = 1)
  k_nys   <- kernel_spec("rbf", bandwidth = 1,
                         approx = "nystrom", approx_rank = 70, approx_seed = 1)

  K_exact <- kernel_matrix(x, kernel = k_exact)
  K_nys   <- kernel_matrix(x, kernel = k_nys)
  expect_lt(relfro(K_exact, K_nys), 0.02)
})

test_that("RFF matches exact RBF kernel as feature count grows", {
  set.seed(2)
  x <- matrix(rnorm(160), 80, 2)
  k_exact <- kernel_spec("rbf", bandwidth = 1)

  err_low_d <- {
    k <- kernel_spec("rbf", bandwidth = 1,
                     approx = "rff", approx_rank = 100, approx_seed = 1)
    relfro(kernel_matrix(x, kernel = k_exact), kernel_matrix(x, kernel = k))
  }
  err_high_d <- {
    k <- kernel_spec("rbf", bandwidth = 1,
                     approx = "rff", approx_rank = 4000, approx_seed = 1)
    relfro(kernel_matrix(x, kernel = k_exact), kernel_matrix(x, kernel = k))
  }
  expect_lt(err_high_d, err_low_d)
  expect_lt(err_high_d, 0.05)
})

# ------------------------------------------------------------------ #
# HSIC statistic parity.
test_that("HSIC statistic agrees across exact / Nystrom / RFF (independent x, y)", {
  set.seed(11)
  n <- 100
  x <- matrix(rnorm(n * 2), n, 2)
  y <- matrix(rnorm(n * 2), n, 2)

  res_exact <- hsic_test(x, y, n_permutations = 50L, seed = 1)
  res_nys <- hsic_test(
    x, y,
    kernel_x = kernel_spec("rbf", approx = "nystrom", approx_rank = 90, approx_seed = 1),
    kernel_y = kernel_spec("rbf", approx = "nystrom", approx_rank = 90, approx_seed = 2),
    n_permutations = 50L, seed = 1
  )
  res_rff <- hsic_test(
    x, y,
    kernel_x = kernel_spec("rbf", approx = "rff", approx_rank = 1500, approx_seed = 1),
    kernel_y = kernel_spec("rbf", approx = "rff", approx_rank = 1500, approx_seed = 2),
    n_permutations = 50L, seed = 1
  )

  rel_nys <- abs(res_nys$statistic - res_exact$statistic) /
    max(abs(res_exact$statistic), 1e-12)
  rel_rff <- abs(res_rff$statistic - res_exact$statistic) /
    max(abs(res_exact$statistic), 1e-12)

  expect_lt(rel_nys, 0.10)
  expect_lt(rel_rff, 0.30)
})

# ------------------------------------------------------------------ #
# MMD statistic parity (two-sample, modest mean shift).
test_that("MMD statistic agrees across exact / Nystrom / RFF", {
  set.seed(12)
  n <- 80
  x <- matrix(rnorm(n * 2), n, 2)
  y <- matrix(rnorm(n * 2, mean = 0.4), n, 2)

  res_exact <- mmd_test(x, y, n_permutations = 50L, seed = 1)
  res_nys <- mmd_test(
    x, y,
    kernel = kernel_spec("rbf", approx = "nystrom", approx_rank = 70, approx_seed = 1),
    n_permutations = 50L, seed = 1
  )
  res_rff <- mmd_test(
    x, y,
    kernel = kernel_spec("rbf", approx = "rff", approx_rank = 1500, approx_seed = 1),
    n_permutations = 50L, seed = 1
  )

  rel_nys <- abs(res_nys$statistic - res_exact$statistic) /
    max(abs(res_exact$statistic), 1e-12)
  rel_rff <- abs(res_rff$statistic - res_exact$statistic) /
    max(abs(res_exact$statistic), 1e-12)

  expect_lt(rel_nys, 0.10)
  expect_lt(rel_rff, 0.30)
})

# ------------------------------------------------------------------ #
# Bd-HSIC integration: with full-rank approximation, weighted statistic
# should be close to the exact path.
test_that("bd-HSIC runs end-to-end with Nystrom kernel and produces sensible output", {
  set.seed(13)
  n <- 100
  x <- matrix(rnorm(n * 2), n, 2)
  z <- matrix(rnorm(n * 2), n, 2)
  y <- 0.5 * z[, 1] + rnorm(n, sd = 0.5)

  # bd_hsic_test cross-fits with a 0.5 split, so each fold has n/2 = 50 rows.
  # `approx_rank` must fit the smaller fold.
  k_nys <- kernel_spec("rbf", approx = "nystrom", approx_rank = 30, approx_seed = 1)

  res <- bd_hsic_test(
    x = x, y = y, z = z,
    kernel_x = k_nys, kernel_y = k_nys,
    n_permutations = 50L, seed = 1
  )

  expect_s3_class(res, "kernel_test_result")
  expect_true(is.numeric(res$statistic))
  expect_true(is.numeric(res$p_value))
  expect_gte(res$p_value, 0)
  expect_lte(res$p_value, 1)
})

# ------------------------------------------------------------------ #
# Low-overlap fixture: under near-perfect propensity separation,
# Nystrom approximation should still produce a finite, non-NA result.
test_that("DR-DATE survives a low-overlap fixture under Nystrom kernel", {
  set.seed(14)
  n <- 120
  x <- matrix(rnorm(n * 2), n, 2)
  ps <- plogis(2.5 * x[, 1])  # strong but not perfect separation
  t  <- stats::rbinom(n, 1, ps)
  y  <- 0.3 * t + 0.5 * x[, 1] + rnorm(n, sd = 0.5)

  # dr_date_test also cross-fits; approx_rank must fit the post-split folds.
  k_nys <- kernel_spec("rbf", approx = "nystrom", approx_rank = 30, approx_seed = 1)

  res <- dr_date_test(
    y = y, treatment = t, covariates = x,
    kernel_y = k_nys,
    n_permutations = 50L, seed = 1
  )

  expect_s3_class(res, "kernel_test_result")
  expect_true(is.finite(res$statistic))
  expect_true(is.finite(res$p_value))
  expect_gte(res$p_value, 0)
  expect_lte(res$p_value, 1)
})
