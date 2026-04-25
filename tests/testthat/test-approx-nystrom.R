test_that("nystrom_factor returns U with the expected shape", {
  set.seed(1)
  x <- matrix(rnorm(100), 50, 2)
  nf <- nystrom_factor(x, kernel = kernel_spec("rbf", bandwidth = 1),
                       m = 20, seed = 1)
  expect_s3_class(nf, "nystrom_factor")
  expect_equal(nrow(nf$U), 50)
  expect_lte(ncol(nf$U), 20)
  expect_equal(nf$m, 20L)
  expect_length(nf$landmarks, 20L)
  expect_true(all(nf$landmarks %in% seq_len(50)))
  expect_false(anyDuplicated(nf$landmarks) > 0)
})

test_that("nystrom reconstruction is symmetric and PSD", {
  set.seed(2)
  x <- matrix(rnorm(120), 60, 2)
  nf <- nystrom_factor(x, kernel = kernel_spec("rbf", bandwidth = 1.5),
                       m = 25, seed = 7)
  K_hat <- nystrom_reconstruct(nf)

  expect_equal(dim(K_hat), c(60, 60))
  expect_true(isSymmetric(K_hat, tol = 1e-10))
  ev <- eigen(K_hat, symmetric = TRUE, only.values = TRUE)$values
  expect_gte(min(ev), -1e-8)
})

test_that("nystrom approximates the exact RBF kernel as m approaches n", {
  set.seed(3)
  n <- 80
  x <- matrix(rnorm(n * 2), n, 2)
  k_exact <- kernel_spec("rbf", bandwidth = 1)
  K <- kernel_matrix(x, kernel = k_exact)

  # With m = n the Nystrom approximation should be exact (up to numerical
  # noise from eigendecomp + reconstruction).
  nf_full <- nystrom_factor(x, kernel = k_exact, m = n, seed = 1)
  K_full <- nystrom_reconstruct(nf_full)
  expect_lt(max(abs(K - K_full)), 1e-8)

  # With m around 0.5n the approximation should be close in Frobenius norm.
  nf_half <- nystrom_factor(x, kernel = k_exact, m = 40, seed = 1)
  K_half <- nystrom_reconstruct(nf_half)
  rel_err <- sqrt(sum((K - K_half)^2) / sum(K^2))
  expect_lt(rel_err, 0.05)
})

test_that("nystrom is consistent with kernel_matrix() approx dispatch", {
  set.seed(4)
  x <- matrix(rnorm(80), 40, 2)
  k <- kernel_spec("rbf", bandwidth = 1.2,
                   approx = "nystrom", approx_rank = 20, approx_seed = 99)

  K_via_dispatch <- kernel_matrix(x, kernel = k)
  nf <- nystrom_factor(x, kernel = k, m = 20, seed = 99)
  K_direct <- nystrom_reconstruct(nf)

  expect_equal(K_via_dispatch, K_direct)
})

test_that("kernel_matrix() falls back to exact for rectangular calls", {
  set.seed(5)
  x <- matrix(rnorm(60), 30, 2)
  y <- matrix(rnorm(40), 20, 2)
  k <- kernel_spec("rbf", bandwidth = 1,
                   approx = "nystrom", approx_rank = 10, approx_seed = 1)

  K_rect <- kernel_matrix(x, y, kernel = k)
  K_exact <- kernel_matrix(x, y, kernel = kernel_spec("rbf", bandwidth = 1))
  expect_equal(K_rect, K_exact)
})

test_that("nystrom_factor seed gives reproducible landmarks and factor", {
  x <- matrix(rnorm(100), 50, 2)
  k <- kernel_spec("rbf", bandwidth = 1)
  nf1 <- nystrom_factor(x, kernel = k, m = 15, seed = 42)
  nf2 <- nystrom_factor(x, kernel = k, m = 15, seed = 42)
  expect_equal(nf1$landmarks, nf2$landmarks)
  expect_equal(nf1$U, nf2$U)
})

test_that("nystrom seed does not leak into the global RNG state", {
  set.seed(11)
  x <- matrix(rnorm(40), 20, 2)
  k <- kernel_spec("rbf", bandwidth = 1)

  set.seed(123)
  before <- runif(1)

  set.seed(123)
  nystrom_factor(x, kernel = k, m = 5, seed = 7)
  after <- runif(1)

  expect_equal(before, after)
})

test_that("nystrom works for matern, linear, and polynomial kernels", {
  set.seed(6)
  x <- matrix(rnorm(60), 30, 2)
  for (k in list(
    kernel_spec("matern", bandwidth = 1, nu = 1.5),
    kernel_spec("matern", bandwidth = 1, nu = 2.5),
    kernel_spec("linear"),
    kernel_spec("polynomial", degree = 2L, offset = 1)
  )) {
    nf <- nystrom_factor(x, kernel = k, m = 20, seed = 1)
    K_hat <- nystrom_reconstruct(nf)
    expect_equal(dim(K_hat), c(30, 30))
    expect_true(isSymmetric(K_hat, tol = 1e-8))
  }
})

test_that("nystrom_factor argument validation", {
  x <- matrix(rnorm(40), 20, 2)
  k <- kernel_spec("rbf", bandwidth = 1)
  expect_error(nystrom_factor(x, k, m = 0), "at least 1")
  expect_error(nystrom_factor(x, k, m = 25), "cannot exceed")
  expect_error(nystrom_factor(x, k, m = 5, tol = -1), "non-negative")
  expect_error(nystrom_factor(matrix(numeric(0), 0, 2), k), "zero rows")
})

test_that("kernel_spec validates approx arguments", {
  expect_error(kernel_spec(approx = "bogus"), "should be one of")
  expect_error(kernel_spec(approx = "nystrom", approx_rank = 0),
               "positive integer")
  expect_error(kernel_spec(approx = "nystrom", approx_rank = -3),
               "positive integer")
  expect_error(kernel_spec(approx = "nystrom", approx_seed = c(1, 2)),
               "single integer")

  k <- kernel_spec("rbf", approx = "nystrom", approx_rank = 50, approx_seed = 7)
  expect_equal(k$approx, "nystrom")
  expect_equal(k$approx_rank, 50L)
  expect_equal(k$approx_seed, 7L)
})

test_that("default_nystrom_rank caps at n", {
  expect_equal(default_nystrom_rank(20), 20L)
  expect_equal(default_nystrom_rank(1000), 100L)
  expect_equal(default_nystrom_rank(40), 40L)
})
