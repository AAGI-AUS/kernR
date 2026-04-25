test_that("rff_features returns Z with the expected shape", {
  set.seed(1)
  x <- matrix(rnorm(200), 100, 2)
  rf <- rff_features(x, kernel = kernel_spec("rbf", bandwidth = 1),
                     D = 50, seed = 1)
  expect_s3_class(rf, "rff_features")
  expect_equal(nrow(rf$Z), 100)
  expect_equal(ncol(rf$Z), 2L * 50L)
  expect_equal(rf$D, 50L)
  expect_equal(dim(rf$frequencies), c(50L, 2L))
})

test_that("rff reconstruction is symmetric and PSD", {
  set.seed(2)
  x <- matrix(rnorm(120), 60, 2)
  rf <- rff_features(x, kernel = kernel_spec("rbf", bandwidth = 1.5),
                     D = 100, seed = 7)
  K_hat <- rff_reconstruct(rf)

  expect_equal(dim(K_hat), c(60, 60))
  expect_true(isSymmetric(K_hat, tol = 1e-10))
  ev <- eigen(K_hat, symmetric = TRUE, only.values = TRUE)$values
  expect_gte(min(ev), -1e-8)
})

test_that("rff approximates the exact RBF kernel as D grows", {
  set.seed(3)
  n <- 80
  x <- matrix(rnorm(n * 2), n, 2)
  k_exact <- kernel_spec("rbf", bandwidth = 1)
  K <- kernel_matrix(x, kernel = k_exact)

  err_small <- {
    rf <- rff_features(x, kernel = k_exact, D = 50, seed = 1)
    sqrt(sum((K - rff_reconstruct(rf))^2) / sum(K^2))
  }
  err_large <- {
    rf <- rff_features(x, kernel = k_exact, D = 2000, seed = 1)
    sqrt(sum((K - rff_reconstruct(rf))^2) / sum(K^2))
  }

  # Larger D should give a tighter approximation.
  expect_lt(err_large, err_small)
  # And with D = 2000 the relative error should be modest.
  expect_lt(err_large, 0.1)
})

test_that("rff is consistent with kernel_matrix() approx dispatch", {
  set.seed(4)
  x <- matrix(rnorm(80), 40, 2)
  k <- kernel_spec("rbf", bandwidth = 1.2,
                   approx = "rff", approx_rank = 100, approx_seed = 99)

  K_via_dispatch <- kernel_matrix(x, kernel = k)
  rf <- rff_features(x, kernel = k, D = 100, seed = 99)
  K_direct <- rff_reconstruct(rf)

  expect_equal(K_via_dispatch, K_direct)
})

test_that("rff falls back to exact for rectangular kernel_matrix calls", {
  set.seed(5)
  x <- matrix(rnorm(60), 30, 2)
  y <- matrix(rnorm(40), 20, 2)
  k <- kernel_spec("rbf", bandwidth = 1,
                   approx = "rff", approx_rank = 50, approx_seed = 1)

  K_rect <- kernel_matrix(x, y, kernel = k)
  K_exact <- kernel_matrix(x, y, kernel = kernel_spec("rbf", bandwidth = 1))
  expect_equal(K_rect, K_exact)
})

test_that("rff seed gives reproducible features and frequencies", {
  x <- matrix(rnorm(100), 50, 2)
  k <- kernel_spec("rbf", bandwidth = 1)
  rf1 <- rff_features(x, kernel = k, D = 30, seed = 42)
  rf2 <- rff_features(x, kernel = k, D = 30, seed = 42)
  expect_equal(rf1$frequencies, rf2$frequencies)
  expect_equal(rf1$Z, rf2$Z)
})

test_that("rff seed does not leak into the global RNG state", {
  set.seed(11)
  x <- matrix(rnorm(40), 20, 2)
  k <- kernel_spec("rbf", bandwidth = 1)

  set.seed(123)
  before <- runif(1)

  set.seed(123)
  rff_features(x, kernel = k, D = 30, seed = 7)
  after <- runif(1)

  expect_equal(before, after)
})

test_that("rff matern (finite nu) approximates the exact matern kernel", {
  set.seed(6)
  n <- 60
  x <- matrix(rnorm(n * 2), n, 2)

  for (nu in c(1.5, 2.5)) {
    k_exact <- kernel_spec("matern", bandwidth = 1, nu = nu)
    K <- kernel_matrix(x, kernel = k_exact)

    rf <- rff_features(x, kernel = k_exact, D = 2000, seed = 1)
    K_hat <- rff_reconstruct(rf)
    rel_err <- sqrt(sum((K - K_hat)^2) / sum(K^2))
    expect_lt(rel_err, 0.15)
  }
})

test_that("rff matern with nu = Inf agrees with rbf", {
  set.seed(7)
  x <- matrix(rnorm(60), 30, 2)

  k_rbf <- kernel_spec("rbf", bandwidth = 1)
  k_matern_inf <- kernel_spec("matern", bandwidth = 1, nu = Inf)

  rf_rbf <- rff_features(x, kernel = k_rbf, D = 100, seed = 1)
  rf_inf <- rff_features(x, kernel = k_matern_inf, D = 100, seed = 1)

  # Same seed → same Gaussian draw → same frequencies → same features.
  expect_equal(rf_rbf$frequencies, rf_inf$frequencies)
  expect_equal(rf_rbf$Z, rf_inf$Z)
})

test_that("rff rejects non-translation-invariant kernels", {
  x <- matrix(rnorm(40), 20, 2)
  expect_error(
    rff_features(x, kernel = kernel_spec("linear"), D = 30),
    "translation-invariant"
  )
  expect_error(
    rff_features(x, kernel = kernel_spec("polynomial"), D = 30),
    "translation-invariant"
  )

  # Same restriction holds for kernel_matrix dispatch.
  k_lin <- kernel_spec("linear", approx = "rff", approx_rank = 30)
  expect_error(kernel_matrix(x, kernel = k_lin), "translation-invariant")
})

test_that("rff_features argument validation", {
  x <- matrix(rnorm(40), 20, 2)
  k <- kernel_spec("rbf", bandwidth = 1)
  expect_error(rff_features(x, k, D = 0), "at least 1")
  expect_error(rff_features(matrix(numeric(0), 0, 2), k), "zero rows")
})

test_that("default_rff_rank matches the Nystrom heuristic shape", {
  expect_equal(default_rff_rank(20), 20L)
  expect_equal(default_rff_rank(1000), 100L)
  expect_equal(default_rff_rank(40), 40L)
})
