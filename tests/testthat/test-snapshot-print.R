test_that("print.kernel_spec is stable", {
  expect_snapshot(print(kernel_spec()))
  expect_snapshot(print(kernel_spec("rbf", bandwidth = 1.5)))
  expect_snapshot(print(kernel_spec("matern", bandwidth = 1, nu = 1.5)))
  expect_snapshot(print(
    kernel_spec("rbf", approx = "nystrom", approx_rank = 50, approx_seed = 1)
  ))
  expect_snapshot(print(
    kernel_spec("rbf", approx = "rff", approx_rank = 200, approx_seed = 7)
  ))
})

test_that("print.nystrom_factor is stable", {
  set.seed(1)
  x <- matrix(rnorm(40), 20, 2)
  nf <- nystrom_factor(x, kernel = kernel_spec("rbf", bandwidth = 1),
                       m = 10, seed = 42)
  expect_snapshot(print(nf))
})

test_that("print.rff_features is stable", {
  set.seed(1)
  x <- matrix(rnorm(40), 20, 2)
  rf <- rff_features(x, kernel = kernel_spec("rbf", bandwidth = 1),
                     D = 25, seed = 42)
  expect_snapshot(print(rf))
})
