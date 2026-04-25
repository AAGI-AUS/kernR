#-----------------------------------------------------------------------
# Tests for adaptive (sequential) permutation procedure (Feature B4).
# Covers the internal helper `adaptive_permutation_pvalue()` and the
# `adaptive = TRUE` branch of the public test functions, plus
# byte-identity of the `adaptive = FALSE` default.
#-----------------------------------------------------------------------

# --- Helpers -----------------------------------------------------------

# Make a deterministic perm_fun whose draws come from a fixed pool.
# Useful to verify monotonicity and identity properties without
# stochasticity from C++ randperm.
make_pool_perm_fun <- function(pool) {
  pos <- 0L
  function(n_perms) {
    out <- pool[(pos + 1L):(pos + n_perms)]
    pos <<- pos + as.integer(n_perms)
    out
  }
}


# --- adaptive_permutation_pvalue: identity and exhaustion --------------

test_that("adaptive_permutation_pvalue matches fixed full run when no early stop", {
  # B_min == B_max forces exhaustion: the budget is consumed on the
  # first batch and no stopping check is consulted before the loop ends.
  # The reported p-value must equal the standard +1 fixed-B formula.
  set.seed(11)
  pool <- rnorm(200)
  observed <- 0
  pf <- make_pool_perm_fun(pool)

  n_extreme <- sum(pool >= observed)
  manual_p <- (1 + n_extreme) / (1 + length(pool))

  res <- kernR:::adaptive_permutation_pvalue(
    observed = observed,
    perm_fun = pf,
    B_min = 200L,
    B_max = 200L,
    batch_size = 200L,
    alpha = 0.05
  )

  expect_equal(res$n_perms_used, 200L)
  # When B_min == B_max the helper still computes the same +1
  # corrected p-value as the fixed run; stop_reason may legitimately
  # be either an early-stop label or "max_reached" depending on the
  # Wilson interval — we only require correctness of `p_value`.
  expect_equal(res$p_value, manual_p)
  expect_length(res$null_distribution, 200L)
})


# --- Stop reasons ------------------------------------------------------

test_that("adaptive stops with `upper_below_alpha` under strong H1", {
  # All null draws far below observed -> Wilson upper rapidly < 0.05.
  observed <- 100
  pf <- make_pool_perm_fun(rep(0, 5000))
  res <- kernR:::adaptive_permutation_pvalue(
    observed = observed,
    perm_fun = pf,
    B_min = 99L,
    B_max = 5000L,
    batch_size = 100L,
    alpha = 0.05
  )
  expect_equal(res$stop_reason, "upper_below_alpha")
  expect_lt(res$n_perms_used, 5000L)
})

test_that("adaptive stops with `lower_above_alpha` under clear H0", {
  # All null draws above observed -> p_hat ~ 1 -> Wilson lower > 0.05.
  observed <- -100
  pf <- make_pool_perm_fun(rep(0, 5000))
  res <- kernR:::adaptive_permutation_pvalue(
    observed = observed,
    perm_fun = pf,
    B_min = 99L,
    B_max = 5000L,
    batch_size = 100L,
    alpha = 0.05
  )
  expect_equal(res$stop_reason, "lower_above_alpha")
  expect_lt(res$n_perms_used, 5000L)
})


# --- Monotonicity ------------------------------------------------------

test_that("monotonicity: increasing B_max does not change p once stopped early", {
  observed <- 100
  pf_a <- make_pool_perm_fun(rep(0, 20000))
  pf_b <- make_pool_perm_fun(rep(0, 20000))

  res_small <- kernR:::adaptive_permutation_pvalue(
    observed = observed, perm_fun = pf_a,
    B_min = 99L, B_max = 1000L,
    batch_size = 100L, alpha = 0.05
  )
  res_big <- kernR:::adaptive_permutation_pvalue(
    observed = observed, perm_fun = pf_b,
    B_min = 99L, B_max = 20000L,
    batch_size = 100L, alpha = 0.05
  )

  expect_equal(res_small$stop_reason, "upper_below_alpha")
  expect_equal(res_big$stop_reason, "upper_below_alpha")
  expect_equal(res_small$p_value, res_big$p_value)
  expect_equal(res_small$n_perms_used, res_big$n_perms_used)
})


# --- Argument validation ----------------------------------------------

test_that("adaptive_permutation_pvalue validates arguments", {
  pf <- function(n) rnorm(n)

  expect_error(
    kernR:::adaptive_permutation_pvalue(0, pf, alpha = 0),
    "alpha"
  )
  expect_error(
    kernR:::adaptive_permutation_pvalue(0, pf, alpha = 1),
    "alpha"
  )
  expect_error(
    kernR:::adaptive_permutation_pvalue(0, pf, alpha = -0.1),
    "alpha"
  )
  expect_error(
    kernR:::adaptive_permutation_pvalue(0, pf, B_min = 1L),
    "B_min"
  )
  expect_error(
    kernR:::adaptive_permutation_pvalue(0, pf, B_min = 100L, B_max = 50L),
    "B_max"
  )
  expect_error(
    kernR:::adaptive_permutation_pvalue(0, pf, batch_size = 0L),
    "batch_size"
  )
  expect_error(
    kernR:::adaptive_permutation_pvalue(0, pf, wilson_level = 0),
    "wilson_level"
  )
  expect_error(
    kernR:::adaptive_permutation_pvalue(0, "not a function"),
    "perm_fun"
  )
  expect_error(
    kernR:::adaptive_permutation_pvalue(c(1, 2), pf),
    "observed"
  )
})

test_that("adaptive_permutation_pvalue errors when perm_fun returns wrong length", {
  bad_pf <- function(n) rnorm(n + 1)
  expect_error(
    kernR:::adaptive_permutation_pvalue(0, bad_pf,
      B_min = 50L, B_max = 200L, batch_size = 50L
    ),
    "perm_fun"
  )
})


# --- HSIC adaptive: rejection rate under H0 ---------------------------

test_that("adaptive HSIC controls type-I error under H0", {
  # 100 reps under H0 (independent X, Y) — rejection rate at alpha=0.05
  # should be roughly nominal. Wide tolerance for 100 reps.
  skip_on_cran()
  reps <- 100
  rejections <- 0L
  for (r in seq_len(reps)) {
    set.seed(1000 + r)
    n <- 60
    x <- rnorm(n)
    y <- rnorm(n)
    res <- hsic_test(x, y,
      adaptive = TRUE, B_max = 999L, batch_size = 100L,
      alpha = 0.05, seed = 1000 + r
    )
    if (res$p_value < 0.05) rejections <- rejections + 1L
  }
  expect_lte(rejections / reps, 0.10)
})


# --- HSIC adaptive: stops early under strong H1 -----------------------

test_that("adaptive HSIC stops early under strong dependence", {
  set.seed(7)
  n <- 200
  x <- rnorm(n)
  y <- x + rnorm(n, sd = 0.2)  # near-perfect dependence
  res <- hsic_test(x, y,
    adaptive = TRUE,
    B_max = 5000L, batch_size = 100L,
    alpha = 0.05, seed = 7
  )
  expect_lt(res$n_perms_used, 5000L)
  expect_true(res$p_value < 0.05)
  expect_equal(res$stop_reason, "upper_below_alpha")
})


# --- Byte identity when adaptive = FALSE ------------------------------

test_that("hsic_test with adaptive = FALSE is byte-identical to default", {
  set.seed(1)
  n <- 80
  x <- rnorm(n)
  y <- x^2 + rnorm(n, sd = 0.5)

  res_default <- hsic_test(x, y, n_permutations = 200, seed = 1)
  res_explicit <- hsic_test(x, y,
    n_permutations = 200, seed = 1,
    adaptive = FALSE
  )

  expect_equal(res_default$statistic, res_explicit$statistic)
  expect_equal(res_default$p_value, res_explicit$p_value)
  expect_equal(res_default$null_distribution, res_explicit$null_distribution)
})

test_that("mmd_test with adaptive = FALSE is byte-identical to default", {
  set.seed(1)
  x <- matrix(rnorm(120), 60, 2)
  y <- matrix(rnorm(120, mean = 0.5), 60, 2)

  res_default <- mmd_test(x, y, n_permutations = 150, seed = 1)
  res_explicit <- mmd_test(x, y,
    n_permutations = 150, seed = 1,
    adaptive = FALSE
  )

  expect_equal(res_default$statistic, res_explicit$statistic)
  expect_equal(res_default$p_value, res_explicit$p_value)
  expect_equal(res_default$null_distribution, res_explicit$null_distribution)
})


# --- MMD adaptive smoke test ------------------------------------------

test_that("mmd_test with adaptive = TRUE returns expected fields", {
  set.seed(3)
  x <- matrix(rnorm(200), 100, 2)
  y <- matrix(rnorm(200, mean = 1.5), 100, 2)
  res <- mmd_test(x, y,
    adaptive = TRUE, B_max = 999L, batch_size = 100L,
    alpha = 0.05, seed = 3
  )
  expect_s3_class(res, "kernel_test_result")
  expect_true(!is.null(res$n_perms_used))
  expect_true(res$stop_reason %in% c("upper_below_alpha", "lower_above_alpha", "max_reached"))
  expect_lt(res$p_value, 0.05)
})


# --- bd-HSIC adaptive smoke test --------------------------------------

test_that("bd_hsic_test with adaptive = TRUE runs and returns adaptive fields", {
  skip_on_cran()
  set.seed(5)
  n <- 150
  z <- matrix(rnorm(n * 2), n, 2)
  x <- z[, 1] + rnorm(n)
  y <- 0.8 * x + z[, 2] + rnorm(n, sd = 0.3)

  res <- bd_hsic_test(x, y, z,
    adaptive = TRUE, B_max = 999L, batch_size = 100L,
    alpha = 0.05, seed = 5
  )
  expect_s3_class(res, "kernel_test_result")
  expect_true(!is.null(res$n_perms_used))
  expect_true(res$stop_reason %in% c("upper_below_alpha", "lower_above_alpha", "max_reached"))
  expect_lte(res$n_perms_used, 999L)
})


# --- DR-DATE adaptive smoke test --------------------------------------

test_that("dr_date_test with adaptive = TRUE returns adaptive fields", {
  skip_on_cran()
  set.seed(9)
  n <- 200
  x <- matrix(rnorm(n * 2), n, 2)
  t <- rbinom(n, 1, plogis(0.5 * x[, 1]))
  y <- t * 1.2 + x[, 1] + rnorm(n, sd = 0.5)

  res <- dr_date_test(y, t, x,
    adaptive = TRUE, B_max = 999L, batch_size = 100L,
    alpha = 0.05, seed = 9, outcome_model = "zero"
  )
  expect_s3_class(res, "kernel_test_result")
  expect_true(!is.null(res$n_perms_used))
  expect_true(res$stop_reason %in% c("upper_below_alpha", "lower_above_alpha", "max_reached"))
})
