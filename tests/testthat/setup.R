## Test-suite setup. Sourced once before any test_that() block.

# Predictable locale and time zone so format()/print() snapshots stay stable.
withr::local_locale(c(LC_COLLATE = "C", LC_CTYPE = "en_AU.UTF-8"),
                    .local_envir = teardown_env())
withr::local_envvar(TZ = "UTC", .local_envir = teardown_env())

# Predictable RNG kind. testthat 3 sets a per-test seed; this only controls
# the underlying sampler, not the seed.
RNGversion(getRversion())

# testthat parallel safety: avoid mucking with the user's options() global.
# Keep cli outputs deterministic in snapshots.
withr::local_options(
  cli.unicode = FALSE,
  cli.dynamic = FALSE,
  cli.num_colors = 1L,
  width = 80L,
  .local_envir = teardown_env()
)
