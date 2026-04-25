# Contributing to kernR

Thanks for your interest in contributing. kernR welcomes bug reports, ideas, and pull requests.

## Quick links

- **Bug reports / feature requests**: open an [issue](https://github.com/AAGI-AUS/kernR/issues) using the appropriate template.
- **Security disclosures**: see [SECURITY.md](SECURITY.md). Do *not* file a public issue for security findings.
- **Code of Conduct**: this project follows the [Contributor Covenant 2.1](CODE_OF_CONDUCT.md).

## Development setup

```r
# 1. Fork and clone the repo, then:
pak::pak("local::.", dependencies = TRUE)

# 2. Iterative development (do NOT R CMD INSTALL while the package is loaded —
#    on Box-cloud-synced paths it can corrupt the lazy-load DB):
devtools::load_all()
devtools::test()

# 3. Lint and document before pushing:
lintr::lint_package()
devtools::document()

# 4. Full check on a tarball (CRAN form). On macOS with a broken FLIBS toolchain
#    you may need R_MAKEVARS_USER=/dev/null and a temporary strip of $(FLIBS)
#    from src/Makevars; restore it before committing.
R CMD build .
R CMD check --as-cran kernR_*.tar.gz
```

## Style conventions

- **R**: snake_case throughout. `data.table` is the default tabular backend; tidyverse is permitted where it improves clarity. Australian/British English.
- **C++**: RcppArmadillo. Keep changes localised; avoid global state.
- **Tests**: `testthat` edition 3, parallel-safe. Add at least one test for every new exported function and every bug fix.
- **Documentation**: roxygen2 with `@returns`, `@examples` (executable, no blanket `\dontrun{}`), `@family` for grouping. `air format` is encouraged but not yet wired into pre-commit.
- **Lint**: `lintr::lint_package()` must be clean before merge. Per-location `# nolint` comments are acceptable for justified cases (with the linter name and a one-line rationale).

## Pull request flow

1. One PR per logical change. Tests must pass on CI (Ubuntu / macOS / Windows × R release / devel / oldrel).
2. Update `NEWS.md` under the dev-version heading for any user-visible change.
3. If you change function signatures, regenerate Rd files (`devtools::document()`) and include them in the PR.
4. The repository uses the `main` branch as the default; please rebase on `main` before requesting review.

## Maintainer competence matrix

| Layer | Primary | Secondary |
|---|---|---|
| R | Max Moldovan | Dino Sejdinovic |
| C++ / RcppArmadillo | Max Moldovan | (seeking second reviewer) |
| Statistical methodology | Dino Sejdinovic | Max Moldovan |

## Reporting research bugs

Numerical correctness is taken seriously. If you find a discrepancy between kernR's output and a reference implementation (the original Python `kgformula` / `DR_distributional_test` repos, or another R kernel-test package), please attach a minimal reproducible example and the comparison output.
