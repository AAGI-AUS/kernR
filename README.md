# kernR

<!-- badges: start -->
[![R-CMD-check](https://github.com/AAGI-AUS/kernR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/AAGI-AUS/kernR/actions/workflows/R-CMD-check.yaml)
[![test-coverage](https://github.com/AAGI-AUS/kernR/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/AAGI-AUS/kernR/actions/workflows/test-coverage.yaml)
[![codecov](https://codecov.io/gh/AAGI-AUS/kernR/graph/badge.svg)](https://app.codecov.io/gh/AAGI-AUS/kernR)
[![lint](https://github.com/AAGI-AUS/kernR/actions/workflows/lint.yaml/badge.svg)](https://github.com/AAGI-AUS/kernR/actions/workflows/lint.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

**Kernel-Based Causal Distributional Testing for R**

kernR implements kernel-based hypothesis tests for causal inference and distributional treatment effects. It detects effects that standard methods (t-tests, linear regression, Double ML) miss entirely — including non-linear causal associations, variance shifts, and distributional changes that leave the mean unchanged.

## Key Features

- **HSIC independence test** — detects non-linear dependencies that Pearson correlation misses
- **MMD two-sample test** — detects any distributional difference, not just mean shifts
- **Backdoor-adjusted HSIC (bd-HSIC)** — tests causal association after adjusting for confounders via density ratio estimation [(Hu, Sejdinovic & Evans, 2024)](https://jmlr.org/beta/papers/v25/21-1409.html)
- **Doubly robust distributional tests (DR-DATE, DR-DETT)** — tests for distributional treatment effects with double robustness [(Fawkes, Hu, Evans & Sejdinovic, 2024)](https://openreview.net/pdf?id=5g5zFVj33K)
- **Hierarchical/nested data** — within-cluster and between-cluster decomposition for clustered data (farms, hospitals, schools)
- **Formula interface** — `kernel_causal_test(y ~ treatment | confounders, data = df)`
- **Fast C++ backend** — RcppArmadillo kernel engine with permutation inference

## Installation

```r
# Install from GitHub (recommended via pak)
# install.packages("pak")
pak::pak("AAGI-AUS/kernR")
```

## Quick Start

```r
library(kernR)

# Detect non-linear dependence (HSIC)
set.seed(42)
x <- rnorm(200)
y <- x^2 + rnorm(200, sd = 0.3)  # quadratic relationship
cor.test(x, y)$p.value            # Pearson: p = 0.85 (misses it!)
hsic_test(x, y, seed = 1)        # HSIC:    p = 0.002 (detects it)

# Test for distributional treatment effect (DR-DATE)
n <- 300
covariates <- matrix(rnorm(n * 2), n, 2)
treatment  <- rbinom(n, 1, plogis(0.3 * covariates[, 1]))
outcome    <- treatment * rnorm(n, sd = 2) +       # treatment changes variance
              (1 - treatment) * rnorm(n, sd = 1) +  # but NOT the mean
              covariates[, 1]

t.test(outcome[treatment == 1], outcome[treatment == 0])  # p = 0.72 (blind)
dr_date_test(outcome, treatment, covariates, seed = 1)    # p = 0.002 (sees it)
```

## When to Use kernR

| Your question | Traditional method | kernR function | Why kernR wins |
|---------------|-------------------|----------------|----------------|
| Are X and Y related? (possibly non-linearly) | `cor.test()` | `hsic_test()` | Detects U-shapes, periodicity, any non-linear pattern |
| Do two groups differ? (beyond means) | `t.test()` | `mmd_test()` | Detects variance, shape, and tail differences |
| Does X causally affect Y? (with confounders) | Double ML | `bd_hsic_test()` | Non-parametric; no functional form assumptions |
| Does treatment change the outcome distribution? | `t.test()` / TMLE | `dr_date_test()` | Full distributional comparison; doubly robust |
| Effect on treated subgroup only? | — | `dr_dett_test()` | One-sided overlap; robust in imperfect settings |
| Clustered/hierarchical data? | Mixed models | `hierarchical_test()` | Non-parametric + proper permutation within clusters |

## Methodological Foundation

kernR implements methods from two peer-reviewed papers:

1. **Hu, R., Sejdinovic, D., & Evans, R. J.** (2024). A kernel test for causal association via noise contrastive backdoor adjustment. *Journal of Machine Learning Research*, 25(160), 1–56. [Paper](https://jmlr.org/beta/papers/v25/21-1409.html) | [Original code](https://github.com/MrHuff/kgformula)

2. **Fawkes, J., Hu, R., Evans, R. J., & Sejdinovic, D.** (2024). Doubly robust kernel statistics for testing distributional treatment effects. *Transactions on Machine Learning Research*. [Paper](https://openreview.net/pdf?id=5g5zFVj33K) | [Original code](https://github.com/Jakefawkes/DR_distributional_test)

The hierarchical extension for nested/clustered data is a novel contribution of kernR.

## Python Companion

A Python translation is available as [**kernP**](https://github.com/AAGI-AUS/kernP) with an equivalent API built on NumPy/SciPy/scikit-learn.

## Vignettes

- `vignette("kernR-quickstart")` — 10-minute introduction
- `vignette("kernR-bdhsic")` — Causal association testing
- `vignette("kernR-drtest")` — Distributional treatment effect tests
- `vignette("kernR-hierarchical")` — Hierarchical/nested data

## Licence

MIT
