# kernR (development version)

## kernR 0.1.0

### New features

* **Low-rank kernel approximations** as opt-in attributes of `kernel_spec()`:
  - `approx = "nystrom"` -- thresholded eigendecomposition of the
    landmark Gram block. Cheap factor exposed as `nystrom_factor()`.
  - `approx = "rff"` -- Random Fourier Features for translation-invariant
    kernels (RBF, Matern via Student-t spectral sampling). Cheap feature
    map exposed as `rff_features()`.
  - Both honoured by `kernel_matrix()` for symmetric Gram matrices;
    rectangular calls fall back to exact computation.
  - `approx_rank` and `approx_seed` arguments give callers explicit
    control; default rank is `max(50, ceiling(0.1 * n))` (capped at n).

* **Adaptive permutation inference** via the new
  `adaptive = TRUE`/`B_max`/`batch_size` arguments on `mmd_test()`,
  `hsic_test()`, `bd_hsic_test()`, `dr_date_test()`, and
  `dr_dett_test()`. Sequential procedure with Wilson-interval stopping
  (Besag & Clifford 1991; Gandy 2009); exits early when the running
  p-value's confidence interval clearly excludes `alpha`. Default
  `adaptive = FALSE` preserves byte-identical behaviour.

* **Two new vignettes**:
  - `kernR-theory.Rmd` -- canonical mathematical reference
    (RKHS, MMD, HSIC, bd-HSIC, conditional mean embeddings,
    DR-DATE/DETT, hierarchical extension).
  - `kernR-performance.Rmd` -- benchmarks and tuning guidance for the
    approximation backends.

* **pkgdown site** scaffolding (`_pkgdown.yml`) with grouped reference
  index. Build locally with `pkgdown::build_site()`.

### Infrastructure

* GitHub Actions workflows added: `R-CMD-check`, `test-coverage`,
  `lint`.
* `.lintr` configuration; package is lint-clean.

## kernR 0.0.0.9000

* Initial development version.
* **Kernel engine**: RBF, Matern, linear, polynomial kernels with RcppArmadillo backend.
* **Bandwidth selection**: Median heuristic (Rcpp), Scott's rule.
* **Base tests**: `hsic_test()` (independence), `mmd_test()` (two-sample).
* **Causal tests**:
  - `bd_hsic_test()`: Backdoor-adjusted HSIC for causal association testing.
  - `dr_date_test()`: Doubly robust distributional average treatment effect.
  - `dr_dett_test()`: Doubly robust distributional effect on the treated.
* **Hierarchical**: `hierarchical_test()` for nested/clustered data with within/between decomposition.
* **Unified interface**: `kernel_causal_test()` with formula syntax `y ~ treatment | confounders`.
* **Density ratio estimation**: Logistic NCE via `glm`/`ranger`/`xgboost`; RuLSIF (kernel-based).
* **Propensity scores**: Cross-fitted estimation with trimming and diagnostics.
* **Diagnostics**: `assess_overlap()`, `plot_weights()`, `effective_sample_size()`.
* **Vignettes**: Quick start, bd-HSIC tutorial, DR-DATE/DETT tutorial, hierarchical data.
