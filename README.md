# stagelife

<!-- badges: start -->
<!-- badges: end -->

**stagelife** implements the **stage life testing (SLT)** model of
[Laumen & Cramer (2019)](https://doi.org/10.1002/nav.21874), *Naval
Research Logistics* 66:632-647 -- a two-stage life-testing model that
unifies progressive censoring with fixed censoring times and step-stress
accelerated life testing under the cumulative-exposure model -- together
with a suite of simulation, estimation and validation tools.

* **12 baseline lifetime distributions**: Weibull, Power Lindley, Exponentiated
  Lindley, Exponentiated Exponential, Marshall-Olkin Exponential, Lomax,
  Exponentiated Teissier, Log-Logistic (Fisk), Chen, Log-normal, Gamma,
  and Power Muth (Power Teissier). Only Weibull, Exponentiated
  Exponential, Marshall-Olkin Exponential, and Gamma reduce to
  Exponential at `alpha = 1` (see `vignette("stagelife-intro")` for the
  full table); `gamma` and `lognormal` are built directly on
  `stats::pgamma()`/`stats::pnorm()` rather than a hand-derived formula.
* **3 identifiable parameter-sharing cases**: `shape1` and `shape11` are
  the paper's own models (Weibull-Exponential, Section 4.2, and
  Exponential-Exponential, Section 4.1); `common` (shared shape) is this
  package's extension. A fourth, fully unrestricted case is deliberately
  not offered -- see `vignette("stagelife-intro")` for why.
* **One data-generating function**: `rslt()` implements the paper's own
  random-withdrawal SLT design (Procedure 2.3), with both the
  `"typeP"`/fixed-proportion and `"typeM"`/fixed-count selection rules;
  the classical simple-step-stress design (every survivor transferred)
  is simply its `prop_transfer = 1` special case (Remark 2.2(4)), so no
  separate function is needed for it.
* **Estimation**: maximum likelihood via a Nelder-Mead/BFGS hybrid with
  multi-start (`fit_slt()`), plus standard `coef()`/`vcov()`/`logLik()`/
  `AIC()`/`BIC()`/`confint()` methods. Verified to reproduce the paper's
  own worked numerical example (Table 2) exactly.
* **Inference**: asymptotic Wald (`normal`/`lognormal`, the latter
  guaranteeing a positive interval) and bootstrap (non-parametric and
  parametric; percentile, bootstrap-*t*, and log-scale bootstrap-*t*,
  the last likewise guaranteeing positivity) confidence intervals
  (`slt_bootstrap_ci()`).
* **Experimental design**: `slt_prob_stage1()` implements the paper's
  Lemma 4.1 (probability of observing at least one elevated-stage
  failure), exactly reproducing the paper's own Table 5.
* **Goodness of fit**: Kolmogorov-Smirnov test with a parametric-bootstrap
  correction for estimated parameters (`slt_gof_ks()`), likelihood-ratio
  tests between nested cases (`slt_lrt()`), and an AIC/BIC comparison
  table across any set of fits (`slt_model_comparison()`).
* **Diagnostics**: fitted CDF/survival, P-P and Q-Q plots.
* **4 bundled real step-stress data sets** for examples.

**Note**: an earlier version of this package included
`slt_pivot_ci_theta1()`, a design-based confidence interval claimed to
generalise across all baseline distributions. It has been **removed**: it
silently assumed the elevated stage becomes Exponential whenever `case`
fixes its shape to 1, which is only true for `weibull`/`expexponential`/
`moexponential`/`gamma` and gives a meaningfully wrong interval otherwise (confirmed by
direct 200,000-replicate simulation). `slt_bootstrap_ci()` is correct
for every family/case combination and is the supported alternative.

**Note**: the paper's own *exact* confidence interval for the
baseline-stage parameter (its Corollary 4.4(1)) is not implemented here,
since it depends on a theorem proved in a separate, inaccessible
companion paper; see `vignette("stagelife-intro")` for a full
explanation and the recommended alternative.

**On simulation studies**: this package deliberately does not include a
bias/coverage/RMSE simulation-study engine. `rslt()`, `fit_slt()`, and
`confint()`/`slt_bootstrap_ci()` are the complete set of building blocks
needed to run one yourself in a plain loop over your own scenario grid.

## Installation

```r
# From the built tarball:
# install.packages("remotes")
remotes::install_local("stagelife_0.1.0.tar.gz")

# Or directly from GitHub:
remotes::install_github("vksharma-bhu/Stagelife")
```

## Quick start

```r
library(stagelife)

theta0 <- c(0.05, 1.5); theta1 <- c(0.20, 1)
tau1 <- slt_tau("weibull", theta0, p = 0.5)
d <- rslt(200, "weibull", theta0, theta1, tau1, prop_transfer = 1)

fit <- fit_slt(d, family = "weibull", case = "shape1", tau1 = tau1)
summary(fit)
confint(fit)
plot_cdf(fit)
```

See `vignette("stagelife-intro", package = "stagelife")` for the full
workflow, including goodness of fit, model comparison, and a real-data
example.
