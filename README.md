# stagelife

<!-- badges: start -->
<!-- badges: end -->

**stagelife** implements the **stage life testing (SLT)** model of
[Laumen & Cramer (2019)](https://doi.org/10.1002/nav.21874), *Naval
Research Logistics* 66:632-647 -- a two-stage life-testing model that
unifies progressive censoring with fixed censoring times and step-stress
accelerated life testing under the cumulative-exposure model -- together
with a suite of simulation, estimation, and validation tools.

## Features

* **12 baseline lifetime distributions**: Weibull, Power Lindley,
  Exponentiated Lindley, Exponentiated Exponential, Marshall-Olkin
  Exponential, Lomax, Exponentiated Teissier, Log-Logistic, Chen,
  Log-normal, Gamma, and Power Muth.
* **3 parameter-sharing cases** for the elevated-stage shape parameter:
  `shape1` and `shape11` (the paper's own models) and `common` (a
  shared-shape extension).
* **Data generation** under the random-withdrawal design (`rslt()`),
  with both fixed-proportion and fixed-count selection rules; the
  classical simple-step-stress design is available as a special case.
* **Estimation**: maximum likelihood via `fit_slt()`, with standard
  `coef()`, `vcov()`, `logLik()`, `AIC()`, `BIC()`, and `confint()`
  methods.
* **Inference**: asymptotic (Wald) and bootstrap (non-parametric and
  parametric) confidence intervals, via `slt_bootstrap_ci()`.
* **Experimental design**: `slt_prob_stage1()` computes the probability
  of observing at least one elevated-stage failure, for planning a test.
* **Goodness of fit and model comparison**: a Kolmogorov-Smirnov test
  with a parametric-bootstrap correction (`slt_gof_ks()`), likelihood-
  ratio tests between nested cases (`slt_lrt()`), and an AIC/BIC
  comparison table across fits (`slt_model_comparison()`).
* **Diagnostics**: fitted CDF, survival, P-P, and Q-Q plots.
* **4 bundled real step-stress data sets** for examples.

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
