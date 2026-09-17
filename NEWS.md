# stagelife 0.1.0

Initial release.

**stagelife** implements the stage life testing (SLT) model of Laumen &
Cramer (2019), "Stage life testing", *Naval Research Logistics*
66(8):632-647, <doi:10.1002/nav.21874> -- a two-stage cumulative-exposure
life-testing model in which units begin on a baseline stress and, at a
fixed change-point, some or all surviving units move to an elevated
stress.

* Twelve baseline lifetime distributions (Weibull, Power Lindley,
  Exponentiated Lindley, Exponentiated Exponential, Marshall-Olkin
  Exponential, Lomax, Exponentiated Teissier, Log-Logistic, Chen,
  Log-normal, Gamma, and Power Muth), each with closed-form density,
  distribution, and quantile functions.
* Three parameter-sharing cases (`common`, `shape1`, `shape11`) for
  the elevated-stage shape parameter.
* Data generation under the random-withdrawal design, with both
  fixed-proportion and fixed-count selection rules, via `rslt()`; the
  classical simple-step-stress design is the `prop_transfer = 1`
  special case.
* Maximum likelihood estimation via `fit_slt()`, with standard
  `print()`, `summary()`, `coef()`, `vcov()`, `logLik()`, `AIC()`, and
  `BIC()` methods.
* Asymptotic (Wald) and bootstrap (non-parametric and parametric;
  percentile, bootstrap-*t*, and log-scale bootstrap-*t*) confidence
  intervals.
* Design-planning tools (`slt_prob_stage1()`), Kolmogorov-Smirnov
  goodness-of-fit testing with a parametric-bootstrap correction
  (`slt_gof_ks()`), likelihood-ratio tests between nested cases
  (`slt_lrt()`), and AIC/BIC model comparison across fits
  (`slt_model_comparison()`).
* Diagnostic plots: fitted CDF, survival, P-P, and Q-Q.
* Four bundled real step-stress data sets.

See `vignette("stagelife-intro")` for a full worked example.
