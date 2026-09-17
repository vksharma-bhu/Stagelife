# stagelife 0.1.0

## Data source corrections

* **Two of the four bundled data sets had incorrect or missing
  citations, now corrected against the original source papers.**
  `slt_kateri_nikolov_2024` (the data values were always correct, but
  the citation was wrong -- a real paper by a partially-overlapping
  author list, not the actual source) is renamed
  **`slt_bobotas_kateri_2015`**: Bobotas, P. and Kateri, M. (2015). The
  step-stress tampered failure rate model under interval monitoring.
  *Statistical Methodology*, 27, 100-122. `slt_data3` (previously
  shipped with no identified source at all) is renamed
  **`slt_wang_fei_2003`**: Wang, R. and Fei, H. (2003). Uniqueness of
  the maximum likelihood estimate of the Weibull distribution tampered
  failure rate model. *Communications in Statistics - Theory and
  Methods*, 32(12), 2321-2338.
* **`slt_zhu_2010`'s change-point documented as a discrepancy, not
  silently resolved either way.** The source describes the nominal
  change-point as 96 hours, but this is inconsistent with the recorded
  data itself: the smallest stage-1 observation is 94.38, which would
  be impossible if the change-point were 96 (a unit cannot fail on
  stage 1 before the stage change happens). `tau1 = 94` -- the only
  value consistent with both stages of the recorded data -- was already
  in use and is kept; the discrepancy is now explicitly documented in
  `?slt_zhu_2010` rather than left unremarked.

## Export surface

* **`slt_expand_params()` and `neg_loglik_slt()` are now internal**
  (`stagelife:::`, not `stagelife::`). Neither has a genuine use case
  outside the package's own fitting code: `fit_slt()`'s return object
  already carries `theta0`/`theta1` directly, so a user never needs to
  expand a raw parameter vector themselves, and `neg_loglik_slt()` is
  purely the `optim()` objective. `slt_cases()` and `slt_reduce_params()`
  remain exported -- the former is a lookup table for the `case`
  argument (the same role `slt_families()` plays for `family`), and the
  latter has a real external use (matching a simulation study's known
  "true" parameters to `fit_slt()`'s output format).
* **`nobs.slt_fit()` removed entirely** (not merely hidden). Checked
  directly against `stats:::BIC.default`'s source before removing it:
  `BIC()`/`AIC()` read the sample size from `attr(logLik(object),
  "nobs")` first, which `logLik.slt_fit()` already sets directly, and
  only fall back to calling `nobs(object)` if that attribute is
  missing -- so `nobs.slt_fit()` was not actually being used by
  anything in the package. Calling `nobs(fit)` directly now falls
  through to `stats::nobs.default`, which gives a clear error ("no
  'nobs' method is available") rather than a silent wrong answer;
  `fit$n` remains the documented, direct way to get the sample size.

## Family roster

* **Gompertz removed; Marshall-Olkin Exponential and Exponentiated
  Lindley added** (net: still growing the catalogue, now 12 families).
  Gompertz's shape parameter multiplies the failure time directly
  inside an exponential, which was verified (by testing whether
  rescaling `x` by a constant `c` can be absorbed into `lambda` alone)
  to make it behave as a second scale parameter rather than a portable,
  unit-invariant shape descriptor -- confirmed both by direct
  dimensional analysis (a joint `lambda`/`alpha` rescaling by the same
  factor `c` reproduces the rescaled distribution exactly; rescaling
  `lambda` alone does not) and by the empirical convergence failures
  this caused at some data scales (see the `.default_init()` entry
  below). **Marshall-Olkin Exponential** (Marshall & Olkin, 1997) adds
  a shape/tilt parameter to the exponential baseline by a different
  mechanism than the exponentiated/power-transform families already in
  the package, and reduces to Exponential exactly at `alpha = 1`.
  **Exponentiated Lindley** (Nadarajah, Bakouch & Tahmasbi, 2011,
  published as the "generalized Lindley distribution") is the
  "exponentiate the cdf" sibling of the existing Power Lindley (exactly
  as Exponentiated Teissier and Power Muth are siblings of the Teissier
  baseline); its quantile function reuses the same Lambert-*W*
  baseline-Lindley inverter Power Lindley already uses. Chen was kept
  despite sharing Gompertz's lack of a portable shape parameter --
  confirmed by the same rescaling test, a fine grid search found no
  value of `lambda` alone that reproduces a rescaled Chen distribution,
  and unlike Gompertz, Chen does not even have a joint `lambda`/`alpha`
  escape hatch that works instead. `alpha = 1` is still a legitimate
  special case of the Chen family in
  its own right (Chen's shape appears as an ordinary exponent on `x`,
  the same structural role a shape parameter plays in Weibull/Lomax/
  etc.), even though it is not a unit-invariant one.

## Correctness

* **Exponentiated Teissier and Power Muth now use their exact
  closed-form (Lambert *W*) quantile function, not a numerical
  inversion.** An earlier version of this package claimed these two
  families have no closed-form quantile at all, "not even via Lambert
  *W*", and inverted their cdfs numerically via `stats::uniroot()`
  instead. This was simply wrong: both have an exact closed form
  (Sharma, Singh & Shekhawat, 2022, Proposition 4.2, for Exponentiated
  Teissier -- the source paper for this family, brought to the
  maintainers' attention directly -- and the analogous identity for
  Power Muth, following Jodra, Jimenez-Gamero & Alba-Fernandez, 2015).
  The error traced back to an incorrect general claim that equations of
  the form `u - e^u = -K` (these two families' shared cumulative-hazard
  equation) cannot be solved via Lambert *W*; verified directly that
  they can, via a different substitution than the more familiar `w e^w
  = z` product form. The closed form was checked against the
  previous numerical quantile (agreement to ~1e-9, i.e. `uniroot()`'s
  own tolerance) before replacing it, and improves round-trip accuracy
  from ~1e-8 to ~1e-15 (machine precision) while being substantially
  faster (no iterative root-finding). The now-unused generic numerical
  quantile inverter (`num_quantile_generic()`) has been removed.

## Documentation

* **Markdown mode enabled** (`Roxygen: list(markdown = TRUE)` in
  `DESCRIPTION`). This had never been set, so every `[func()]`-style
  cross-reference written anywhere in the package's documentation was
  rendering as literal bracket text rather than an actual hyperlink;
  enabling it retroactively repaired links throughout the entire
  manual (verified by re-rendering the full PDF manual and inspecting
  pages written well before this fix).
* **Full mathematical reference for all twelve baseline distributions**
  added to `slt_dist()`: pdf/cdf in proper notation, literature
  citation for each, and a "Numerical notes" section explaining
  precisely which families need the Lambert *W* function or numerical
  root-finding for their quantile and why (verified numerically, not
  just asserted).
* **The two-stage population cdf/pdf, the cumulative-exposure-model
  construction, and their exact relationship to the estimation
  log-likelihood** are now derived with full notation in
  `slt_equiv_age()` and `pslt()`/`dslt()`/`qslt()`, including when the
  two agree exactly and when they don't (partial-withdrawal designs).
* **The three parameter-sharing cases are now derived mathematically**
  in `slt_cases()` (the parameter-index mapping each case implements,
  why the fully unrestricted case is excluded, and the nesting
  relationships behind `slt_lrt()`), backed by a fresh end-to-end
  re-verification that all three cases recover correct,
  correctly-labelled parameter estimates.
* **Five previously-undocumented S3 methods**
  (`print.slt_fit`/`coef.slt_fit`/`vcov.slt_fit`/`logLik.slt_fit`/
  `print.summary.slt_fit`) now have their own proper
  documentation topics.
* **All four bundled data sets renamed** to a uniform
  `slt_<author>_<year>` convention, with researched (not assumed)
  citations: `han_kundu_2014` -> `slt_han_kundu_2014`, `zhu_2010` ->
  `slt_zhu_2010`. The other two were subsequently corrected again once
  the actual source paper became available (see the "Data source
  corrections" entry above): `kateri_nikolov` -> `slt_bobotas_kateri_2015`
  and `slt_dataset3` -> `slt_wang_fei_2003`.
* Title and Description revised for accuracy; `Imports`/`Suggests`
  alphabetised.

## Initial release

Implements the **stage life testing (SLT)** model of
Laumen & Cramer (2019), "Stage life testing", *Naval Research Logistics*
66:632-647: units begin on a baseline stage and, at a fixed stage-change
time `tau1`, `R*_1` of the still-surviving units are randomly withdrawn
and moved to an elevated stage, linked by the cumulative-exposure model.

## Model and data generation

* `rslt()` is the package's single data-generating function,
  implementing the paper's own Procedure 2.3, with both of the paper's
  selection rules: `design = "typeP"` (fixed target proportion `pi1`)
  and `design = "typeM"` (fixed target count `R0_1`). The classical
  "simple step-stress" design (every survivor transferred) is simply the
  `prop_transfer = 1` special case (the paper's Remark 2.2(4)) -- no
  separate function is needed for it.
* Twelve baseline lifetime distributions: Weibull, Power Lindley, Exponentiated
  Lindley, Exponentiated Exponential, Marshall-Olkin Exponential, Lomax,
  Exponentiated Teissier, Log-Logistic (Fisk), Chen, Log-normal, Gamma,
  and Power Muth (Power Teissier; Jodra, Gomez, Jimenez-Gamero &
  Alba-Fernandez, 2017). `weibull`, `expexponential`, `lomax`,
  `expteissier`, and `powerlindley` are the paper's own baseline-stage
  extension points; the other seven are this package's additions,
  chosen to cover hazard shapes the original five do not (a
  non-monotonic hazard for `loglogistic`; a bathtub-shaped hazard for
  `chen`; two of the most standard lifetime distributions, `gamma` and
  `lognormal`, built directly on `stats::pgamma()`/`stats::pnorm()`
  rather than a hand-derived formula to avoid an entire class of
  algebra/overflow bugs; `powermuth`/`explindley`, the "power the
  argument" siblings of `expteissier`/`powerlindley`'s "exponentiate
  the cdf", each pair generalising the same baseline distribution --
  verified to agree with its sibling exactly at `alpha = 1`;
  `moexponential`, Marshall & Olkin's (1997) "tilt parameter"
  construction applied to the exponential baseline). Only `weibull`,
  `expexponential`, `moexponential`, and `gamma` reduce to Exponential
  at `alpha = 1`; see `vignette("stagelife-intro")` for the full table
  and why this distinction matters.
* Three identifiable parameter-sharing cases (`common`, `shape1`,
  `shape11`); the fully unrestricted ("full") case is intentionally not
  supported due to empirical non-identifiability. `shape1` and `shape11`
  are the paper's own models (Sections 4.2 and 4.1 respectively);
  `common` is this package's extension.

## Estimation and inference

* Maximum likelihood estimation (`fit_slt()`) with a Nelder-Mead/BFGS
  hybrid optimiser and multi-start.
* Asymptotic Wald confidence intervals (`confint()`): `type =
  "lognormal"` (default, always positive) and `type = "normal"`.
* Bootstrap confidence intervals (`slt_bootstrap_ci()`, non-parametric
  or parametric): percentile (always positive), bootstrap-*t*, and
  log-scale bootstrap-*t* (`log_boot_t`, always positive by the same
  log-and-exponentiate construction as `type = "lognormal"` above) --
  the plain bootstrap-*t* interval, like `type = "normal"`, is an
  ordinary linear pivot and can occasionally produce a negative lower
  bound for a positive parameter at small `n`.
* Parametric-bootstrap-corrected Kolmogorov-Smirnov goodness-of-fit test
  (`slt_gof_ks()`), likelihood-ratio tests (`slt_lrt()`), and AIC/BIC
  model comparison (`slt_model_comparison()`).
* `slt_prob_stage1()` implements the paper's Lemma 4.1 (probability of
  observing at least one elevated-stage failure), for experimental
  design/planning.
* Fitted CDF/survival/P-P/Q-Q diagnostic plots.
* Four bundled real step-stress data sets.

## API naming

Every exported function name is consistent with the package being about
exactly one model (the SLT model of the paper): the population
density/CDF/quantile functions follow R's usual `d`/`p`/`q`/`r` convention
(`dslt()`, `pslt()`, `qslt()`, `rslt()`), estimation/inference/GOF
functions are prefixed `slt_*` or take an `"slt_fit"` object, and the
single stage-change-time argument is named `tau1` everywhere (matching
the paper's own notation, including the field name in the four bundled
real data sets).

## Fidelity to Laumen & Cramer (2019)

Verified by directly reproducing the paper's own numerical results (see
`vignette("stagelife-intro")`, Section 5):

* `fit_slt()` reproduces the paper's Section 5 worked example (Table 2,
  Type-P design) MLEs exactly (theta0_hat = 32.35, theta1_hat = 27.55).
* `slt_prob_stage1()` exactly reproduces all four rows of the paper's
  Table 5, including a correction of what appears to be an
  OCR/typesetting rendering issue in the printed equation (4.1): the
  verified formula is `n* = floor(n - 1/pi1)` for the `"typeP"` design,
  not `floor((n-1)/pi1)` as the available scan renders it.
* The paper's *exact* confidence interval for the baseline-stage
  parameter (its Corollary 4.4(1)) is **not** implemented: it depends on
  a finite-sample distribution theorem proved in a separate, inaccessible
  companion paper (Laumen & Cramer, 2019, "Progressive censoring with
  fixed censoring times", *Statistics*, 53, 569-600). A best-effort
  transcription of the printed formula was tried and failed Monte Carlo
  validation, so it was not shipped. `slt_bootstrap_ci(method =
  "parametric")` with `prop_transfer` set to your actual design is the
  recommended alternative.
* A design-based confidence interval for the elevated-stage parameter
  (a claimed generalisation of the paper's Theorem 4.2(2)/4.7) was
  developed and shipped in an earlier development version, then
  **removed**: it silently relied on the elevated stage becoming
  Exponential whenever `case` fixes its shape to 1, which is only true
  for the three families whose `alpha = 1` case actually is Exponential
  (`weibull`, `expexponential`, `gamma`) and gave a meaningfully wrong
  interval for the others (confirmed by direct 200,000-replicate
  simulation, differing from the true survival function by 0.08-0.30 in
  probability). `slt_bootstrap_ci()` is correct for every family/case
  combination and is the supported alternative.
* `slt_gof_ks()` and the four diagnostic plot functions assume the
  simple-step-stress special case (`prop_transfer = 1`): their pooled
  reference distribution (`pslt()`) is only the correct marginal CDF of
  the combined `Y`/`Z` sample in that case, since a non-transferred
  survivor under genuine partial withdrawal keeps failing on stage 0
  unmodified, which `pslt()` does not represent. A fully general version
  (separate `Y`-vs-`F0` and `Z`-vs-shifted-`F1` tests, both of which are
  design-invariant) would be a reasonable future extension but is not
  currently implemented.

## Numerical robustness

* **Power Lindley's quantile function** is hardened against extreme
  `lambda`: the Lambert-*W*-based closed form is numerically unsound
  (`NaN`, `Inf`, or silently wrong-signed) once `1 + lambda` approaches
  the point where `exp(-(1+lambda))` underflows (around `lambda > 745`),
  and loses precision well before that. A large-`lambda` asymptotic
  branch takes over once `1 + lambda > 400`, and `lambda = Inf`/`NaN`
  probability inputs are handled explicitly.
* **The Teissier/Power Muth shared baseline** (`.teiss_G0`/`.teiss_g0`,
  used by `expteissier` and `powermuth`) is hardened against three
  nested layers of floating-point cancellation for small `lambda*x`:
  the naive `exp(z)-1-z` is wrong (even wrong-*signed*) by `z=1e-8`;
  `expm1(z)-z` fixes that down to about `z=1e-8` but itself fails below
  that, since the `z^2/2` correction term becomes smaller than the
  ~52-bit relative precision with which `z` itself is representable (a
  hardware/IEEE-754 floor, not an implementation bug) -- fixed with a
  direct Taylor-series fallback for `z < 1e-8`. A second, independent
  instance of the same pattern (`1-exp(-H0)` for tiny `H0`) produced
  spurious `Inf` when combined with a shape `< 1`; fixed with
  `-expm1(-H0)`. The same `exp(z)-1`-type cancellation risk in
  `chen` is fixed defensively the same way.
* **A single numerically difficult replicate no longer aborts an entire
  bootstrap or goodness-of-fit run.** `slt_bootstrap_ci()` and
  `slt_gof_ks()` wrap every individual replicate in its own error
  handler and discard failed replicates (including
  `parallel::mclapply()`'s `"try-error"` objects, which would otherwise
  slip past a `NULL`-only filter and corrupt downstream aggregation).
* **`slt_bootstrap_ci()` does not require more successful replicates
  than were requested.** The minimum-convergence check is capped at
  `min(B, 30)` rather than a hardcoded floor of 30.
* **`fit_slt()`/`.fit_slt_once()` never report `converged = TRUE` while
  stuck at the internal optimiser penalty sentinel.** `optim()` can
  report `convergence = 0` ("successful") while every point it visited
  during the search returned the same penalty value (confirmed for
  `family = "gompertz"` on real data, where a poorly-scaled starting
  region made a shifted elevated-stage observation numerically `Inf`
  everywhere the optimiser looked); a fit whose log-likelihood is
  anywhere near that sentinel is now always reported as
  non-convergence rather than a false success.

## Input validation

* **Every family rejects an invalid `theta`** (`NA`/`NaN`, negative, or
  zero `lambda`/`alpha`) as `NA` instead of a silently wrong but
  numerically plausible value (verified case-by-case: e.g. `weibull`
  previously returned a *negative* CDF for `lambda < 0`, `expexponential`
  a CDF *greater than 1* for `alpha < 0`). `lambda = +Inf` (a legitimate
  degenerate boundary for `powerlindley`) is deliberately still accepted.
* **`confint()` and `slt_bootstrap_ci()` validate `level`.** `level`
  outside `(0, 1)` previously produced either a cryptic warning or, for
  some values, a silently backwards interval (lower bound greater than
  upper bound) with no warning at all.
* **`fit_slt()` rejects non-positive/non-finite observed failure times**
  (`data$Y`/`data$Z`), and requires at least one observation on *each*
  stage: with `Z` empty, for example, the stage-1 log-likelihood is
  identically 0 for every value of `lambda1` (a perfectly flat
  objective), so `optim()` previously returned the untouched starting
  value labelled as a converged MLE.

## Scope note

This package deliberately does not include a simulation-study engine
(bias/coverage/RMSE evaluation across a scenario grid): `rslt()` (data
simulation), `fit_slt()` (estimation) and `confint()`/`slt_bootstrap_ci()`
(inference) are the complete set of building blocks needed to run one
yourself, and composing them in a plain loop over your own scenario grid
is clearer and more flexible than a maintained built-in engine.
