#' stagelife: Random Number Generation and Inferences for Two-Stage Life-Testing Data
#'
#' Tools for simulating, fitting, and validating the stage life testing
#' (SLT) model of Laumen & Cramer (2019).
#'
#' @details
#' Implements a two-stage cumulative-exposure model for stage life
#' testing data, supplied here across twelve baseline lifetime
#' distributions and three parameter-sharing cases (see
#' [slt_families()] and [slt_cases()]). See the package vignette
#' (`vignette("stagelife-intro", package = "stagelife")`) for a full
#' worked workflow.
#'
#' [rslt()] (data generation), [fit_slt()] (estimation), and
#' [confint()] / [slt_bootstrap_ci()] (inference) are the building
#' blocks needed to compose a simulation study yourself, in a plain
#' loop over your own scenario grid; this package does not include a
#' built-in simulation-study engine.
#'
#' @section No right-censoring:
#' The underlying model (Laumen & Cramer's Procedure 2.3, and every
#' function in this package built on it) assumes every one of the `n`
#' units is eventually observed to fail -- either on stage 0 (an early
#' failure, or a survivor never selected for transfer) or on stage 1 (a
#' transferred survivor). There is no mechanism here for a unit that is
#' still functioning when a finite-duration experiment ends
#' (right-censoring); every `Y`/`Z` value supplied to [fit_slt()] or
#' [rslt()]'s output is treated as an exact, observed failure time.
#' Real step-stress experiments often do have such unfailed units at
#' the end of testing, and some of the bundled real data sets
#' illustrate this directly: the original Han & Kundu (2014) experiment
#' placed `n = 35` units on test, but only 31 failure times are
#' reported in the source paper (`?slt_han_kundu_2014`); the original
#' Zhu (2010) experiment used two sets of 32 bulbs (`n = 64`), but only
#' 53 failure times are bundled here (`?slt_zhu_2010`) -- consistent
#' with (though not confirmed by this package as definitely caused by)
#' units that survived to the end of testing being excluded, rather
#' than contributing a censored-likelihood term, when these data sets
#' were adapted to a framework with no censoring mechanism.
#'
#' @section Main functions:
#' \describe{
#'   \item{Population distribution}{[dslt()], [pslt()], [qslt()]}
#'   \item{Simulation}{[rslt()]}
#'   \item{Estimation}{[fit_slt()], plus S3 methods
#'     [print.slt_fit()], [summary.slt_fit()], [coef.slt_fit()],
#'     [vcov.slt_fit()], [logLik.slt_fit()]}
#'   \item{Inference}{[confint.slt_fit()], [slt_bootstrap_ci()]}
#'   \item{Experimental design}{[slt_prob_stage1()]}
#'   \item{Goodness-of-fit / model comparison}{[slt_gof_ks()], [slt_lrt()],
#'     [slt_model_comparison()]}
#'   \item{Diagnostic plots}{[plot_cdf()], [plot_survival()],
#'     [plot_pp()], [plot_qq()], [plot.slt_fit()]}
#' }
#'
#' @references
#' Laumen, B. and Cramer, E. (2019). Stage life testing. *Naval Research
#' Logistics*, 66(8), 632-647. \doi{10.1002/nav.21874}
#'
#' See [slt_dist()] for the literature source of each of the twelve
#' baseline lifetime distributions.
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom stats optim qnorm quantile runif sd rnorm pchisq ks.test ppoints
#' @importFrom stats confint coef vcov logLik AIC BIC
#' @importFrom graphics plot
#' @importFrom parallel mclapply
#' @importFrom compiler cmpfun
## usethis namespace: end
NULL

utils::globalVariables(c(".data"))
