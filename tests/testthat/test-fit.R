case_npar <- function(case) switch(case, common = 3L, shape1 = 3L, shape11 = 2L)

test_that("fit_slt recovers known parameters at moderate n, with sane SEs/logLik/AIC", {
  set.seed(10)
  theta0 <- c(0.08, 1.4); theta1 <- c(0.25, 1)
  tau1 <- slt_tau("weibull", theta0, 0.5)
  d <- rslt(500, "weibull", theta0, theta1, tau1, prop_transfer = 1)
  fit <- fit_slt(d, "weibull", "shape1", tau1 = tau1)
  true_par <- slt_reduce_params(theta0, theta1, "shape1")
  expect_true(fit$converged)
  expect_equal(unname(fit$par), unname(true_par), tolerance = 0.15)
  expect_true(all(fit$se > 0) && is.finite(fit$loglik) &&
                isTRUE(all.equal(fit$aic, -2 * fit$loglik + 2 * fit$npar)))
})

test_that("fit_slt never reports converged=TRUE while stuck at the internal penalty sentinel", {
  # Regression test for an audit finding: optim() can report
  # convergence=0 ("successful") while every point it visited returned
  # neg_loglik_slt()'s penalty value (.Machine$double.xmax/2). A "fit"
  # whose reported log-likelihood is anywhere near -1e307 is never real
  # for any actual dataset, so it must be reported as non-convergence,
  # not success. chen with case="shape1" on this real, large-scale
  # (tau1 ~ 500) data set is a genuine, structural instance of this,
  # not merely a bad starting value: case="shape1" fixes the
  # elevated-stage shape at exactly alpha1=1, and Chen's shape is an
  # exponent on x (x^alpha) feeding directly into an exponential, so
  # alpha1=1 makes the elevated-stage density numerically zero at the
  # observed failure times for *any* lambda1 -- no initial value can
  # fix this, since it is the model itself that is degenerate here, not
  # the optimiser's starting point. See the test below for the fix that
  # *is* possible: case="common", where the shape is shared and freely
  # estimated rather than fixed at 1.
  d_real <- list(Y = slt_han_kundu_2014$y0, Z = slt_han_kundu_2014$y1)
  fit_c <- fit_slt(d_real, "chen", "shape1", tau1 = slt_han_kundu_2014$tau1, n_starts = 3)
  expect_false(fit_c$converged)
  # and, for every family/case where a fit IS reported, its logLik must
  # be a real, bounded value, never anywhere near the penalty sentinel
  combos <- expand.grid(family = slt_families(), case = c("common", "shape1", "shape11"),
                         stringsAsFactors = FALSE)
  sane <- mapply(function(fam, case) {
    fit <- tryCatch(fit_slt(d_real, fam, case, tau1 = slt_han_kundu_2014$tau1, n_starts = 2),
                     error = function(e) NULL)
    is.null(fit) || !isTRUE(fit$converged) || abs(fit$loglik) < 1e6
  }, combos$family, combos$case)
  expect_true(all(sane))
})

test_that(".default_init() calibrates the initial shape against the data, avoiding cdf saturation at large scales (regression test)", {
  # Regression test for a fix: a blind alpha = 1 starting guess left
  # the cdf already saturated (F0(tau1) == 1 exactly) for shape
  # parameters that appear as an exponent feeding an exponential
  # (Chen's x^alpha) on data at this scale, trapping every multi-start
  # attempt at the penalty sentinel before optimisation could even
  # begin (the log-scale perturbation used for extra starting values is
  # far too narrow to escape this). .default_init() now calibrates the
  # initial shape via a data-driven root-find instead. This is only
  # possible where the shape is free (case %in% c("common", "shape1"));
  # it cannot rescue case="shape1" itself here, since alpha1 = 1 is
  # degenerate at this scale regardless of lambda1 (see the test above)
  # -- but case="common" (shared, freely-estimated shape) converges
  # cleanly.
  d_real <- list(Y = slt_han_kundu_2014$y0, Z = slt_han_kundu_2014$y1)
  tau1 <- slt_han_kundu_2014$tau1
  init <- stagelife:::.default_init(d_real, tau1, "shape1", "chen")
  fam <- stagelife:::get_family("chen")
  F0_tau1 <- fam$cdf(tau1, exp(init[1:2]))
  expect_true(F0_tau1 > 0.01 && F0_tau1 < 0.99)

  fit <- fit_slt(d_real, "chen", "common", tau1 = tau1, n_starts = 3)
  expect_true(fit$converged)
})

test_that("fit_slt converges with the right number of parameters for every family x case combination", {
  set.seed(11)
  combos <- expand.grid(family = slt_families(), case = c("common", "shape1", "shape11"),
                         stringsAsFactors = FALSE)
  ok <- mapply(function(fam, case) {
    theta0 <- switch(case, common = c(0.08, 1.4), shape1 = c(0.08, 1.4), shape11 = c(0.08, 1))
    theta1 <- switch(case, common = c(0.25, 1.4), shape1 = c(0.25, 1), shape11 = c(0.25, 1))
    tau1 <- slt_tau(fam, theta0, 0.5)
    d <- rslt(300, fam, theta0, theta1, tau1, prop_transfer = 1)
    fit <- fit_slt(d, fam, case, tau1 = tau1, n_starts = 3)
    fit$converged && length(fit$par) == case_npar(case)
  }, combos$family, combos$case)
  names(ok) <- paste(combos$family, combos$case)
  expect_true(all(ok), info = paste("Failed for:", paste(names(ok)[!ok], collapse = ", ")))
})

test_that("fit_slt validates its inputs", {
  d <- list(Y = c(1, 2, 3), Z = c(5, 6))
  expect_error(fit_slt(d, "weibull", "shape1", tau1 = -1), "tau1")
  expect_error(fit_slt(list(Y = c(1, 2)), "weibull", "shape1", tau1 = 3), "Y and Z")
  expect_error(fit_slt(d, "weibull", "shape1", tau1 = 4, init = c(1, 2)), "length")
  expect_error(fit_slt(d, "weibull", "not_a_case", tau1 = 4))
  expect_error(fit_slt(d, "not_a_family", "shape1", tau1 = 4), "Unknown family")
  expect_error(fit_slt(list(Y = numeric(0), Z = numeric(0)), "weibull", "shape1", tau1 = 5), "observations")
})

test_that("fit_slt rejects negative/zero/non-finite observed failure times", {
  # Regression test for an audit finding: a data-entry error (a negative
  # or zero "failure time", impossible by construction since both stages
  # have positive support) previously converged silently with no error.
  expect_error(fit_slt(list(Y = c(1, 2, -3), Z = c(10, 11)), "weibull", "shape1", tau1 = 5),
               "data\\$Y")
  expect_error(fit_slt(list(Y = c(1, 2, 0), Z = c(10, 11)), "weibull", "shape1", tau1 = 5),
               "data\\$Y")
  expect_error(fit_slt(list(Y = c(1, 2, 3), Z = c(10, -11)), "weibull", "shape1", tau1 = 5),
               "data\\$Z")
  expect_error(fit_slt(list(Y = c(1, 2, NA), Z = c(10, 11)), "weibull", "shape1", tau1 = 5),
               "data\\$Y")
})

test_that("fit_slt rejects an empty Y or Z rather than silently 'converging' to an arbitrary, unidentified estimate", {
  # Regression test for a HIGH-SEVERITY audit finding: with Z empty, the
  # stage-1 log-likelihood contribution is identically 0 for every value
  # of lambda1 (a perfectly flat objective), so optim() never moves and
  # simply reports back whatever the (arbitrary) starting value was as
  # if it were a real MLE, with converged = TRUE. This previously passed
  # the total-count check (n_obs >= cinfo$npar) because that check does
  # not verify each *stage* individually has data. Verified directly:
  # neg_loglik_slt() was bit-for-bit identical across lambda1 spanning
  # 4 orders of magnitude when Z was empty.
  d_empty_z <- list(Y = c(1, 2, 3, 4, 5), Z = numeric(0))
  expect_error(fit_slt(d_empty_z, "weibull", "shape11", tau1 = 5, n_starts = 1),
               "at least one observed")

  d_empty_y <- list(Y = numeric(0), Z = c(6, 7, 8, 9, 10))
  expect_error(fit_slt(d_empty_y, "weibull", "shape11", tau1 = 5, n_starts = 1),
               "at least one observed")
})

test_that("print.slt_fit reports non-convergence clearly, and multi-start never hurts", {
  fake <- structure(list(converged = FALSE, family = "weibull", family_name = "Weibull",
                          case = "shape1", tau1 = 5, data = list(Y = numeric(0), Z = numeric(0))),
                     class = "slt_fit")
  expect_output(print(fake), "DID NOT CONVERGE")

  set.seed(12)
  theta0 <- c(0.08, 1.4); theta1 <- c(0.25, 1)
  tau1 <- slt_tau("weibull", theta0, 0.5)
  d <- rslt(200, "weibull", theta0, theta1, tau1, prop_transfer = 1)
  fit1 <- fit_slt(d, "weibull", "shape1", tau1 = tau1, n_starts = 1)
  fit5 <- fit_slt(d, "weibull", "shape1", tau1 = tau1, n_starts = 5)
  expect_true(fit5$loglik >= fit1$loglik - 1e-6)
})

test_that("S3 methods (print, summary, coef, vcov, logLik, AIC, BIC) all dispatch correctly", {
  set.seed(1); theta0 <- c(0.08, 1.4); theta1 <- c(0.25, 1)
  tau1 <- slt_tau("weibull", theta0, 0.5)
  d <- rslt(300, "weibull", theta0, theta1, tau1, prop_transfer = 1)
  fit <- fit_slt(d, "weibull", "shape1", tau1 = tau1)

  expect_output(print(fit), "Weibull")
  s <- summary(fit)
  expect_s3_class(s, "summary.slt_fit")
  expect_output(print(s), "Coefficients")

  cf <- coef(fit)
  expect_equal(names(cf), c("lambda0", "alpha0", "lambda1"))

  vc <- vcov(fit)
  expect_true(all(dim(vc) == c(3, 3)) && all(diag(vc) > 0))
  expect_equal(diag(vc), (fit$par * sqrt(diag(fit$vcov_log)))^2, tolerance = 1e-8)

  ll <- logLik(fit)
  expect_true(isTRUE(all.equal(as.numeric(ll), fit$loglik)) && attr(ll, "df") == 3L)
  # There is no nobs.slt_fit() method (removed: not needed by anything --
  # AIC()/BIC() read sample size from logLik()'s own "nobs" attribute,
  # checked directly here; a bare nobs(fit) call correctly falls through
  # to stats::nobs.default() and errors, since fit$n (not fit$nobs) is
  # this object's field name).
  expect_equal(attr(ll, "nobs"), fit$n)
  expect_error(nobs(fit), "no 'nobs' method")
  expect_equal(AIC(fit), fit$aic, tolerance = 1e-8)
  expect_equal(BIC(fit), fit$bic, tolerance = 1e-8)
})
