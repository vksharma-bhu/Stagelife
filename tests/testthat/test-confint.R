fit_example_ci <- function(seed = 21, n = 300) {
  set.seed(seed)
  theta0 <- c(0.08, 1.4); theta1 <- c(0.25, 1)
  tau1 <- slt_tau("weibull", theta0, 0.5)
  d <- rslt(n, "weibull", theta0, theta1, tau1, prop_transfer = 1)
  fit_slt(d, "weibull", "shape1", tau1 = tau1)
}

test_that("confint.slt_fit gives valid, comparable Wald intervals and supports subsetting", {
  fit <- fit_example_ci()
  ci_ln <- confint(fit, type = "lognormal")
  ci_n <- confint(fit, type = "normal")
  expect_true(all(ci_ln[, 1] > 0) && all(ci_ln[, 1] <= fit$par & fit$par <= ci_ln[, 2]))
  expect_equal(ci_ln, ci_n, tolerance = 0.05)
  expect_equal(rownames(confint(fit, parm = c("lambda0", "lambda1"))), c("lambda0", "lambda1"))
})

test_that("confint.slt_fit errors informatively on a non-converged fit", {
  fake <- structure(list(converged = FALSE), class = "slt_fit")
  expect_error(confint(fake), "did not converge")
})

test_that("confint.slt_fit rejects an invalid level rather than silently returning a backwards CI", {
  # Regression test for an audit finding: level outside (0,1) previously
  # produced either a cryptic "NaNs produced" warning or, worse, a
  # *silently* backwards interval (lower bound > upper bound) with no
  # warning at all, for e.g. a fat-fingered level = 95 instead of 0.95.
  fit <- fit_example_ci()
  for (bad_level in list(-0.2, 0, 1, 1.5, 95, NA, "0.95", c(0.9, 0.95))) {
    expect_error(confint(fit, level = bad_level), "level")
  }
})

test_that("slt_bootstrap_ci rejects an invalid level", {
  fit <- fit_example_ci()
  expect_error(slt_bootstrap_ci(fit, B = 10, level = -0.2), "level")
  expect_error(slt_bootstrap_ci(fit, B = 10, level = 95), "level")
})

test_that("slt_bootstrap_ci does not spuriously require more successes than B replicates", {
  skip_on_cran()
  fit <- fit_example_ci()
  set.seed(101)
  b <- slt_bootstrap_ci(fit, B = 20, method = "nonparametric")  # B < the old hardcoded floor of 30
  expect_s3_class(b, "slt_boot_ci")
  expect_true(b$B_converged <= 20 && b$B_converged >= 1)
})

test_that("slt_bootstrap_ci returns valid, ordered intervals under both methods and both rslt() designs", {
  skip_on_cran()
  fit <- fit_example_ci()
  configs <- list(
    list(method = "nonparametric", prop_transfer = 1),
    list(method = "parametric",    prop_transfer = 1),
    list(method = "parametric",    prop_transfer = 0.5)
  )
  set.seed(100)
  ok <- vapply(configs, function(cfg) {
    b <- slt_bootstrap_ci(fit, B = 60, method = cfg$method, prop_transfer = cfg$prop_transfer)
    inherits(b, "slt_boot_ci") && all(b$percentile$lower <= b$percentile$upper) &&
      all(b$boot_t$lower <= b$boot_t$upper) &&
      all(b$log_boot_t$lower <= b$log_boot_t$upper) &&
      all(b$log_boot_t$lower > 0) && b$B_converged > 0
  }, logical(1))
  expect_true(all(ok))

  # design = "typeM" path too
  set.seed(102)
  b_m <- slt_bootstrap_ci(fit, B = 60, method = "parametric", design = "typeM", R0_1 = 50)
  expect_s3_class(b_m, "slt_boot_ci")
  expect_true(all(b_m$percentile$lower <= b_m$percentile$upper))
})

test_that("a rare per-replicate error (e.g. a Power Lindley numerical edge case) is skipped, not fatal, on both the single- and multi-core code paths", {
  skip_on_cran()
  fit <- fit_example_ci()
  orig_rslt <- stagelife::rslt
  count <- 0
  # NOTE: patched() intentionally keeps its *original* closure (the local
  # `count`/`orig_rslt` above) -- assignInNamespace() only needs to change
  # which function object the *name* `rslt` is bound to inside the package
  # namespace; it must not touch patched()'s own environment, or its free
  # variables would no longer resolve to the `count`/`orig_rslt` defined
  # here.
  patched <- function(...) {
    count <<- count + 1
    if (count %% 5 == 0) stop("simulated rare numerical edge-case failure")
    orig_rslt(...)
  }
  assignInNamespace("rslt", patched, ns = "stagelife")
  on.exit(assignInNamespace("rslt", orig_rslt, ns = "stagelife"), add = TRUE)

  for (nc in c(1L, 2L)) {
    count <- 0
    b <- tryCatch(
      slt_bootstrap_ci(fit, B = 40, method = "parametric", n_cores = nc),
      error = function(e) e
    )
    expect_false(inherits(b, "error"), info = paste("n_cores =", nc))
    expect_true(inherits(b, "slt_boot_ci"), info = paste("n_cores =", nc))
    expect_true(b$B_converged < 40, info = paste("n_cores =", nc))  # some replicates were skipped
  }
})
