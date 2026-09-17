# One test per *concern*, each checking all 12 families / 3 cases at once
# (via vapply + a single aggregated expectation) rather than one
# expectation per family/case. This keeps every family and case exercised
# without a combinatorial explosion of near-identical test cases; failures
# still name the offending family/case via `info =`.

test_that("slt_families() lists all supported families", {
  expect_setequal(slt_families(), c("weibull", "powerlindley", "explindley", "expexponential",
                                     "moexponential", "lomax", "expteissier", "loglogistic",
                                     "chen", "lognormal", "gamma", "powermuth"))
})

test_that("get_family() matches case-insensitively and errors on unknown keys", {
  expect_equal(stagelife:::get_family("Weibull")$key, "weibull")
  expect_equal(stagelife:::get_family("WEIBULL")$key, "weibull")
  expect_error(stagelife:::get_family("not_a_family"), "Unknown family")
})

test_that("every baseline family has a valid pdf, an inverting quantile function, and sane boundaries", {
  theta <- c(0.7, 1.4)
  p <- c(0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99)
  check_one <- function(fam) {
    d <- slt_dist(fam)
    # tolerance = 1e-4 deliberately matches stats::integrate()'s own
    # default rel.tol (.Machine$double.eps^0.25 ~= 1.22e-4), not an
    # arbitrary looseness: this is integrate()'s contractual guarantee,
    # so testing tighter would occasionally flake on an integrand
    # integrate() doesn't out-perform its own target on, and testing
    # looser would fail to catch a real regression. (Empirically,
    # achieved error is ~1e-10 to ~1e-14 for all 5 families at this
    # theta -- see audit notes -- so this tolerance has large margin in
    # practice without relying on that margin to pass.)
    integral <- stats::integrate(function(x) d$d(x, lambda = theta[1], alpha = theta[2]),
                                  0, Inf, subdivisions = 500L)$value
    x <- d$q(p, lambda = theta[1], alpha = theta[2])
    p_back <- d$p(x, lambda = theta[1], alpha = theta[2])
    isTRUE(all.equal(integral, 1, tolerance = 1e-4)) &&
      all(diff(x) > 0) &&
      isTRUE(all.equal(p_back, p, tolerance = 1e-6)) &&
      d$p(0, lambda = theta[1], alpha = theta[2]) == 0 &&
      d$q(0, lambda = theta[1], alpha = theta[2]) == 0 &&
      is.infinite(d$q(1, lambda = theta[1], alpha = theta[2]))
  }
  ok <- vapply(slt_families(), check_one, logical(1))
  expect_true(all(ok), info = paste("Failed for:", paste(names(ok)[!ok], collapse = ", ")))
})

test_that("pslt/dslt/qslt behave correctly at the exact p=0/p=1/x=0 boundaries, for every family", {
  theta0 <- c(0.05, 1.3); theta1 <- c(0.2, 1.3); tau1 <- 5
  ok <- vapply(slt_families(), function(fam) {
    qslt(0, fam, theta0, theta1, tau1) == 0 &&
      is.infinite(qslt(1, fam, theta0, theta1, tau1)) &&
      pslt(0, fam, theta0, theta1, tau1) == 0
  }, logical(1))
  expect_true(all(ok), info = paste("Failed for:", paste(names(ok)[!ok], collapse = ", ")))
})

test_that("shape = 1 gives the documented reduced special-case family", {
  lambda <- 0.6; x <- c(0.1, 1, 3, 8)
  wei <- slt_dist("weibull"); ee <- slt_dist("expexponential")
  expect_equal(wei$p(x, lambda = lambda, alpha = 1), 1 - exp(-lambda * x))
  expect_equal(ee$p(x, lambda = lambda, alpha = 1), 1 - exp(-lambda * x))
})

test_that("Exponentiated Teissier handles large x without NaN/Inf (regression test)", {
  fam <- stagelife:::get_family("expteissier")
  x <- c(1, 10, 100, 1000, 1e6)
  pdf_vals <- fam$pdf(x, c(0.7, 1.3)); cdf_vals <- fam$cdf(x, c(0.7, 1.3))
  expect_true(all(is.finite(pdf_vals)) && all(is.finite(cdf_vals)) &&
                isTRUE(all.equal(cdf_vals[length(cdf_vals)], 1, tolerance = 1e-8)))
})

test_that("every family rejects invalid theta as NA rather than a plausible-looking wrong value", {
  # Regression test for an audit finding: un-guarded closed forms could
  # silently return a numerically plausible but meaningless value for
  # invalid theta (e.g. weibull gave a *negative* cdf for lambda < 0;
  # expexponential/lomax gave a cdf > 1 for alpha < 0), which no
  # is.finite()-based downstream safety net would catch.
  bad_thetas <- list(c(-1, 1.4), c(0.7, -1), c(0, 1.4), c(0.7, 0), c(NA, 1.4), c(0.7, NaN))
  ok <- vapply(slt_families(), function(fam_name) {
    fam <- stagelife:::get_family(fam_name)
    all(vapply(bad_thetas, function(theta) {
      is.na(fam$cdf(2.5, theta)) && is.na(fam$pdf(2.5, theta)) && is.na(fam$quantile(0.5, theta))
    }, logical(1)))
  }, logical(1))
  expect_true(all(ok), info = paste("Failed for:", paste(names(ok)[!ok], collapse = ", ")))
})

test_that("theta = Inf lambda remains a valid, non-rejected boundary for Power Lindley", {
  # Must NOT be caught by the invalid-theta guard: lambda -> Inf is a
  # legitimate degenerate (point-mass-at-0) limit that .pl_cdf/.pl_pdf/
  # .qlindley are specifically hardened to handle (see test above this
  # one for the *rejected* cases, and the "Exponentiated Teissier
  # handles large x" test for the analogous case in that family).
  fam <- stagelife:::get_family("powerlindley")
  expect_equal(fam$cdf(5, c(Inf, 1.4)), 1)
  expect_equal(fam$pdf(5, c(Inf, 1.4)), 0)
  expect_equal(fam$quantile(0.5, c(Inf, 1.4)), 0)
})
test_that("Power Lindley handles lambda -> Inf and NaN inputs without corruption (regression test)", {
  fam <- stagelife:::get_family("powerlindley")
  # (1) lambda = Inf used to hit an Inf/Inf indeterminate form in the raw
  #     cdf formula (NaN); it must now saturate to 1, matching the
  #     mathematically correct point-mass-at-0 limit (as Weibull already does).
  expect_equal(fam$cdf(5, c(Inf, 1.4)), 1)
  expect_equal(fam$pdf(5, c(Inf, 1.4)), 0)
  # (2) a NaN probability used to be silently mapped to a quantile of 0
  #     (masking an invalid parameter region as valid); it must propagate NA.
  expect_true(is.na(fam$quantile(NaN, c(0.2, 1.4))))
  # (3) quantile/cdf stay finite and well-ordered across a wide lambda range.
  ok <- vapply(c(1e-6, 1, 1e6, 1e150), function(lambda) {
    x <- fam$quantile(c(0.1, 0.5, 0.9), c(lambda, 1.4))
    all(is.finite(x)) && all(diff(x) > 0)
  }, logical(1))
  expect_true(all(ok))
})

test_that("the two-stage cdf/pdf/quantile round-trip, are continuous at tau1, and integrate to 1 for every family", {
  theta0 <- c(0.05, 1.3); theta1 <- c(0.2, 1.3); tau1 <- 5
  u <- c(0.05, 0.2, 0.4, 0.6, 0.8, 0.95)
  check_one <- function(fam) {
    t <- qslt(u, fam, theta0, theta1, tau1)
    u2 <- pslt(t, fam, theta0, theta1, tau1)
    eps <- 1e-6
    Fm <- pslt(tau1 - eps, fam, theta0, theta1, tau1)
    Fp <- pslt(tau1 + eps, fam, theta0, theta1, tau1)
    dens_integral <- stats::integrate(function(x) dslt(x, fam, theta0, theta1, tau1),
                                       0, Inf, subdivisions = 1000L)$value
    isTRUE(all.equal(u2, u, tolerance = 1e-6)) &&
      isTRUE(all.equal(Fm, Fp, tolerance = 1e-4)) &&
      isTRUE(all.equal(dens_integral, 1, tolerance = 1e-3))
  }
  ok <- vapply(slt_families(), check_one, logical(1))
  expect_true(all(ok), info = paste("Failed for:", paste(names(ok)[!ok], collapse = ", ")))
})

test_that("slt_equiv_age()/slt_tau() return finite, non-negative values for every family", {
  theta0 <- c(0.05, 1.3); theta1 <- c(0.2, 1.3)
  ok <- vapply(slt_families(), function(fam) {
    tau1 <- slt_tau(fam, theta0, p = 0.5)
    nu1 <- slt_equiv_age(fam, tau1, theta0, theta1)
    tau1 > 0 && is.finite(nu1) && nu1 >= 0
  }, logical(1))
  expect_true(all(ok))
})

test_that("slt_cases() exposes exactly the three identifiable cases, and 'full' is never exposed", {
  cs <- slt_cases()
  expect_setequal(cs$case, c("common", "shape1", "shape11"))
  expect_equal(cs$n_par[match(c("common", "shape1", "shape11"), cs$case)], c(3L, 3L, 2L))
  expect_error(stagelife:::slt_expand_params(c(1, 2, 3, 4), case = "full"))
  expect_error(stagelife:::case_info("full"))
})

test_that("slt_expand_params() / slt_reduce_params() are mutual inverses for every case", {
  ok <- vapply(c("common", "shape1", "shape11"), function(case) {
    par <- switch(case, common = c(0.1, 0.2, 1.5), shape1 = c(0.1, 1.5, 0.2), shape11 = c(0.1, 0.2))
    pars <- stagelife:::slt_expand_params(par, case)
    back <- slt_reduce_params(pars$theta0, pars$theta1, case)
    isTRUE(all.equal(unname(back), unname(par)))
  }, logical(1))
  expect_true(all(ok))
})

test_that("slt_reduce_params() warns when 'common' inputs do not actually share a shape", {
  expect_warning(slt_reduce_params(c(0.1, 1.2), c(0.2, 1.5), case = "common"), "shared shape")
})

test_that("Teissier/Power Muth cumulative-hazard helpers are stable across the full precision range (regression test)", {
  # Regression test for an audit finding: the naive exp(z)-1-z formula
  # for the shared Teissier/Power Muth baseline cumulative hazard loses
  # essentially all precision for small z (even wrong-*signed* by
  # z=1e-8), and a second, independent cancellation in 1-exp(-H0)
  # produced spurious Inf once alpha < 1 combined with H0 below machine
  # epsilon. Both were fixed with expm1()-based rewrites, with a direct
  # Taylor-series fallback below z=1e-8 where even expm1(z)-z cannot
  # resolve the z^2/2 correction term in double precision. Checked here
  # via both families that share this code path, across a range
  # spanning the exact z where the naive formula was wrong and the
  # extreme z=1e-16 range where even expm1() alone was insufficient.
  for (fam_name in c("expteissier", "powermuth")) {
    fam <- stagelife:::get_family(fam_name)
    vals <- vapply(c(1e-2, 1e-6, 1e-10, 1e-14, 1e-16), function(x) {
      fam$pdf(x, c(1e-4, 0.1))
    }, numeric(1))
    ok <- all(is.finite(vals)) && all(vals > 0) && all(diff(vals) > 0)
    expect_true(ok, label = fam_name, info = paste(fam_name, ":", paste(vals, collapse = ", ")))
  }
})

test_that("Exponentiated Teissier and Power Muth quantiles use the exact closed Lambert-W form (regression test)", {
  # Regression test for a correctness finding: an earlier version of
  # this package claimed these two families have no closed-form
  # quantile "not even via Lambert W" and inverted their cdfs
  # numerically via uniroot() (correct-valued, but only accurate to
  # uniroot's own tolerance, ~1e-8, and far slower). This was wrong --
  # both have an exact closed form (Sharma, Singh & Shekhawat, 2022,
  # Proposition 4.2, for Exponentiated Teissier; the analogous identity
  # for Power Muth) -- verified directly against the published formula
  # before replacing the numerical version. The round trip below should
  # now be accurate to close to machine precision (~1e-15), not merely
  # to uniroot's ~1e-8 tolerance; a regression back to a numerical
  # quantile would show up here as a much larger error.
  p <- c(0.001, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999)
  for (fam_name in c("expteissier", "powermuth")) {
    d <- slt_dist(fam_name)
    x <- d$q(p, lambda = 0.7, alpha = 1.4)
    p_back <- d$p(x, lambda = 0.7, alpha = 1.4)
    expect_lt(max(abs(p_back - p)), 1e-10)
  }
})

test_that("Power Muth (Power Teissier) at alpha=1 exactly matches Exponentiated Teissier's baseline at alpha=1", {
  # Both families generalise the SAME baseline Teissier(lambda)
  # distribution (power-the-argument vs exponentiate-the-cdf); verified
  # analytically that H(x;lambda,alpha=1) is identical for both.
  pmuth <- stagelife:::get_family("powermuth")
  et <- stagelife:::get_family("expteissier")
  x <- c(0.1, 0.5, 1, 2, 5)
  expect_equal(pmuth$cdf(x, c(0.7, 1)), et$cdf(x, c(0.7, 1)))
  expect_equal(pmuth$pdf(x, c(0.7, 1)), et$pdf(x, c(0.7, 1)))
})

test_that("Gamma and Log-normal match their base-R references directly (sanity anchor)", {
  gamma_fam <- stagelife:::get_family("gamma")
  x <- c(0.5, 1, 2, 5)
  expect_equal(gamma_fam$cdf(x, c(0.7, 1.4)), pgamma(x, shape = 1.4, rate = 0.7))
  expect_equal(gamma_fam$pdf(x, c(0.7, 1.4)), dgamma(x, shape = 1.4, rate = 0.7))
  expect_equal(gamma_fam$quantile(c(0.25, 0.75), c(0.7, 1.4)), qgamma(c(0.25, 0.75), shape = 1.4, rate = 0.7))

  lnorm_fam <- stagelife:::get_family("lognormal")
  expect_equal(lnorm_fam$cdf(x, c(0.7, 1.4)), pnorm(log(0.7 * x) / 1.4))
})

test_that("families reduce to exactly Exponential at alpha=1 exactly when documented (exponential_at_alpha1 flag)", {
  lambda <- 0.7; x <- c(0.5, 1, 2, 5)
  exp_cdf <- 1 - exp(-lambda * x)
  ok <- vapply(slt_families(), function(fam_name) {
    fam <- stagelife:::get_family(fam_name)
    is_exp <- isTRUE(all.equal(fam$cdf(x, c(lambda, 1)), exp_cdf, tolerance = 1e-8))
    is_exp == isTRUE(fam$exponential_at_alpha1)
  }, logical(1))
  expect_true(all(ok), info = paste("Failed for:", paste(names(ok)[!ok], collapse = ", ")))
})

test_that("dslt()/pslt()/qslt() handle NA in the input without crashing, for every family (regression test)", {
  # Regression test for a critical finding from a rigorous audit against
  # the source paper: `ind0 <- x <= tau1; ind1 <- !ind0` is NA wherever
  # x is NA, and `out[NA_logical] <- value` errors ("NAs are not
  # allowed in subscripted assignments") -- confirmed to crash even for
  # weibull, the simplest, most-tested family, so this was not specific
  # to any one distribution; it affected all three of pslt()/dslt()/
  # qslt() for every one of the 12 families. Fixed by computing an
  # explicit ok <- !is.na(x) mask before any branching. A mixed vector
  # (not an all-NA one) is used deliberately, since an all-NA scalar
  # input can misleadingly appear safe due to unrelated lazy-evaluation
  # behaviour elsewhere in the codebase (see the test below).
  for (fam_name in slt_families()) {
    d <- tryCatch(dslt(c(1, NA, 3), fam_name, c(0.5, 1.3), c(0.5, 1.3), tau1 = 2),
                  error = function(e) NULL)
    p <- tryCatch(pslt(c(1, NA, 3), fam_name, c(0.5, 1.3), c(0.5, 1.3), tau1 = 2),
                  error = function(e) NULL)
    q <- tryCatch(qslt(c(0.1, NA, 0.5), fam_name, c(0.5, 1.3), c(0.5, 1.3), tau1 = 2),
                  error = function(e) NULL)
    ok <- !is.null(d) && is.na(d[2]) && all(is.finite(d[-2])) &&
      !is.null(p) && is.na(p[2]) && all(is.finite(p[-2])) &&
      !is.null(q) && is.na(q[2]) && all(is.finite(q[-2]))
    expect_true(ok, label = fam_name)
  }
})

test_that(".teiss_G0()/.teiss_g0() (shared by expteissier and powermuth) handle NA without crashing (regression test)", {
  # Regression test for a second, independent critical finding in the
  # same audit: `if (any(mid))` crashes outright when `mid` contains NA
  # (which it does whenever the input contains NA), the same failure
  # mode as the test above but at a different call site, affecting
  # exactly the two families that call these shared Teissier-baseline
  # helpers directly (expteissier, powermuth). A scalar all-NA test
  # alone is not sufficient to catch this: `ifelse(x<=0, 0, EXPR)`
  # never forces EXPR when x is a length-1 NA (an accident of
  # ifelse()'s internal implementation, not a deliberate safeguard),
  # which let expteissier's cdf/pdf appear NA-safe for that specific
  # input shape while still crashing on any realistic mixed vector.
  for (fam_name in c("expteissier", "powermuth")) {
    fam <- stagelife:::get_family(fam_name)
    scalar_na <- tryCatch(fam$cdf(NA, c(0.5, 1.3)), error = function(e) "CRASH")
    mixed <- tryCatch(fam$cdf(c(1, NA, 3), c(0.5, 1.3)), error = function(e) "CRASH")
    mixed_pdf <- tryCatch(fam$pdf(c(1, NA, 3), c(0.5, 1.3)), error = function(e) "CRASH")
    expect_true(is.na(scalar_na), label = fam_name)
    expect_true(is.numeric(mixed) && is.na(mixed[2]) && all(is.finite(mixed[-2])), label = fam_name)
    expect_true(is.numeric(mixed_pdf) && is.na(mixed_pdf[2]) && all(is.finite(mixed_pdf[-2])), label = fam_name)
  }
})

test_that("weibull pdf(Inf) is exactly 0, not NaN (regression test)", {
  # Regression test for a finding from the same audit: .wei_pdf() was
  # missing the is.finite() guard present in every other family,
  # producing an unguarded Inf*0 = NaN at the literal value x = Inf
  # (every finite, even very large, x already underflowed to exactly 0
  # correctly; only the literal Inf case was wrong). Checked across all
  # 12 families, not just weibull, since this is exactly the kind of
  # gap that could recur in any single family without affecting the
  # others.
  for (fam_name in slt_families()) {
    fam <- stagelife:::get_family(fam_name)
    v <- fam$pdf(Inf, c(0.5, 1.3))
    expect_true(isTRUE(all.equal(v, 0)), label = fam_name)
  }
})
