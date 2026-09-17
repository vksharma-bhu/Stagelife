############################################################################
# Baseline two-parameter lifetime distributions theta = (lambda, alpha):
#   lambda > 0 : rate/scale-type parameter
#   alpha  > 0 : shape/power parameter
#
# Every distribution here has a closed-form pdf/cdf and quantile function;
# the quantile functions of Power Lindley, Exponentiated Lindley, Exponentiated
# Teissier and Power Muth all require the Lambert W function. Only some
# distributions reduce to a *named* special case at alpha = 1 (see each
# distribution's `reduces_to` field and its
# `exponential_at_alpha1` flag in the registry below -- both are
# documentation metadata only, not enforced anywhere in this package).
############################################################################

## ---- Weibull --------------------------------------------------------
.wei_cdf <- function(x, theta) {
  lambda <- theta[1]; alpha <- theta[2]
  ifelse(x <= 0, 0, 1 - exp(-lambda * x^alpha))
}
.wei_pdf <- function(x, theta) {
  lambda <- theta[1]; alpha <- theta[2]
  val <- lambda * alpha * x^(alpha - 1) * exp(-lambda * x^alpha)
  ifelse(x <= 0, 0, ifelse(is.finite(val), val, 0))
}
.wei_qtl <- function(p, theta) {
  lambda <- theta[1]; alpha <- theta[2]
  (-log(1 - p) / lambda)^(1 / alpha)
}

## ---- Power Lindley (Ghitany et al. 2013) -----------------------------
## Quantile of the (ordinary) Lindley distribution via the Lambert W_{-1}
## branch: F(x) = p => (1 + lambda + lambda*x) * exp(-lambda*x) = (1+lambda)(1-p)
##
## Hardened against three numerical failure modes that are specific to
## this family (and that otherwise show up as occasional, hard-to-explain
## crashes deep inside bootstrap/simulation loops -- see
## vignette("stagelife-intro"), section "Power Lindley and numerical
## robustness"):
##   1. lambda -> Inf (via exp(log_par) overflow during an optimiser's
##      wide Nelder-Mead/BFGS search) makes the *raw* cdf formula an
##      Inf/Inf indeterminate form (NaN), whereas the mathematically
##      correct limit is a cdf of 1 (point mass at 0). Unlike Weibull's
##      cdf, which saturates to 1 cleanly in the same limit, the Power
##      Lindley formula's extra (1+lambda+lambda*x)/(1+lambda) term does
##      not simplify away on its own.
##   2. Silently returning 0 for a non-finite (NaN) input probability
##      (rather than propagating NaN) would let an invalid parameter
##      region masquerade as a valid, finite likelihood contribution,
##      instead of being correctly penalised.
##   3. Floating-point noise in an optimiser's finite-difference probes
##      can push the Lambert-W argument z fractionally outside its valid
##      domain [-1/e, 0); lamW::lambertWm1() returns NaN (silently, no
##      error) just outside that domain, which must be prevented rather
##      than allowed to propagate.
.qlindley <- function(p, lambda) {
  res <- rep(NA_real_, length(p))
  finite_p <- is.finite(p)
  res[finite_p & p <= 0] <- 0
  res[finite_p & p >= 1] <- Inf
  idx <- which(finite_p & p > 0 & p < 1)
  if (!length(idx)) return(res)

  if (is.infinite(lambda) && lambda > 0) {
    res[idx] <- 0  # lambda -> Inf collapses the Lindley distribution to a point mass at 0
    return(res)
  }
  if (!is.finite(lambda) || lambda <= 0) return(res)  # invalid region: leave as NA

  a <- 1 + lambda
  if (a > 400) {
    # For lambda beyond about 400, exp(-a) is already many orders of
    # magnitude smaller than double precision can usefully resolve
    # relative to the O(a) terms below, and it underflows to exact 0.0
    # once a exceeds ~745 -- at which point the exact closed form below
    # collapses to NaN/Inf outright (verified: it is already NaN at
    # a = 745, Inf immediately after) rather than merely imprecise. This
    # branch switches, well before that cliff, to the leading-order
    # large-lambda asymptotic of the defining equation
    # (1-p) = [1 + lambda*x/(1+lambda)] * exp(-lambda*x), namely
    # (1-p) ~= exp(-lambda*x), i.e. x ~= -log(1-p)/lambda -- already
    # accurate to a fraction of a percent at a = 400 and improving as
    # lambda grows further. This branch only ever triggers for lambda
    # values far outside any realistic fitted-parameter range (it exists
    # to keep transient optimiser/bootstrap excursions finite and
    # correctly signed rather than to be precise there).
    res[idx] <- -log(1 - p[idx]) / lambda
    return(res)
  }

  z <- -a * (1 - p[idx]) * exp(-a)
  z <- pmin(z, 0)          # guard against fp noise pushing z fractionally positive
  z <- pmax(z, -exp(-1))   # guard against fp noise pushing z fractionally below -1/e
  w <- lamW::lambertWm1(z)
  res[idx] <- -(a + w) / lambda
  res
}
.pl_cdf <- function(x, theta) {
  lambda <- theta[1]; alpha <- theta[2]
  xp <- ifelse(x <= 0, 0, x^alpha)
  val <- 1 - ((1 + lambda + lambda * xp) / (1 + lambda)) * exp(-lambda * xp)
  val <- ifelse(is.finite(val), val, 1)  # non-finite (e.g. lambda -> Inf) saturates to 1
  ifelse(x <= 0, 0, pmin(pmax(val, 0), 1))
}
.pl_pdf <- function(x, theta) {
  lambda <- theta[1]; alpha <- theta[2]
  val <- ifelse(x <= 0, 0,
                (alpha * lambda^2 / (1 + lambda)) * (1 + x^alpha) * x^(alpha - 1) *
                  exp(-lambda * x^alpha))
  ifelse(is.finite(val), val, 0)  # overflow/underflow combinations -> density 0
}
.pl_qtl <- function(p, theta) {
  lambda <- theta[1]; alpha <- theta[2]
  (.qlindley(p, lambda))^(1 / alpha)
}

## ---- Exponentiated Lindley ---------------------------------------------
## F(x;lambda,alpha) = [F0(x;lambda)]^alpha, where F0/f0 are the baseline
## Lindley cdf/pdf (lambda as in Power Lindley above). This is the
## "exponentiate the cdf" sibling of Power Lindley's "power the
## argument" -- exactly as Exponentiated Teissier and Power Muth are the
## two siblings of the Teissier baseline. Reduces to baseline
## Lindley(lambda) at alpha = 1 (not Exponential; see `?slt_dist`).
## Quantile function reuses the same Lambert-W baseline-Lindley inverter
## as Power Lindley (`.qlindley()` above), just inverting a different
## target probability (p^(1/alpha) here, vs p directly for Power
## Lindley's outer power-of-x transform) -- so it inherits that
## function's extreme-lambda hardening automatically.
.el_G0 <- function(x, lambda) {
  ifelse(x <= 0, 0, 1 - ((1 + lambda + lambda * x) / (1 + lambda)) * exp(-lambda * x))
}
.el_g0 <- function(x, lambda) {
  ifelse(x <= 0, 0, (lambda^2 / (1 + lambda)) * (1 + x) * exp(-lambda * x))
}
.el_cdf <- function(x, theta) {
  lambda <- theta[1]; alpha <- theta[2]
  G0 <- .el_G0(x, lambda)
  val <- G0^alpha
  ifelse(x <= 0, 0, ifelse(is.finite(val), val, 1))
}
.el_pdf <- function(x, theta) {
  lambda <- theta[1]; alpha <- theta[2]
  G0 <- .el_G0(x, lambda); g0 <- .el_g0(x, lambda)
  val <- alpha * g0 * G0^(alpha - 1)
  ifelse(x <= 0, 0, ifelse(is.finite(val), val, 0))
}
.el_qtl <- function(p, theta) {
  lambda <- theta[1]; alpha <- theta[2]
  .qlindley(p^(1 / alpha), lambda)
}

## ---- Exponentiated Exponential (Gupta & Kundu 1999) -------------------
.ee_cdf <- function(x, theta) {
  lambda <- theta[1]; alpha <- theta[2]
  ifelse(x <= 0, 0, (1 - exp(-lambda * x))^alpha)
}
.ee_pdf <- function(x, theta) {
  lambda <- theta[1]; alpha <- theta[2]
  ifelse(x <= 0, 0,
         alpha * lambda * exp(-lambda * x) * (1 - exp(-lambda * x))^(alpha - 1))
}
.ee_qtl <- function(p, theta) {
  lambda <- theta[1]; alpha <- theta[2]
  -log(1 - p^(1 / alpha)) / lambda
}

## ---- Marshall-Olkin Exponential (Marshall & Olkin, 1997) ---------------
## F(x;lambda,alpha) = 1 - alpha*e^{-lambda x} / (1 - (1-alpha)*e^{-lambda x}).
## Reduces to Exponential(lambda) exactly at alpha = 1. Constructed by a
## different mechanism than the exponentiated/power-transform families
## above (Marshall & Olkin's "add a tilt parameter" generalisation,
## applicable to any baseline survival function); alpha is the tilt
## parameter. The denominator 1-(1-alpha)*e^{-lambda x} is always
## strictly positive for alpha, lambda, x > 0 (since (1-alpha)*e^{-lambda
## x} < 1 whenever alpha > 0), so unlike several families above this one
## needs no overflow/cancellation guarding beyond the usual x <= 0
## boundary.
.moe_cdf <- function(x, theta) {
  lambda <- theta[1]; alpha <- theta[2]
  e <- exp(-lambda * x)
  val <- 1 - alpha * e / (1 - (1 - alpha) * e)
  ifelse(x <= 0, 0, val)
}
.moe_pdf <- function(x, theta) {
  lambda <- theta[1]; alpha <- theta[2]
  e <- exp(-lambda * x)
  val <- lambda * alpha * e / (1 - (1 - alpha) * e)^2
  ifelse(x <= 0, 0, val)
}
.moe_qtl <- function(p, theta) {
  lambda <- theta[1]; alpha <- theta[2]
  log((1 - p + p * alpha) / (1 - p)) / lambda
}

## ---- Lomax (Pareto Type II) -------------------------------------------
.lx_cdf <- function(x, theta) {
  lambda <- theta[1]; alpha <- theta[2]
  ifelse(x <= 0, 0, 1 - (1 + lambda * x)^(-alpha))
}
.lx_pdf <- function(x, theta) {
  lambda <- theta[1]; alpha <- theta[2]
  ifelse(x <= 0, 0, alpha * lambda * (1 + lambda * x)^(-alpha - 1))
}
.lx_qtl <- function(p, theta) {
  lambda <- theta[1]; alpha <- theta[2]
  ((1 - p)^(-1 / alpha) - 1) / lambda
}

## ---- Exponentiated Teissier --------------------------------------------
## Baseline Teissier (1934) hazard h0(x;lambda) = lambda*(exp(lambda*x) - 1),
## cumulative hazard H0(x) = exp(lambda*x) - 1 - lambda*x, so that the
## baseline survival/CDF is S0(x) = exp(-H0(x)), G0(x) = 1 - S0(x).
## The exponentiated family raises G0 to the power alpha; closed-form pdf,
## cdf AND quantile function (via the Lambert W function) all exist --
## see .et_qtl() below and Sharma, Singh & Shekhawat (2022), Proposition 4.2.
## Guarded against overflow of exp(lambda*x) for large x: once
## lambda*x > 700, exp(lambda*x) overflows double precision, but the
## cumulative hazard H0 is already astronomically large at that point, so
## the true survival/density are (numerically) exactly 0/1 there; we
## short-circuit rather than let Inf*0 propagate as NaN (which would
## abort numerical integration/optimisation). Separately guarded against
## catastrophic cancellation for SMALL lambda*x: the naive
## `exp(z) - 1 - z` loses essentially all significant digits once
## z is small enough (verified: wrong by z=1e-8, even wrong-*signed*),
## since it subtracts two nearly-equal O(1) and O(z) quantities to
## recover an O(z^2) result. Using `expm1(z) - z` (stats::expm1 is the
## numerically stable exp(x)-1) matches the true Taylor series
## (z^2/2 + z^3/6 + ...) down to about z = 1e-12 -- but for z smaller
## than about 1e-8, EVEN `expm1(z) - z` fails: the z^2/2 correction term
## is by then smaller than the ~52-bit relative precision with which z
## itself is representable, so `expm1(z)` cannot resolve it at all and
## the subtraction rounds to exactly 0 (verified directly: z = 1e-16
## gives exactly 0 instead of the true 5e-33). This is a hardware/
## IEEE-754 double-precision floor, not fixable by any expm1()
## implementation; for z below 1e-8 (where the series itself is already
## accurate to > 15 digits, so no accuracy is sacrificed by switching)
## the Taylor series is used directly, which involves no subtraction of
## comparable-magnitude quantities at all. A second, independent
## instance of the same z-cancellation pattern occurs one step later:
## `1 - exp(-H0)` for the resulting tiny H0 rounds to *exactly* 0 once
## H0 is below machine epsilon, which combined with an exponentiated
## shape alpha < 1 (a negative power on G0) produced spurious `Inf`
## from `0^(negative)`; fixed the same way, with `-expm1(-H0)` in place
## of `1 - exp(-H0)` (this one only needs the expm1() swap, not a
## further series fallback, since it is not a *repeated* application of
## the same cancellation).
.h0_stable <- function(z) {
  out <- numeric(length(z))
  out[is.na(z)] <- NA_real_
  ok <- !is.na(z)
  small <- ok & abs(z) < 1e-8
  if (any(small)) {
    zs <- z[small]
    out[small] <- zs^2 / 2 + zs^3 / 6 + zs^4 / 24
  }
  big <- ok & !small
  if (any(big)) out[big] <- expm1(z[big]) - z[big]
  out
}
.teiss_G0 <- function(x, lambda) {
  # NA-safety regression: `pos <- x > 0` is NA wherever x is NA, and the
  # original `if (any(mid))` crashed outright ("missing value where
  # TRUE/FALSE needed") the moment any NA reached it -- discovered via a
  # direct rigorous audit against the paper, not a user report: every
  # other family in this package returns NA cleanly for NA input
  # (matching stats::dnorm()-style conventions), but this shared
  # Teissier-baseline helper (used directly by both expteissier and
  # powermuth) did not, for any input containing NA, not just an
  # all-NA vector. `ok` isolates the non-NA positions before any
  # logical test that feeds a control-flow `if`/`any()`, and NA
  # positions are set explicitly rather than left at numeric()'s
  # default 0.
  out <- numeric(length(x))
  out[is.na(x)] <- NA_real_
  ok <- !is.na(x)
  pos <- ok & x > 0
  big <- pos & (lambda * x > 700)
  out[big] <- 1
  mid <- pos & !big
  if (any(mid)) {
    xm <- x[mid]
    H0 <- .h0_stable(lambda * xm)
    out[mid] <- -expm1(-H0)  # = 1 - exp(-H0), stable for tiny H0 (see comment above)
  }
  out
}
.teiss_g0 <- function(x, lambda) {
  out <- numeric(length(x))
  out[is.na(x)] <- NA_real_
  ok <- !is.na(x)
  pos <- ok & x > 0
  big <- pos & (lambda * x > 700)
  out[big] <- 0
  mid <- pos & !big
  if (any(mid)) {
    xm <- x[mid]
    H0 <- .h0_stable(lambda * xm)
    h0 <- lambda * expm1(lambda * xm)
    out[mid] <- h0 * exp(-H0)
  }
  out
}
.et_cdf <- function(x, theta) {
  lambda <- theta[1]; alpha <- theta[2]
  ifelse(x <= 0, 0, .teiss_G0(x, lambda)^alpha)
}
.et_pdf <- function(x, theta) {
  lambda <- theta[1]; alpha <- theta[2]
  G0 <- .teiss_G0(x, lambda); g0 <- .teiss_g0(x, lambda)
  ifelse(x <= 0, 0, alpha * g0 * G0^(alpha - 1))
}
.et_qtl <- function(p, theta) {
  lambda <- theta[1]; alpha <- theta[2]
  # Closed form (Sharma, Singh & Shekhawat, 2022, Proposition 4.2), NOT
  # a numerical inversion: correcting an earlier version of this package
  # that claimed no closed form exists for this family. Verified
  # directly against the paper's equation (5) and cross-checked against
  # this package's own previous (correct-valued, but needlessly
  # numerical) uniroot-based quantile to agree to ~1e-9 (uniroot's own
  # tolerance) before replacing it.
  z <- (p^(1 / alpha) - 1) / exp(1)
  (1 / lambda) * log(-lamW::lambertWm1(z))
}

## ---- Log-Logistic (Fisk) ------------------------------------------------
## F(x;lambda,alpha) = (lambda x)^alpha / (1 + (lambda x)^alpha). The only
## lifetime distribution here whose hazard is non-monotonic (rises then
## falls for alpha > 1), complementing the monotonic-hazard distributions
## above. Guarded against (lambda x)^alpha overflowing to Inf, which
## would otherwise make the cdf's Inf/Inf and the pdf's Inf/Inf^2 both
## evaluate to NaN instead of the correct limits of 1 and 0.
.llogis_cdf <- function(x, theta) {
  lambda <- theta[1]; alpha <- theta[2]
  v <- ifelse(x <= 0, 0, (lambda * x)^alpha)
  val <- v / (1 + v)
  ifelse(x <= 0, 0, ifelse(is.finite(val), val, 1))
}
.llogis_pdf <- function(x, theta) {
  lambda <- theta[1]; alpha <- theta[2]
  v <- ifelse(x <= 0, 0, (lambda * x)^alpha)
  val <- alpha * lambda * (lambda * x)^(alpha - 1) / (1 + v)^2
  ifelse(x <= 0, 0, ifelse(is.finite(val), val, 0))
}
.llogis_qtl <- function(p, theta) {
  lambda <- theta[1]; alpha <- theta[2]
  (p / (1 - p))^(1 / alpha) / lambda
}

## ---- Chen (2000) ---------------------------------------------------------
## F(x;lambda,alpha) = 1 - exp{lambda(1 - e^{x^alpha})}. alpha appears as
## an exponent on x (the same structural role as a shape parameter in
## Weibull/Lomax/etc.), with alpha = 1 a legitimate special case of the
## family (Chen's own boundary member at shape 1, not a separately named
## distribution -- see `?slt_dist`). Hazard is increasing for alpha >= 1
## and bathtub-shaped for alpha < 1, popular in reliability for exactly
## that flexibility. `1 - e^{x^alpha}` is rewritten as `-expm1(x^alpha)`
## (see the Teissier/Power Muth comment above for why: the naive form
## loses precision for small x^alpha). Also guarded against overflow of
## e^{x^alpha} for large x (cdf needs no guard; pdf does, an Inf*0
## pattern).
.chen_cdf <- function(x, theta) {
  lambda <- theta[1]; alpha <- theta[2]
  H <- lambda * expm1(x^alpha)  # cumulative hazard, stable for small x^alpha
  ifelse(x <= 0, 0, -expm1(-H))  # = 1 - exp(-H), stable for small H too
}
.chen_pdf <- function(x, theta) {
  lambda <- theta[1]; alpha <- theta[2]
  H <- lambda * expm1(x^alpha)
  val <- lambda * alpha * x^(alpha - 1) * exp(x^alpha) * exp(-H)
  ifelse(x <= 0, 0, ifelse(is.finite(val), val, 0))
}
.chen_qtl <- function(p, theta) {
  lambda <- theta[1]; alpha <- theta[2]
  (log(1 - log(1 - p) / lambda))^(1 / alpha)
}

## ---- Log-normal -----------------------------------------------------------
## Parameterised here as F(x;lambda,alpha) = Phi(log(lambda x)/alpha),
## i.e. log(X) ~ N(-log(lambda), alpha^2): 1/lambda is the MEDIAN (not
## the mean) and alpha is the log-scale standard deviation, keeping
## "lambda = a rate/scale-type parameter, alpha = a shape parameter" for
## consistency with every other family here. Built directly on
## stats::pnorm()/dnorm()/qnorm(), which are already robust at extreme
## arguments, so no additional overflow guarding is needed.
.lnorm_cdf <- function(x, theta) {
  lambda <- theta[1]; alpha <- theta[2]
  ifelse(x <= 0, 0, stats::pnorm(log(lambda * x) / alpha))
}
.lnorm_pdf <- function(x, theta) {
  lambda <- theta[1]; alpha <- theta[2]
  ifelse(x <= 0, 0, stats::dnorm(log(lambda * x) / alpha) / (alpha * x))
}
.lnorm_qtl <- function(p, theta) {
  lambda <- theta[1]; alpha <- theta[2]
  exp(alpha * stats::qnorm(p)) / lambda
}

## ---- Gamma -----------------------------------------------------------------
## F(x;lambda,alpha) = pgamma(x, shape = alpha, rate = lambda); reduces
## to Exponential(lambda) exactly at alpha = 1. Built directly on
## stats::pgamma()/dgamma()/qgamma() rather than a hand-derived formula,
## deliberately: those are extensively tested elsewhere and remove an
## entire class of the algebra/overflow bugs found and fixed in this
## package's other, hand-derived families during earlier audits.
.gamma_cdf <- function(x, theta) {
  lambda <- theta[1]; alpha <- theta[2]
  stats::pgamma(x, shape = alpha, rate = lambda)
}
.gamma_pdf <- function(x, theta) {
  lambda <- theta[1]; alpha <- theta[2]
  stats::dgamma(x, shape = alpha, rate = lambda)
}
.gamma_qtl <- function(p, theta) {
  lambda <- theta[1]; alpha <- theta[2]
  stats::qgamma(p, shape = alpha, rate = lambda)
}

## ---- Power Muth (Jodra, Gomez, Jimenez-Gamero & Alba-Fernandez, 2017) /
## Power Teissier ------------------------------------------------------------
## Cumulative hazard H(x;lambda,alpha) = e^{lambda x^alpha} - 1 -
## lambda x^alpha -- structurally IDENTICAL to the Exponentiated
## Teissier baseline above with x replaced by x^alpha (equivalently,
## Jodra et al.'s scale parameterisation with beta = lambda^{-1/alpha}).
## "Power Muth" and "Power Teissier" are the same distribution under
## this correspondence; both names are kept since both are used in the
## literature. This is the "power the argument" sibling of Exponentiated
## Teissier's "exponentiate the cdf" -- exactly as Power Lindley and a
## (hypothetical) Exponentiated Lindley would be two different
## generalisations of the same baseline Lindley(lambda). Verified: at
## alpha = 1 this reduces to exactly the same baseline Teissier(lambda)
## that Exponentiated Teissier also reduces to (checked bit-for-bit).
## Reuses the already-overflow-hardened .teiss_G0()/.teiss_g0() helpers
## above; the only new overflow path is the chain-rule factor
## alpha*x^(alpha-1) in the pdf, guarded the same way as the other new
## families here.
.pmuth_cdf <- function(x, theta) {
  lambda <- theta[1]; alpha <- theta[2]
  xa <- ifelse(x <= 0, 0, x^alpha)
  val <- .teiss_G0(xa, lambda)
  ifelse(x <= 0, 0, val)
}
.pmuth_pdf <- function(x, theta) {
  lambda <- theta[1]; alpha <- theta[2]
  xa <- ifelse(x <= 0, 0, x^alpha)
  g0 <- .teiss_g0(xa, lambda)
  val <- g0 * alpha * x^(alpha - 1)
  ifelse(x <= 0, 0, ifelse(is.finite(val), val, 0))
}
.pmuth_qtl <- function(p, theta) {
  lambda <- theta[1]; alpha <- theta[2]
  # Closed form via the same Lambert W identity as .et_qtl() above
  # (Jodra, Gomez, Jimenez-Gamero & Alba-Fernandez, 2017, use the same
  # approach for the Muth/power Muth quantile): u = lambda*x^alpha
  # solves the same "baseline Teissier" equation as .et_qtl()'s
  # argument does (with its own lambda folded to 1, since it is already
  # absorbed into u), so the same substitution applies with p in place
  # of p^(1/alpha) there. Correcting the same earlier "no closed form"
  # error as .et_qtl(); verified to agree with the previous numerical
  # quantile to ~1e-9 before replacing it.
  u <- log(-lamW::lambertWm1((p - 1) / exp(1)))
  (u / lambda)^(1 / alpha)
}

############################################################################
## Registry
############################################################################

## Guards every family's pdf/cdf/quantile against invalid `theta`
## (wrong length, NA/NaN, non-positive). This closes a validation gap
## found while auditing the package: several of the closed forms above
## have no reason to reject a nonsensical theta on their own (they are
## ordinary arithmetic expressions), so an invalid theta -- e.g. a
## negative lambda, or an NA/NaN reaching here from a corrupted
## optimiser step -- can silently produce a numerically *plausible but
## meaningless* value (a "probability" outside [0, 1], or a specific
## family value inside [0, 1] purely by coincidence) instead of an error
## or NA. That would not be caught by any is.finite()-based safety net
## downstream (all of this package's estimation code treats non-finite
## as invalid, so a finite-but-wrong value is the dangerous case, not
## the safe one). Verified case-by-case before adding this guard:
## `weibull` returns a *negative* cdf for lambda < 0; `expexponential`
## and `lomax` return a cdf > 1 for alpha < 0; `powerlindley` returns
## exactly 1 for any non-finite lambda (an accidental side effect of the
## lambda -> Inf saturation logic below); `expteissier` happens to throw
## an uncontrolled "missing value where TRUE/FALSE needed" error for
## NA/NaN lambda (a control-flow accident, not deliberate validation).
##
## `+Inf` is deliberately NOT rejected here: it is a legitimate
## (degenerate point-mass-at-0) limit that `powerlindley`'s functions in
## particular are specifically hardened to return correct, finite
## results for (see the "extreme lambda" notes above), and rejecting it
## outright would discard that.
.with_theta_guard <- function(fn) {
  force(fn)  # evaluate fn (e.g. fam$cdf) NOW, not lazily on first call -- otherwise,
             # since callers do `fam$cdf <- .with_theta_guard(fam$cdf)`, a lazily
             # forced `fn` would resolve to the *already-reassigned* fam$cdf (the
             # wrapper itself), causing infinite recursion on the first real call.
  function(x, theta) {
    if (length(theta) < 2L || anyNA(theta) || any(theta <= 0)) {
      return(rep(NA_real_, length(x)))
    }
    fn(x, theta)
  }
}

.slt_family_table <- lapply(list(
  weibull = list(
    key = "weibull", name = "Weibull",
    pdf = .wei_pdf, cdf = .wei_cdf, quantile = .wei_qtl,
    reduces_to = "Exponential (rate = lambda) when alpha = 1",
    exponential_at_alpha1 = TRUE,
    param_symbols = c("lambda", "alpha")
  ),
  powerlindley = list(
    key = "powerlindley", name = "Power Lindley",
    pdf = .pl_pdf, cdf = .pl_cdf, quantile = .pl_qtl,
    reduces_to = "Lindley (rate = lambda) when alpha = 1",
    exponential_at_alpha1 = FALSE,
    param_symbols = c("lambda", "alpha")
  ),
  explindley = list(
    key = "explindley", name = "Exponentiated Lindley",
    pdf = .el_pdf, cdf = .el_cdf, quantile = .el_qtl,
    reduces_to = "Lindley (rate = lambda) when alpha = 1",
    exponential_at_alpha1 = FALSE,
    param_symbols = c("lambda", "alpha")
  ),
  expexponential = list(
    key = "expexponential", name = "Exponentiated Exponential",
    pdf = .ee_pdf, cdf = .ee_cdf, quantile = .ee_qtl,
    reduces_to = "Exponential (rate = lambda) when alpha = 1",
    exponential_at_alpha1 = TRUE,
    param_symbols = c("lambda", "alpha")
  ),
  moexponential = list(
    key = "moexponential", name = "Marshall-Olkin Exponential",
    pdf = .moe_pdf, cdf = .moe_cdf, quantile = .moe_qtl,
    reduces_to = "Exponential (rate = lambda) when alpha = 1",
    exponential_at_alpha1 = TRUE,
    param_symbols = c("lambda", "alpha")
  ),
  lomax = list(
    key = "lomax", name = "Lomax",
    pdf = .lx_pdf, cdf = .lx_cdf, quantile = .lx_qtl,
    reduces_to = "Lomax with alpha = 1 (boundary member of the family)",
    exponential_at_alpha1 = FALSE,
    param_symbols = c("lambda", "alpha")
  ),
  expteissier = list(
    key = "expteissier", name = "Exponentiated Teissier",
    pdf = .et_pdf, cdf = .et_cdf, quantile = .et_qtl,
    reduces_to = "(baseline) Teissier (1934) when alpha = 1",
    exponential_at_alpha1 = FALSE,
    param_symbols = c("lambda", "alpha")
  ),
  loglogistic = list(
    key = "loglogistic", name = "Log-Logistic (Fisk)",
    pdf = .llogis_pdf, cdf = .llogis_cdf, quantile = .llogis_qtl,
    reduces_to = "no special reduction at alpha = 1 (still Log-Logistic, shape 1)",
    exponential_at_alpha1 = FALSE,
    param_symbols = c("lambda", "alpha")
  ),
  chen = list(
    key = "chen", name = "Chen (2000)",
    pdf = .chen_pdf, cdf = .chen_cdf, quantile = .chen_qtl,
    reduces_to = "no special reduction at alpha = 1 (still Chen, shape 1)",
    exponential_at_alpha1 = FALSE,
    param_symbols = c("lambda", "alpha")
  ),
  lognormal = list(
    key = "lognormal", name = "Log-normal",
    pdf = .lnorm_pdf, cdf = .lnorm_cdf, quantile = .lnorm_qtl,
    reduces_to = "no special reduction at alpha = 1 (still Log-normal, log-sd 1)",
    exponential_at_alpha1 = FALSE,
    param_symbols = c("lambda", "alpha")
  ),
  gamma = list(
    key = "gamma", name = "Gamma",
    pdf = .gamma_pdf, cdf = .gamma_cdf, quantile = .gamma_qtl,
    reduces_to = "Exponential (rate = lambda) when alpha = 1",
    exponential_at_alpha1 = TRUE,
    param_symbols = c("lambda", "alpha")
  ),
  powermuth = list(
    key = "powermuth", name = "Power Muth (Power Teissier)",
    pdf = .pmuth_pdf, cdf = .pmuth_cdf, quantile = .pmuth_qtl,
    reduces_to = "(baseline) Teissier (1934) when alpha = 1 (same limit as Exponentiated Teissier)",
    exponential_at_alpha1 = FALSE,
    param_symbols = c("lambda", "alpha")
  )
), function(fam) {
  fam$pdf <- .with_theta_guard(fam$pdf)
  fam$cdf <- .with_theta_guard(fam$cdf)
  fam$quantile <- .with_theta_guard(fam$quantile)
  fam
})

#' Available baseline lifetime distributions
#'
#' Lists the two-parameter baseline lifetime distributions, theta = (lambda,
#' alpha), that can be used as the `family` argument throughout the
#' package. See [slt_dist()] for the pdf/cdf and literature source of
#' each.
#'
#' @return A character vector of family keys.
#' @seealso [slt_dist()] for the mathematical definition and reference
#'   of each distribution; [slt_cases()] for the parameter-sharing cases
#'   these distributions are combined under.
#' @export
#' @examples
#' slt_families()
slt_families <- function() {
  names(.slt_family_table)
}

#' Look up a baseline lifetime distribution object
#'
#' @param family One of [slt_families()], matched partially and
#'   case-insensitively (e.g. `"weibull"`, `"Weibull"`, `"pl"` will not
#'   match but `"powerlindley"`/`"PowerLindley"` will).
#' @return A list with elements `key`, `name`, `pdf`, `cdf`, `quantile`,
#'   `reduces_to`, `param_symbols`.
#' @keywords internal
get_family <- function(family) {
  if (is.list(family) && all(c("pdf", "cdf", "quantile") %in% names(family))) {
    return(family)
  }
  key <- tolower(as.character(family))
  key <- gsub("[^a-z0-9]", "", key)
  if (!key %in% names(.slt_family_table)) {
    stop(sprintf(
      "Unknown family '%s'. Available families: %s.",
      family, paste(slt_families(), collapse = ", ")
    ), call. = FALSE)
  }
  .slt_family_table[[key]]
}

#' Baseline density, distribution and quantile functions for a distribution
#'
#' Returns the single-stage `d`/`p`/`q` functions, each of `(x_or_p,
#' lambda, alpha)`, for a baseline lifetime distribution. Useful for exploring a
#' distribution in isolation from the two-stage step-stress structure. In every
#' distribution, `lambda > 0` is a rate/scale-type parameter and `alpha > 0` is
#' a shape/power parameter; see Details for the pdf/cdf of each and its
#' literature source.
#'
#' @details
#' Below, \eqn{F(x;\lambda,\alpha)} and \eqn{f(x;\lambda,\alpha) =
#' F'(x;\lambda,\alpha)} denote the cdf and pdf for \eqn{x > 0} (both are
#' 0 for \eqn{x \le 0}). \eqn{\Phi}/\eqn{\phi} are the standard normal
#' cdf/pdf, and \eqn{\Gamma(\alpha, x) = \int_x^\infty t^{\alpha-1}
#' e^{-t}\,dt} is the upper incomplete gamma function.
#'
#' \strong{Weibull} (Weibull, 1951):
#' \deqn{F(x;\lambda,\alpha) = 1 - \exp(-\lambda x^\alpha), \qquad
#'       f(x;\lambda,\alpha) = \lambda\alpha x^{\alpha-1}\exp(-\lambda x^\alpha).}
#' Reduces to Exponential(rate = \eqn{\lambda}) at \eqn{\alpha = 1}.
#'
#' \strong{Power Lindley} (Ghitany, Al-Mutairi, Balakrishnan & Al-Enezi, 2013):
#' \deqn{F(x;\lambda,\alpha) = 1 - \frac{1+\lambda+\lambda x^\alpha}{1+\lambda}
#'       \exp(-\lambda x^\alpha), \qquad
#'       f(x;\lambda,\alpha) = \frac{\alpha\lambda^2}{1+\lambda}
#'       (1+x^\alpha)\,x^{\alpha-1}\exp(-\lambda x^\alpha).}
#' Reduces to (ordinary) Lindley(\eqn{\lambda}) at \eqn{\alpha = 1} (not
#' Exponential). The quantile function inverts a Lindley-type equation
#' via the \eqn{W_{-1}} branch of the Lambert \eqn{W} function
#' ([lamW::lambertWm1()]): the survival function's numerator-times-
#' exponential (product) structure reduces to the \eqn{w e^w = z} form
#' \eqn{W} solves.
#'
#' \strong{Exponentiated Lindley} (Nadarajah, Bakouch & Tahmasbi, 2011):
#' with baseline Lindley cdf/pdf \eqn{G_0(x;\lambda) = 1 -
#' \frac{1+\lambda+\lambda x}{1+\lambda}\exp(-\lambda x)} and
#' \eqn{g_0(x;\lambda) = \frac{\lambda^2}{1+\lambda}(1+x)\exp(-\lambda x)}
#' (Power Lindley's \eqn{\alpha = 1} case above),
#' \deqn{F(x;\lambda,\alpha) = G_0(x;\lambda)^\alpha, \qquad
#'       f(x;\lambda,\alpha) = \alpha\,g_0(x;\lambda)\,G_0(x;\lambda)^{\alpha-1}.}
#' Reduces to (ordinary) Lindley(\eqn{\lambda}) at \eqn{\alpha = 1} (not
#' Exponential) -- the same limit as \strong{Power Lindley} above, of
#' which this is the "exponentiate the cdf" sibling (exactly as
#' Exponentiated Teissier and Power Muth are the two siblings of the
#' Teissier baseline). Its quantile function reuses the same
#' Lambert-\eqn{W} baseline-Lindley inverter as Power Lindley, just
#' inverting a different target probability.
#'
#' \strong{Exponentiated Exponential} (Gupta & Kundu, 1999):
#' \deqn{F(x;\lambda,\alpha) = \left(1-\exp(-\lambda x)\right)^\alpha, \qquad
#'       f(x;\lambda,\alpha) = \alpha\lambda\exp(-\lambda x)
#'       \left(1-\exp(-\lambda x)\right)^{\alpha-1}.}
#' Reduces to Exponential(\eqn{\lambda}) at \eqn{\alpha = 1}.
#'
#' \strong{Marshall-Olkin Exponential} (Marshall & Olkin, 1997):
#' \deqn{F(x;\lambda,\alpha) = 1 - \frac{\alpha\,e^{-\lambda x}}
#'       {1-(1-\alpha)e^{-\lambda x}}, \qquad
#'       f(x;\lambda,\alpha) = \frac{\lambda\alpha\,e^{-\lambda x}}
#'       {\left(1-(1-\alpha)e^{-\lambda x}\right)^2}.}
#' Reduces to Exponential(\eqn{\lambda}) at \eqn{\alpha = 1}. Constructed
#' by a different mechanism than the exponentiated/power-transform
#' families here -- Marshall & Olkin's "tilt parameter" construction,
#' applicable to any baseline survival function -- with \eqn{\alpha} as
#' the tilt parameter (an odds-ratio-type parameter, interpretable via
#' the family's proportional-odds property).
#'
#' \strong{Lomax} (Lomax, 1954):
#' \deqn{F(x;\lambda,\alpha) = 1 - (1+\lambda x)^{-\alpha}, \qquad
#'       f(x;\lambda,\alpha) = \alpha\lambda(1+\lambda x)^{-\alpha-1}.}
#' At \eqn{\alpha = 1} this is simply the boundary member of the Lomax
#' family itself (not a different named distribution).
#'
#' \strong{Exponentiated Teissier} (Sharma, Singh & Shekhawat, 2022):
#' with baseline Teissier (1934) cdf/pdf
#' \eqn{G_0(x;\lambda) = 1-\exp\{-(e^{\lambda x}-1-\lambda x)\}} and
#' \eqn{g_0(x;\lambda) = \lambda(e^{\lambda x}-1)\exp\{-(e^{\lambda x}-1-\lambda x)\}},
#' \deqn{F(x;\lambda,\alpha) = G_0(x;\lambda)^\alpha, \qquad
#'       f(x;\lambda,\alpha) = \alpha\,g_0(x;\lambda)\,G_0(x;\lambda)^{\alpha-1},}
#' with quantile function (Sharma, Singh & Shekhawat, 2022, Proposition 4.2)
#' \deqn{Q(p;\lambda,\alpha) = \frac{1}{\lambda}
#'       \log\!\left(-W_{-1}\!\left(\frac{p^{1/\alpha}-1}{e}\right)\right),}
#' where \eqn{W_{-1}} is the lower (\eqn{W_{-1}}) branch of the Lambert
#' \eqn{W} function ([lamW::lambertWm1()]). Reduces to baseline
#' Teissier(\eqn{\lambda}) at \eqn{\alpha = 1} (the same limit as
#' \strong{Power Muth} below, of which this is the "exponentiate the
#' cdf" sibling). The Teissier distribution is also known as the Muth
#' distribution; see Power Muth below, both verified against the cited
#' paper's own closed-form quantile.
#'
#' \strong{Log-Logistic / Fisk} (Fisk, 1961):
#' \deqn{F(x;\lambda,\alpha) = \frac{(\lambda x)^\alpha}{1+(\lambda x)^\alpha}, \qquad
#'       f(x;\lambda,\alpha) = \frac{\alpha\lambda(\lambda x)^{\alpha-1}}
#'       {\left(1+(\lambda x)^\alpha\right)^2}.}
#' The only family here with a non-monotonic (rise-then-fall) hazard
#' rate for \eqn{\alpha > 1}; at \eqn{\alpha = 1} it is simply
#' Log-Logistic with shape 1 (not a different named distribution).
#'
#' \strong{Chen} (Chen, 2000):
#' \deqn{F(x;\lambda,\alpha) = 1 - \exp\left(\lambda\left(1-e^{x^\alpha}\right)\right), \qquad
#'       f(x;\lambda,\alpha) = \lambda\alpha x^{\alpha-1} e^{x^\alpha}
#'       \exp\left(\lambda\left(1-e^{x^\alpha}\right)\right).}
#' Hazard rate is increasing for \eqn{\alpha \ge 1} and bathtub-shaped
#' for \eqn{\alpha < 1}.
#'
#' \strong{Log-normal}: parameterised here as
#' \deqn{F(x;\lambda,\alpha) = \Phi\!\left(\frac{\log(\lambda x)}{\alpha}\right), \qquad
#'       f(x;\lambda,\alpha) = \frac{1}{\alpha x}\phi\!\left(\frac{\log(\lambda x)}{\alpha}\right),}
#' i.e. \eqn{\log X \sim N(-\log\lambda,\ \alpha^2)}: \eqn{1/\lambda} is
#' the \strong{median} (not the mean) and \eqn{\alpha} is the log-scale
#' standard deviation, keeping "\eqn{\lambda} a rate/scale parameter,
#' \eqn{\alpha} a shape parameter" consistent with every other family
#' here. Built directly on [stats::plnorm()]-equivalent expressions
#' ([stats::pnorm()]/[stats::dnorm()]/[stats::qnorm()]) rather than a
#' hand-derived formula.
#'
#' \strong{Gamma}: standard shape-rate parameterisation,
#' \deqn{F(x;\lambda,\alpha) = 1 - \frac{\Gamma(\alpha,\lambda x)}{\Gamma(\alpha)}, \qquad
#'       f(x;\lambda,\alpha) = \frac{\lambda^\alpha}{\Gamma(\alpha)}
#'       x^{\alpha-1} e^{-\lambda x}.}
#' Reduces to Exponential(\eqn{\lambda}) at \eqn{\alpha = 1}. Built
#' directly on [stats::pgamma()]/[stats::dgamma()]/[stats::qgamma()]
#' rather than a hand-derived formula.
#'
#' \strong{Power Muth} (also known as \strong{Power Teissier}; Jodra,
#' Gomez, Jimenez-Gamero & Alba-Fernandez, 2017, generalising the
#' baseline Muth/Teissier distribution of Muth, 1977, and Teissier, 1934,
#' and itself studied via the Lambert \eqn{W} function by Jodra,
#' Jimenez-Gamero & Alba-Fernandez, 2015):
#' \deqn{F(x;\lambda,\alpha) = 1-\exp\left(-\left(e^{\lambda x^\alpha}-1-\lambda x^\alpha\right)\right), \qquad
#'       f(x;\lambda,\alpha) = \lambda\alpha x^{\alpha-1}\left(e^{\lambda x^\alpha}-1\right)
#'       \exp\left(-\left(e^{\lambda x^\alpha}-1-\lambda x^\alpha\right)\right),}
#' with quantile function
#' \deqn{Q(p;\lambda,\alpha) = \left(\frac{1}{\lambda}
#'       \log\!\left(-W_{-1}\!\left(\frac{p-1}{e}\right)\right)\right)^{1/\alpha},}
#' using the same Lambert-\eqn{W} identity as Exponentiated Teissier
#' above (\eqn{p} in place of \eqn{p^{1/\alpha}} there, since the power
#' transform is applied to \eqn{x} here rather than to the cdf). Reduces
#' to baseline Teissier/Muth(\eqn{\lambda}) at \eqn{\alpha = 1} -- the
#' same limit as \strong{Exponentiated Teissier} above, of which this is
#' the "power the argument" sibling: both generalise the same baseline
#' distribution, via \eqn{x \mapsto x^\alpha} here versus \eqn{F_0
#' \mapsto F_0^\alpha} there (exactly as Power Lindley and Exponentiated
#' Exponential are, respectively, the "power the argument" and
#' "exponentiate the cdf" generalisations of their own baselines).
#'
#' @param family A family key from [slt_families()].
#' @return A list with functions `d`, `p`, `q`.
#' @references
#' Chen, Z. (2000). A new two-parameter lifetime distribution with
#' bathtub shape or increasing failure rate function. *Statistics &
#' Probability Letters*, 49(2), 155-161.
#'
#' Fisk, P. R. (1961). The graduation of income distributions.
#' *Econometrica*, 29(2), 171-185.
#'
#' Ghitany, M. E., Al-Mutairi, D. K., Balakrishnan, N., and Al-Enezi,
#' L. J. (2013). Power Lindley distribution and associated inference.
#' *Computational Statistics & Data Analysis*, 64, 20-33.
#'
#' Gupta, R. D. and Kundu, D. (1999). Generalized exponential
#' distributions. *Australian & New Zealand Journal of Statistics*,
#' 41(2), 173-188.
#'
#' Jodra, P., Gomez, H. W., Jimenez-Gamero, M. D., and Alba-Fernandez,
#' M. V. (2017). The power Muth distribution. *Mathematical Modelling
#' and Analysis*, 22(2), 186-201.
#'
#' Jodra, P., Jimenez-Gamero, M. D., and Alba-Fernandez, M. V. (2015).
#' On the Muth distribution. *Mathematical Modelling and Analysis*,
#' 20(3), 291-310.
#'
#' Lomax, K. S. (1954). Business failures: Another example of the
#' analysis of failure data. *Journal of the American Statistical
#' Association*, 49(268), 847-852.
#'
#' Marshall, A. W. and Olkin, I. (1997). A new method for adding a
#' parameter to a family of distributions with application to the
#' exponential and Weibull families. *Biometrika*, 84(3), 641-652.
#'
#' Nadarajah, S., Bakouch, H. S., and Tahmasbi, R. (2011). A generalized
#' Lindley distribution. *Sankhya B*, 73(2), 331-359.
#'
#' Muth, E. J. (1977). Reliability models with positive memory derived
#' from the mean residual life function. In *The Theory and
#' Applications of Reliability*, Vol. 2 (C. P. Tsokos and I. N. Shimi,
#' eds.), 401-435. Academic Press, New York.
#'
#' Sharma, V. K., Singh, S. V., and Shekhawat, K. (2022). Exponentiated
#' Teissier distribution with increasing, decreasing and bathtub hazard
#' functions. *Journal of Applied Statistics*, 49(2), 371-393.
#' \doi{10.1080/02664763.2020.1813694}
#'
#' Teissier, G. (1934). Recherches sur le vieillissement et sur les lois
#' de la mortalite. *Annales de Physiologie et de Physicochimie
#' Biologique*, 10, 237-284.
#'
#' Weibull, W. (1951). A statistical distribution function of wide
#' applicability. *Journal of Applied Mechanics*, 18(3), 293-297.
#' @seealso [slt_families()] for the list of family keys; [dslt()],
#'   [pslt()], [qslt()], [rslt()] for the two-stage (not baseline)
#'   population functions built from these; [slt_equiv_age()] for how a
#'   baseline distribution's cdf/quantile combine across the stage change.
#' @export
#' @examples
#' fam <- slt_dist("weibull")
#' fam$p(1, lambda = 1, alpha = 2)
#' fam$q(0.5, lambda = 1, alpha = 2)
slt_dist <- function(family) {
  fam <- get_family(family)
  list(
    d = function(x, lambda, alpha) fam$pdf(x, c(lambda, alpha)),
    p = function(x, lambda, alpha) fam$cdf(x, c(lambda, alpha)),
    q = function(p, lambda, alpha) fam$quantile(p, c(lambda, alpha))
  )
}
