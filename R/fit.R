############################################################################
# Maximum likelihood estimation: Nelder-Mead warm start (derivative-free,
# so it never triggers the finite-difference edge cases that can occur
# when a boundary of a family's support is probed) followed by BFGS for a
# fast, accurate final optimum and numerical Hessian. Multiple starting
# values are tried and the converged fit with the highest log-likelihood
# is kept, guarding against local optima.
############################################################################

.default_init <- function(data, tau1, case, family) {
  Y <- data$Y; Z <- data$Z
  mY <- mean(Y)
  if (!is.finite(mY) || mY <= 0) mY <- 1
  mZ <- mean(pmax(Z - tau1, 1e-3))
  if (!is.finite(mZ) || mZ <= 0) mZ <- mY
  lambda0 <- 1 / mY
  lambda1 <- 1 / (mZ + 1 / lambda0)

  # Calibrate the initial shape guess against the data via a bounded
  # 1-D root-find, rather than a blind alpha = 1 guess: for families
  # where alpha appears as an exponent on x feeding directly into an
  # exponential (e.g. Chen's x^alpha), alpha = 1 can leave the cdf
  # already saturated at tau1 for large-scale data (verified: Chen's
  # F0(tau1) ~ 1 at alpha0 = 1 for data with tau1 ~ 500), silently
  # trapping every multi-start attempt before optimisation even
  # begins -- the log-scale perturbation used for extra starting
  # values is far too narrow to escape this. Calibrating alpha so that
  # cdf(tau1; lambda0, alpha) matches the empirical fraction of
  # Y <= tau1 is family-agnostic and fixes this without special-casing
  # any one family; falls back to 1 if the root-find fails (e.g.
  # degenerate data) or the family has no free shape at all (shape11).
  alpha0 <- 1
  if (case %in% c("common", "shape1")) {
    fam <- get_family(family)
    target <- mean(Y <= tau1)
    target <- min(max(target, 0.05), 0.95)  # stay off the 0/1 boundary itself
    alpha0 <- tryCatch({
      root <- stats::uniroot(function(a) fam$cdf(tau1, c(lambda0, a)) - target,
                              lower = 1e-6, upper = 1e6, extendInt = "yes")$root
      if (is.finite(root) && root > 0) root else 1
    }, error = function(e) 1)
  }
  alpha1 <- 1
  par <- switch(case,
    common  = c(lambda0, lambda1, alpha0),
    shape1  = c(lambda0, alpha0, lambda1),
    shape11 = c(lambda0, lambda1)
  )
  log(pmax(par, 1e-3))
}

.fit_slt_once <- function(data, family, case, tau1, init_par, bootstrap, control) {
  # Same identifiability guard as fit_slt() (see the comment there): if
  # either stage has zero observations, the corresponding lambda is
  # completely unidentified (the objective is exactly flat in that
  # direction) and optim() would silently report back an arbitrary
  # "estimate". Unlike fit_slt(), this internal, per-replicate entry
  # point (used directly by slt_bootstrap_ci()/slt_gof_ks()) returns
  # NULL rather than erroring, so it is treated as an ordinary failed
  # replicate by the existing drop_failed_replicates()/B_converged
  # machinery in those callers -- exactly the outcome a degenerate
  # bootstrap/GOF replicate should have.
  if (length(data$Y) < 1L || length(data$Z) < 1L) return(NULL)
  ctrl_nm <- utils::modifyList(
    list(maxit = if (bootstrap) 300L else 2000L, reltol = 1e-8), control
  )
  fit_nm <- tryCatch(
    stats::optim(par = init_par, fn = neg_loglik_slt, method = "Nelder-Mead",
                 family = family, case = case, data = data, tau1 = tau1,
                 control = ctrl_nm),
    error = function(e) NULL
  )
  start2 <- if (!is.null(fit_nm)) fit_nm$par else init_par
  ctrl_bfgs <- utils::modifyList(
    list(maxit = if (bootstrap) 200L else 1000L,
         reltol = if (bootstrap) 1e-5 else 1e-7),
    control
  )
  fit <- tryCatch(
    stats::optim(par = start2, fn = neg_loglik_slt, method = "BFGS",
                 family = family, case = case, data = data, tau1 = tau1,
                 hessian = TRUE, control = ctrl_bfgs),
    error = function(e) NULL
  )
  if (is.null(fit) || fit$convergence != 0L || !is.finite(fit$value)) return(NULL)
  # Audit finding: optim() can report convergence = 0 ("successful")
  # while sitting exactly at (or extremely near) neg_loglik_slt()'s
  # penalty sentinel (.Machine$double.xmax / 2) -- e.g. found for
  # family = "gompertz" on real data where every parameter value tried
  # made a shifted elevated-stage observation numerically Inf, so every
  # point BFGS visited returned the same penalty and it "converged" by
  # having nowhere left to move, not by finding a real optimum. A
  # genuine negative log-likelihood, for any real dataset, is never
  # remotely close to this value (the threshold below is still ~150
  # orders of magnitude above any log-likelihood a real fit could
  # produce), so this is a safe, general guard against reporting that
  # kind of stuck optimiser run as a successful fit.
  if (fit$value >= .Machine$double.xmax / 4) return(NULL)
  fit
}

#' Fit a two-stage life testing model by maximum likelihood
#'
#' Fit the stage life testing model for given data using the
#' cumulative-exposure method for one of 12 baseline distributions and
#' one of three parameter cases.
#'
#' @details
#' Estimation uses a Nelder-Mead warm start followed by BFGS (with a
#' numerical Hessian for inference), repeated from several perturbed
#' starting values to guard against local optima; the converged fit
#' with the highest log-likelihood is returned.
#'
#' @param data A list with numeric vectors `Y` (stage-0 observed failure
#'   times) and `Z` (stage-1 observed failure times), as produced by
#'   [rslt()] (or assembled directly from real data).
#' @param family One of [slt_families()].
#' @param case One of `"common"`, `"shape1"`, `"shape11"` (see
#'   [slt_cases()]).
#' @param tau1 The stage-change time, used to shift `Z` back to the
#'   stage-1 origin.
#' @param init Optional numeric vector of starting values **on the
#'   natural scale**, in the order given by `case_info(case)$par_names`
#'   (see [slt_cases()] for the parameter order for each case). If `NULL`
#'   a data-driven default is used.
#' @param n_starts Number of (perturbed) starting values to try; the
#'   converged fit with the highest log-likelihood is returned.
#' @param bootstrap If `TRUE`, use looser/faster optimiser control
#'   settings appropriate for repeated bootstrap refits.
#' @param control Optional list overriding [stats::optim()] control
#'   parameters for both the Nelder-Mead and BFGS stages.
#' @return An object of class `"slt_fit"`, a list with elements including
#'   `converged`, `family`, `case`, `tau1`, `par` (natural-scale
#'   estimates), `se`, `theta0`, `theta1`, `vcov_log`, `loglik`, `aic`,
#'   `bic`. If no starting value converges, `converged` is `FALSE` and
#'   most other elements are absent.
#' @export
#' @examples
#' set.seed(1)
#' d <- rslt(150, "weibull", c(0.05, 1.5), c(0.2, 1.5), tau1 = 5, prop_transfer = 1)
#' fit <- fit_slt(d, family = "weibull", case = "shape1", tau1 = 5)
#' summary(fit)
fit_slt <- function(data, family = "weibull",
                       case = c("common", "shape1", "shape11"),
                       tau1, init = NULL, n_starts = 4L,
                       bootstrap = FALSE, control = list()) {
  case <- match.arg(case)
  fam <- get_family(family)
  cinfo <- case_info(case)
  if (missing(tau1) || !is.numeric(tau1) || length(tau1) != 1 || tau1 <= 0) {
    stop("tau1 must be a single positive number.", call. = FALSE)
  }
  if (is.null(data$Y) || is.null(data$Z)) {
    stop("'data' must be a list with numeric elements Y and Z.", call. = FALSE)
  }
  # Audit finding: an empty Y or Z previously passed the n_obs >=
  # cinfo$npar check below (which only counts the TOTAL across both
  # stages) and the optimiser would report converged = TRUE with a
  # specific-looking estimate for the affected stage's lambda -- but
  # that estimate is completely arbitrary. Every case has at least one
  # free parameter tied to each stage (a stage-0 lambda from Y, a
  # stage-1 lambda from Z), so if e.g. Z is empty the stage-1
  # log-likelihood contribution is identically 0 for *every* value of
  # lambda1 (verified: bit-for-bit identical neg_loglik_slt() output
  # across lambda1 spanning 4 orders of magnitude) -- a perfectly flat
  # objective in that direction, so optim() never moves away from
  # whatever the (arbitrary) starting value was and reports it back as
  # if it were a real estimate. Both stages must have at least one
  # observed failure for every parameter to be identified at all.
  if (length(data$Y) < 1L || length(data$Z) < 1L) {
    stop("'data$Y' and 'data$Z' must each contain at least one observed ",
         "failure time: every case has a stage-0 lambda (needs Y) and a ",
         "stage-1 lambda (needs Z), and the likelihood is exactly flat in ",
         "(hence cannot identify) a parameter whose stage has zero ",
         "observations.", call. = FALSE)
  }
  check_positive(data$Y, "data$Y")
  check_positive(data$Z, "data$Z")
  n_obs <- length(data$Y) + length(data$Z)
  if (n_obs < cinfo$npar) {
    stop(sprintf(
      "Not enough observations (%d) to estimate the %d free parameters of case '%s'.",
      n_obs, cinfo$npar, case
    ), call. = FALSE)
  }

  init_par <- if (!is.null(init)) {
    if (length(init) != cinfo$npar) {
      stop(sprintf("'init' must have length %d for case '%s' (%s).",
                    cinfo$npar, case, paste(cinfo$par_names, collapse = ", ")),
           call. = FALSE)
    }
    check_positive(init, "init")
    log(init)
  } else {
    .default_init(data, tau1, case, family)
  }

  n_starts <- max(1L, as.integer(n_starts))
  starts <- vector("list", n_starts)
  starts[[1]] <- init_par
  if (n_starts > 1L) {
    for (i in 2:n_starts) {
      starts[[i]] <- init_par + stats::rnorm(length(init_par), sd = 0.4 * (i - 1))
    }
  }

  fits <- lapply(starts, function(sp) {
    .fit_slt_once(data, fam$key, case, tau1, sp, bootstrap, control)
  })
  fits <- Filter(Negate(is.null), fits)

  if (!length(fits)) {
    return(structure(
      list(converged = FALSE, family = fam$key, family_name = fam$name,
           case = case, tau1 = tau1, data = data, call = match.call()),
      class = "slt_fit"
    ))
  }

  best <- fits[[which.min(vapply(fits, `[[`, numeric(1), "value"))]]
  log_par <- best$par
  natural_par <- exp(log_par)
  names(natural_par) <- cinfo$par_names
  pars <- slt_expand_params(natural_par, case)
  vc <- safe_vcov(best$hessian)
  se_log <- if (!is.null(vc)) sqrt(diag(vc)) else rep(NA_real_, cinfo$npar)
  se_natural <- natural_par * se_log
  names(se_natural) <- cinfo$par_names
  loglik <- -best$value
  n <- length(data$Y) + length(data$Z)

  structure(list(
    converged = TRUE,
    family = fam$key, family_name = fam$name, case = case,
    tau1 = tau1, data = data, n = n,
    log_par = log_par, par = natural_par, se = se_natural,
    theta0 = pars$theta0, theta1 = pars$theta1,
    vcov_log = vc, hessian = best$hessian,
    loglik = loglik, npar = cinfo$npar,
    aic = -2 * loglik + 2 * cinfo$npar,
    bic = -2 * loglik + log(n) * cinfo$npar,
    optim_fit = best,
    call = match.call()
  ), class = "slt_fit")
}
