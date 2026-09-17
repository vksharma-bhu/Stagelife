############################################################################
# Likelihood for the stage life testing (SLT) model (rslt()): Y observed
# under the theta0 density; Z observed under the theta1 density after
# being shifted back, via the cumulative-exposure equivalent age, to the
# theta1 time origin.
############################################################################

#' Two-stage cumulative-exposure log-likelihood
#'
#' @inheritParams slt_equiv_age
#' @param data A list with numeric vectors `Y` (stage-0 observed failure
#'   times) and `Z` (stage-1 observed failure times), as returned by
#'   [rslt()].
#' @param tau1 The stage-change time used to shift `Z` back to the
#'   stage-1 time origin.
#' @return The log-likelihood (a scalar), or `-Inf` for parameter values
#'   outside the support of the model.
#' @export
#' @examples
#' set.seed(1)
#' d <- rslt(80, "weibull", c(0.05, 1.5), c(0.2, 1.5), tau1 = 5, prop_transfer = 1)
#' loglik_slt(d, "weibull", c(0.05, 1.5), c(0.2, 1.5), tau1 = 5)
loglik_slt <- function(data, family, theta0, theta1, tau1) {
  fam <- get_family(family)
  Y <- data$Y; Z <- data$Z
  nu1 <- fam$quantile(fam$cdf(tau1, theta0), theta1)
  Z1 <- Z + nu1 - tau1
  if (any(Y <= 0) || (length(Z1) && any(Z1 <= 0))) return(-Inf)
  lf0 <- log(fam$pdf(Y, theta0))
  if (any(!is.finite(lf0))) return(-Inf)
  lf1 <- if (length(Z1)) {
    v <- log(fam$pdf(Z1, theta1))
    if (any(!is.finite(v))) return(-Inf)
    v
  } else numeric(0L)
  sum(lf0) + sum(lf1)
}
loglik_slt <- compiler::cmpfun(loglik_slt)

#' Case-aware negative log-likelihood (optimiser objective)
#'
#' Maps a log-scale, case-specific parameter vector to `(theta0, theta1)`
#' via `slt_expand_params()` and returns the negative log-likelihood, with any
#' internal error (e.g. from the Lambert-W-based Power Lindley quantile
#' straying outside `[0, 1]` during a finite-difference gradient step) or
#' non-finite result converted to a large, finite penalty so that
#' [stats::optim()] never aborts.
#'
#' @param log_par Numeric vector, `log()` of the natural-scale case
#'   parameters (ensures positivity without constrained optimisation).
#' @inheritParams loglik_slt
#' @param case One of `"common"`, `"shape1"`, `"shape11"`.
#' @return A finite scalar (never `-Inf`/`Inf`/`NaN`), suitable as an
#'   [stats::optim()] objective.
#' @keywords internal
neg_loglik_slt <- function(log_par, family, case, data, tau1) {
  ll <- tryCatch({
    pars <- slt_expand_params(exp(log_par), case)
    loglik_slt(data, family, pars$theta0, pars$theta1, tau1)
  }, error = function(e) -Inf)
  if (!is.finite(ll)) return(.Machine$double.xmax / 2)
  -ll
}
neg_loglik_slt <- compiler::cmpfun(neg_loglik_slt)
