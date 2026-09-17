#' Print a fitted two-stage life testing model
#'
#' @param x An `"slt_fit"` object from [fit_slt()].
#' @param digits Number of digits to round to.
#' @param ... Currently unused.
#' @return `x`, invisibly.
#' @seealso [summary.slt_fit()] for a more detailed summary with
#'   confidence limits.
#' @export
print.slt_fit <- function(x, digits = 4, ...) {
  cat("Two-stage lifetime model fit ('slt_fit')\n")
  cat(sprintf("  Family : %s\n", x$family_name %||% x$family))
  cat(sprintf("  Case   : %s\n", x$case))
  cat(sprintf("  tau1   : %s\n", format(x$tau1)))
  if (!isTRUE(x$converged)) {
    cat("  Status : DID NOT CONVERGE (all starting values failed)\n")
    return(invisible(x))
  }
  cat(sprintf("  n      : %d  (n0 = %d, n1 = %d)\n",
              x$n, length(x$data$Y), length(x$data$Z)))
  cat(sprintf("  logLik : %s   AIC: %s   BIC: %s\n",
              format(round(x$loglik, digits)), format(round(x$aic, digits)),
              format(round(x$bic, digits))))
  cat("\nCoefficients (natural scale):\n")
  tab <- data.frame(Estimate = x$par, `Std. Error` = x$se, check.names = FALSE)
  print(round_df(tab, digits))
  invisible(x)
}

#' Summarise a fitted two-stage life testing model
#'
#' @param object An `"slt_fit"` object from [fit_slt()].
#' @param level Confidence level for the Wald (log-normal) confidence
#'   limits reported in the coefficient table.
#' @param ... Currently unused.
#' @return An object of class `"summary.slt_fit"`.
#' @export
summary.slt_fit <- function(object, level = 0.95, ...) {
  if (!isTRUE(object$converged)) {
    return(structure(list(converged = FALSE, family = object$family,
                           case = object$case), class = "summary.slt_fit"))
  }
  z <- stats::qnorm(1 - (1 - level) / 2)
  se_log <- if (!is.null(object$vcov_log)) sqrt(diag(object$vcov_log)) else NA
  lower <- exp(object$log_par - z * se_log)
  upper <- exp(object$log_par + z * se_log)
  tab <- data.frame(
    Estimate = object$par, `Std. Error` = object$se,
    Lower = lower, Upper = upper, check.names = FALSE
  )
  rownames(tab) <- names(object$par)
  structure(list(
    converged = TRUE, family = object$family_name, case = object$case,
    tau1 = object$tau1, n = object$n, loglik = object$loglik,
    aic = object$aic, bic = object$bic, level = level, coefficients = tab,
    theta0 = object$theta0, theta1 = object$theta1
  ), class = "summary.slt_fit")
}

#' Print a summary of a fitted two-stage life testing model
#'
#' @param x A `"summary.slt_fit"` object from [summary.slt_fit()].
#' @param digits Number of digits to round to.
#' @param ... Currently unused.
#' @return `x`, invisibly.
#' @export
print.summary.slt_fit <- function(x, digits = 4, ...) {
  cat("Summary of two-stage lifetime model fit\n")
  cat(sprintf("  Family : %s\n  Case   : %s\n", x$family, x$case))
  if (!isTRUE(x$converged)) {
    cat("  Status : DID NOT CONVERGE\n")
    return(invisible(x))
  }
  cat(sprintf("  n = %d   logLik = %s   AIC = %s   BIC = %s\n",
              x$n, round(x$loglik, digits), round(x$aic, digits), round(x$bic, digits)))
  cat(sprintf("\nCoefficients (natural scale, %.0f%% Wald log-normal CI):\n", 100 * x$level))
  print(round_df(x$coefficients, digits))
  cat(sprintf("\ntheta0 (stage 0) = (%s)\ntheta1 (stage 1) = (%s)\n",
              paste(round(x$theta0, digits), collapse = ", "),
              paste(round(x$theta1, digits), collapse = ", ")))
  invisible(x)
}

#' Extract coefficients from a fitted two-stage life testing model
#'
#' @param object An `"slt_fit"` object from [fit_slt()].
#' @param ... Currently unused.
#' @return The named numeric vector `object$par` (natural, positive
#'   scale; names and length depend on `object$case`, see
#'   [slt_cases()]).
#' @export
coef.slt_fit <- function(object, ...) object$par

#' Variance-covariance matrix of a fitted two-stage life testing model's coefficients
#'
#' Computes the natural-scale variance-covariance matrix by the delta
#' method, from the log-scale one fitted internally.
#'
#' @details
#' `object$vcov_log` (the inverse of the numerical Hessian of the
#' negative log-likelihood in `log(par)`) is transformed to the natural
#' scale by the delta method: with \eqn{\eta = \log(\theta)},
#' \eqn{d\theta/d\eta = \mathrm{diag}(\theta)}, so
#' \eqn{\mathrm{Var}(\theta) \approx \mathrm{diag}(\theta)\,
#' \mathrm{Var}(\eta)\,\mathrm{diag}(\theta)}.
#'
#' @param object An `"slt_fit"` object from [fit_slt()].
#' @param ... Currently unused.
#' @return The natural-scale variance-covariance matrix, or `NULL` if
#'   the Hessian was not invertible at the optimum.
#' @seealso [confint.slt_fit()], which uses this (or the log scale
#'   directly) to build a confidence interval.
#' @export
vcov.slt_fit <- function(object, ...) {
  if (is.null(object$vcov_log)) return(NULL)
  # delta method: Var(exp(eta)) approx diag(exp(eta)) %*% Var(eta) %*% diag(exp(eta))
  J <- diag(object$par, nrow = length(object$par))
  vc <- J %*% object$vcov_log %*% J
  dimnames(vc) <- list(names(object$par), names(object$par))
  vc
}

#' Extract the log-likelihood of a fitted two-stage life testing model
#'
#' @param object An `"slt_fit"` object from [fit_slt()].
#' @param ... Currently unused.
#' @return An object of class `"logLik"` (the fitted log-likelihood,
#'   `NA_real_` if `object` did not converge, with `"df"`/`"nobs"`
#'   attributes set for use by [stats::AIC()]/[stats::BIC()]).
#' @export
logLik.slt_fit <- function(object, ...) {
  val <- if (isTRUE(object$converged)) object$loglik else NA_real_
  structure(val, df = object$npar, nobs = object$n, class = "logLik")
}
