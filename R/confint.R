#' Asymptotic (Wald) confidence intervals for a fitted two-stage life testing model
#'
#' Two Wald-type intervals are available: `"normal"` and `"lognormal"`.
#'
#' @details
#' `"lognormal"` (the default) builds the interval on the `log(theta)`
#' scale (where the MLE is asymptotically normal by construction) and
#' exponentiates, which guarantees a positive interval and is generally
#' preferred for these positive shape/scale parameters. `"normal"`
#' applies the delta-method standard error directly on the natural
#' scale and can (rarely) produce a negative lower limit for small
#' samples.
#'
#' @param object An `"slt_fit"` object from [fit_slt()].
#' @param parm Optional character vector of parameter names to include
#'   (default: all).
#' @param level Confidence level.
#' @param type `"lognormal"` (default) or `"normal"`.
#' @param ... Currently unused.
#' @return A numeric matrix with one row per parameter and two columns
#'   (lower, upper).
#' @export
#' @examples
#' set.seed(1)
#' d <- rslt(200, "weibull", c(0.05, 1.5), c(0.2, 1.5), tau1 = 5, prop_transfer = 1)
#' fit <- fit_slt(d, "weibull", "shape1", tau1 = 5)
#' confint(fit)
confint.slt_fit <- function(object, parm, level = 0.95,
                             type = c("lognormal", "normal"), ...) {
  type <- match.arg(type)
  check_level(level)
  if (!isTRUE(object$converged)) {
    stop("Model did not converge; no confidence intervals available.", call. = FALSE)
  }
  if (is.null(object$vcov_log)) {
    stop("Hessian was not invertible at the MLE; asymptotic CIs unavailable. ",
         "Try slt_bootstrap_ci().", call. = FALSE)
  }
  alpha <- 1 - level
  z <- stats::qnorm(1 - alpha / 2)
  se_log <- sqrt(diag(object$vcov_log))
  if (type == "lognormal") {
    lower <- exp(object$log_par - z * se_log)
    upper <- exp(object$log_par + z * se_log)
  } else {
    lower <- object$par - z * object$se
    upper <- object$par + z * object$se
  }
  out <- cbind(lower, upper)
  rownames(out) <- names(object$par)
  colnames(out) <- c(sprintf("%.1f%%", 100 * alpha / 2),
                      sprintf("%.1f%%", 100 * (1 - alpha / 2)))
  if (!missing(parm) && !is.null(parm)) out <- out[parm, , drop = FALSE]
  out
}

#' Bootstrap confidence intervals for a fitted two-stage life testing model
#'
#' Computes percentile, bootstrap-t, and log-scale bootstrap-t confidence
#' intervals via either non-parametric or a parametric bootstrap.
#'
#' @details
#' **Non-parametric** bootstrap resamples (case resampling) the observed
#' `Y`/`Z` failure times; **parametric** bootstrap simulates fresh data
#' from the fitted model at the MLE. The parametric bootstrap is
#' generally more efficient when the fitted model is (approximately)
#' correctly specified -- see [slt_gof_ks()] -- but is more sensitive to
#' misspecification; the non-parametric bootstrap is robust to
#' misspecification but, for small samples, can only reproduce
#' combinations of already-observed values.
#'
#' The plain (natural-scale) bootstrap-t interval, like the asymptotic
#' Wald interval from `confint(..., type = "normal")`, is an ordinary
#' linear pivot (`estimate -/+ quantile * se`) and can occasionally
#' produce a negative lower bound for a positive parameter, especially
#' at small `n` or for a skewed bootstrap distribution. The log-scale
#' bootstrap-t interval applies the same pivot construction on
#' `log(theta)` and exponentiates the result, which guarantees a
#' positive interval by construction -- directly analogous to why
#' `confint(..., type = "lognormal")` (the default) is generally
#' preferred over `type = "normal"` for the asymptotic interval. The
#' percentile interval is always positive already (it is built from
#' quantiles of the always-positive bootstrap estimates themselves), so
#' only bootstrap-t needed a log-scale counterpart.
#'
#' @param fit A converged `"slt_fit"` object.
#' @param B Number of bootstrap replicates.
#' @param level Confidence level.
#' @param method `"nonparametric"` (default) or `"parametric"`.
#' @param design For `method = "parametric"`, the design used to simulate
#'   new data via [rslt()]: `"typeP"` (fixed target proportion, the
#'   default) or `"typeM"` (fixed target count).
#' @param prop_transfer Proportion transferred, for `design = "typeP"`.
#'   Defaults to `1` (every survivor transferred, the classical
#'   simple-step-stress design, matching a typical published step-stress
#'   data set); set to the actual value used in your experiment if it was
#'   a genuine partial-withdrawal design.
#' @param R0_1 Target number transferred, for `design = "typeM"`.
#' @param n_cores Number of cores for [parallel::mclapply()] (ignored on
#'   Windows, where replicates run sequentially).
#' @return An object of class `"slt_boot_ci"`: a list with `percentile`,
#'   `boot_t`, and `log_boot_t` data.frames (columns `lower`, `upper`,
#'   one row per parameter), `B`, `B_converged`, `method`, `design`,
#'   `level`, and the raw `boot_estimates` matrix.
#' @export
#' @examples
#' \donttest{
#' set.seed(1)
#' d <- rslt(150, "weibull", c(0.05, 1.5), c(0.2, 1.5), tau1 = 5, prop_transfer = 1)
#' fit <- fit_slt(d, "weibull", "shape1", tau1 = 5)
#' slt_bootstrap_ci(fit, B = 100)
#' }
slt_bootstrap_ci <- function(fit, B = 1000, level = 0.95,
                          method = c("nonparametric", "parametric"),
                          design = c("typeP", "typeM"),
                          prop_transfer = 1, R0_1 = NULL,
                          n_cores = 1L) {
  method <- match.arg(method); design <- match.arg(design)
  check_level(level)
  if (!isTRUE(fit$converged)) stop("'fit' did not converge.", call. = FALSE)
  if (is.null(fit$vcov_log)) {
    stop("Hessian was not invertible at the original fit; bootstrap-t requires it.",
         call. = FALSE)
  }
  alpha <- 1 - level
  cinfo <- case_info(fit$case)
  orig_log <- fit$log_par
  orig_est <- fit$par
  orig_se <- fit$se
  orig_se_log <- sqrt(diag(fit$vcov_log))

  gen_one <- function(b) {
    tryCatch({
      bd <- if (method == "nonparametric") {
        list(Y = if (length(fit$data$Y)) sample(fit$data$Y, length(fit$data$Y), replace = TRUE) else numeric(0),
             Z = if (length(fit$data$Z)) sample(fit$data$Z, length(fit$data$Z), replace = TRUE) else numeric(0))
      } else {
        rslt(fit$n, fit$family, fit$theta0, fit$theta1, tau1 = fit$tau1,
             design = design, prop_transfer = prop_transfer, R0_1 = R0_1)
      }
      fb <- .fit_slt_once(bd, fit$family, fit$case, fit$tau1, orig_log,
                             bootstrap = TRUE, control = list())
      if (is.null(fb)) return(NULL)
      vcb <- safe_vcov(fb$hessian)
      if (is.null(vcb)) return(NULL)
      est_b <- exp(fb$par); names(est_b) <- cinfo$par_names
      se_log_b <- sqrt(diag(vcb))
      list(est = est_b, se = est_b * se_log_b, log_par = fb$par, se_log = se_log_b)
    }, error = function(e) NULL)
  }

  res <- if (n_cores > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(seq_len(B), gen_one, mc.cores = n_cores)
  } else {
    lapply(seq_len(B), gen_one)
  }
  res <- drop_failed_replicates(res)
  n_conv <- length(res)
  min_required <- min(B, 30L)  # can't require more successes than replicates requested
  if (n_conv < min_required) {
    stop(sprintf("Only %d/%d bootstrap replicates converged without a numerical error; increase B. ",
                  n_conv, B),
         "This is more common for family = 'powerlindley' at extreme parameter ",
         "values -- see vignette(\"stagelife-intro\").", call. = FALSE)
  }

  boot_est <- do.call(rbind, lapply(res, `[[`, "est"))
  boot_se  <- do.call(rbind, lapply(res, `[[`, "se"))
  boot_log_par <- do.call(rbind, lapply(res, `[[`, "log_par"))
  boot_se_log  <- do.call(rbind, lapply(res, `[[`, "se_log"))

  perc_lower <- apply(boot_est, 2, stats::quantile, probs = alpha / 2, na.rm = TRUE)
  perc_upper <- apply(boot_est, 2, stats::quantile, probs = 1 - alpha / 2, na.rm = TRUE)

  t_star  <- sweep(boot_est, 2, orig_est, "-") / boot_se
  tq_low  <- apply(t_star, 2, stats::quantile, probs = alpha / 2, na.rm = TRUE)
  tq_high <- apply(t_star, 2, stats::quantile, probs = 1 - alpha / 2, na.rm = TRUE)
  bt_lower <- orig_est - tq_high * orig_se
  bt_upper <- orig_est - tq_low * orig_se

  # Log-scale bootstrap-t: the same pivot construction as bt_lower/
  # bt_upper above, but built from log(theta) (where boot_log_par,
  # boot_se_log already live) and exponentiated at the end. This
  # guarantees a positive interval by construction, addressing the same
  # small-sample "negative lower bound" issue that motivates
  # confint(..., type = "lognormal") for the asymptotic Wald interval --
  # the plain (natural-scale) bootstrap-t above has no such guarantee,
  # since orig_est - tq_high*orig_se is an ordinary linear combination
  # that can go negative for a skewed/small-n bootstrap distribution.
  t_star_log  <- sweep(boot_log_par, 2, orig_log, "-") / boot_se_log
  tq_low_log  <- apply(t_star_log, 2, stats::quantile, probs = alpha / 2, na.rm = TRUE)
  tq_high_log <- apply(t_star_log, 2, stats::quantile, probs = 1 - alpha / 2, na.rm = TRUE)
  lbt_lower <- exp(orig_log - tq_high_log * orig_se_log)
  lbt_upper <- exp(orig_log - tq_low_log * orig_se_log)

  structure(list(
    method = method, design = if (method == "parametric") design else NA_character_,
    B = B, B_converged = n_conv, level = level,
    percentile = data.frame(lower = perc_lower, upper = perc_upper,
                             row.names = names(orig_est)),
    boot_t = data.frame(lower = bt_lower, upper = bt_upper,
                         row.names = names(orig_est)),
    log_boot_t = data.frame(lower = lbt_lower, upper = lbt_upper,
                             row.names = names(orig_est)),
    boot_estimates = boot_est
  ), class = "slt_boot_ci")
}

#' @export
print.slt_boot_ci <- function(x, digits = 4, ...) {
  cat(sprintf("Bootstrap confidence intervals (%s%s, B = %d, %d converged, %.0f%% level)\n",
              x$method, if (!is.na(x$design)) paste0(", design = ", x$design) else "",
              x$B, x$B_converged, 100 * x$level))
  cat("\nPercentile:\n"); print(round_df(x$percentile, digits))
  cat("\nBootstrap-t:\n"); print(round_df(x$boot_t, digits))
  cat("\nLog-scale bootstrap-t (always positive):\n"); print(round_df(x$log_boot_t, digits))
  invisible(x)
}
