#' Kolmogorov-Smirnov goodness-of-fit test with parametric-bootstrap correction for a fitted two-stage life testing model
#'
#' Performs the Kolmogorov-Smirnov goodness-of-fit test and reports both
#' the naive p-value (for reference) and the standard
#' parametric-bootstrap-corrected p-value.
#'
#' @details
#' The classical KS null distribution assumes a fully specified (not
#' estimated) reference distribution; using it when the parameters were
#' estimated from the same data makes the naive p-value anti-conservative
#' (too small), a problem first noted by Lilliefors. The parametric
#' bootstrap correction simulates `B` fresh data sets of the same size
#' from the fitted model via [rslt()] (at the MLE), refits each with the
#' same family/case, and takes the corrected p-value to be the
#' proportion of bootstrap KS statistics at least as extreme as the one
#' observed on the real data.
#'
#' This test assumes every survivor was transferred at `tau1` (the
#' classical simple-step-stress design, `prop_transfer = 1` in [rslt()]),
#' because the single pooled reference distribution it uses ([pslt()]) is
#' only the correct marginal CDF of the *combined* `Y`/`Z` sample in that
#' case. Under genuine partial withdrawal, units that survive `tau1` but
#' are *not* transferred keep failing on stage 0 unmodified, which
#' [pslt()] does not represent -- pooling them into one KS statistic
#' against [pslt()] would be a biased test. If `fit$data` came from a
#' partial-withdrawal [rslt()] call (`prop_transfer < 1` or
#' `design = "typeM"`), this function is not appropriate; a Y-only and a
#' Z-only KS test against their own (design-invariant) marginals would be
#' the correct generalisation, but is not currently implemented.
#'
#' @param fit A converged `"slt_fit"` object, fit to data from a
#'   simple-step-stress ([rslt()] with `prop_transfer = 1`) design.
#' @param B Number of parametric bootstrap replicates.
#' @param n_cores Number of cores for [parallel::mclapply()].
#' @return An object of class `"slt_gof"` with elements `D`
#'   (KS statistic), `p_naive`, `p_bootstrap`, `B`, `B_converged`.
#' @export
#' @examples
#' \donttest{
#' set.seed(1)
#' d <- rslt(150, "weibull", c(0.05, 1.5), c(0.2, 1.5), tau1 = 5, prop_transfer = 1)
#' fit <- fit_slt(d, "weibull", "shape1", tau1 = 5)
#' slt_gof_ks(fit, B = 100)
#' }
slt_gof_ks <- function(fit, B = 1000, n_cores = 1L) {
  if (!isTRUE(fit$converged)) stop("'fit' did not converge.", call. = FALSE)
  Tobs <- c(fit$data$Y, fit$data$Z)
  cdf_fn <- function(t) pslt(t, fit$family, fit$theta0, fit$theta1, fit$tau1)
  ks_obs <- suppressWarnings(stats::ks.test(Tobs, cdf_fn))
  D_obs <- unname(ks_obs$statistic)

  one_rep <- function(b) {
    tryCatch({
      bd <- rslt(fit$n, fit$family, fit$theta0, fit$theta1, fit$tau1,
                 design = "typeP", prop_transfer = 1)
      fb <- .fit_slt_once(bd, fit$family, fit$case, fit$tau1, fit$log_par,
                             bootstrap = TRUE, control = list())
      if (is.null(fb)) return(NA_real_)
      natural <- exp(fb$par)
      pars <- slt_expand_params(natural, fit$case)
      cdf_b <- function(t) pslt(t, fit$family, pars$theta0, pars$theta1, fit$tau1)
      Tb <- c(bd$Y, bd$Z)
      unname(suppressWarnings(stats::ks.test(Tb, cdf_b))$statistic)
    }, error = function(e) NA_real_)
  }
  D_boot_list <- if (n_cores > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(seq_len(B), one_rep, mc.cores = n_cores)
  } else {
    lapply(seq_len(B), one_rep)
  }
  D_boot <- suppressWarnings(as.numeric(unlist(drop_failed_replicates(D_boot_list))))
  D_boot <- D_boot[is.finite(D_boot)]

  structure(list(
    D = D_obs, p_naive = unname(ks_obs$p.value),
    p_bootstrap = mean(D_boot >= D_obs), B = B, B_converged = length(D_boot)
  ), class = "slt_gof")
}

#' @export
print.slt_gof <- function(x, digits = 4, ...) {
  cat("Kolmogorov-Smirnov goodness-of-fit (two-stage cumulative-exposure model)\n")
  cat(sprintf("  D statistic         : %s\n", round(x$D, digits)))
  cat(sprintf("  p-value (naive)     : %s   [anti-conservative: parameters estimated]\n",
              round(x$p_naive, digits)))
  cat(sprintf("  p-value (bootstrap) : %s   (%d/%d replicates converged)\n",
              round(x$p_bootstrap, digits), x$B_converged, x$B))
  invisible(x)
}

#' Likelihood-ratio test between nested parameter-sharing cases
#'
#' Tests whether a reduced (nested) parameter-sharing case is adequate
#' against a less-restricted alternative, via the likelihood-ratio
#' statistic.
#'
#' @details
#' The only nested pairs among the three supported cases are
#' `("shape1", "shape11")` and `("common", "shape11")` (both fix the
#' remaining free shape parameter(s) to 1); see [slt_cases()]. This
#' function tests `H0`: the reduced (`"shape11"`) model is adequate,
#' against the corresponding less-restricted alternative.
#'
#' @param fit_full,fit_reduced Converged `"slt_fit"` objects for the same
#'   `family`, `tau1` and `data`, with `fit_reduced` nested in `fit_full`.
#' @return An object of class `"slt_lrt"` with `statistic`, `df`,
#'   `p.value`.
#' @export
#' @examples
#' set.seed(1)
#' d <- rslt(200, "weibull", c(0.05, 1), c(0.2, 1), tau1 = 5, prop_transfer = 1)
#' f1 <- fit_slt(d, "weibull", "shape1", tau1 = 5)
#' f0 <- fit_slt(d, "weibull", "shape11", tau1 = 5)
#' slt_lrt(f1, f0)
slt_lrt <- function(fit_full, fit_reduced) {
  if (!isTRUE(fit_full$converged) || !isTRUE(fit_reduced$converged)) {
    stop("Both models must have converged.", call. = FALSE)
  }
  if (fit_full$family != fit_reduced$family) {
    stop("'fit_full' and 'fit_reduced' must use the same family.", call. = FALSE)
  }
  valid <- (fit_full$case == "shape1" && fit_reduced$case == "shape11") ||
    (fit_full$case == "common" && fit_reduced$case == "shape11")
  if (!valid) {
    stop("slt_lrt() supports (fit_full$case, fit_reduced$case) = ",
         "(\"shape1\",\"shape11\") or (\"common\",\"shape11\") only; ",
         "these are the only nested pairs among the supported cases.",
         call. = FALSE)
  }
  df <- fit_full$npar - fit_reduced$npar
  LR <- max(2 * (fit_full$loglik - fit_reduced$loglik), 0)
  p_value <- stats::pchisq(LR, df = df, lower.tail = FALSE)
  structure(list(statistic = LR, df = df, p.value = p_value,
                 model_full = fit_full$case, model_reduced = fit_reduced$case,
                 family = fit_full$family_name),
            class = "slt_lrt")
}

#' @export
print.slt_lrt <- function(x, digits = 4, ...) {
  cat(sprintf("Likelihood-ratio test (%s): H1 = '%s' vs. H0 = '%s'\n",
              x$family, x$model_full, x$model_reduced))
  cat(sprintf("  LR statistic = %s, df = %d, p-value = %s\n",
              round(x$statistic, digits), x$df, format.pval(x$p.value, digits = digits)))
  invisible(x)
}

#' Compare several fitted models by information criteria
#'
#' Builds a ranked comparison table across any number of fitted models,
#' by information criteria.
#'
#' @details
#' By default sorted by AIC, across any number of `"slt_fit"` objects --
#' typically different `family`/`case` combinations fit to the same
#' data set -- for exploratory model selection on real data. This is a
#' lightweight diagnostic table for a single data set, not a Monte
#' Carlo simulation study.
#'
#' @param fits A (preferably named) list of `"slt_fit"` objects.
#' @param sort_by One of `"aic"` (default), `"bic"`, `"loglik"`.
#' @return A data.frame with columns `model`, `family`, `case`, `n_par`,
#'   `logLik`, `AIC`, `BIC`, `converged`.
#' @export
#' @examples
#' set.seed(1)
#' d <- rslt(200, "weibull", c(0.05, 1.5), c(0.2, 1.5), tau1 = 5, prop_transfer = 1)
#' fits <- list(
#'   "Weibull / shape1"  = fit_slt(d, "weibull", "shape1", tau1 = 5),
#'   "Weibull / shape11" = fit_slt(d, "weibull", "shape11", tau1 = 5),
#'   "Weibull / common"  = fit_slt(d, "weibull", "common", tau1 = 5)
#' )
#' slt_model_comparison(fits)
slt_model_comparison <- function(fits, sort_by = c("aic", "bic", "loglik")) {
  sort_by <- match.arg(sort_by)
  if (is.null(names(fits)) || any(names(fits) == "")) {
    names(fits) <- paste0("model", seq_along(fits))
  }
  rows <- lapply(names(fits), function(nm) {
    f <- fits[[nm]]
    conv <- isTRUE(f$converged)
    data.frame(
      model = nm, family = f$family_name %||% f$family, case = f$case,
      n_par = if (conv) f$npar else NA_integer_,
      logLik = if (conv) f$loglik else NA_real_,
      AIC = if (conv) f$aic else NA_real_,
      BIC = if (conv) f$bic else NA_real_,
      converged = conv, stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  col <- switch(sort_by, aic = "AIC", bic = "BIC", loglik = "logLik")
  out <- if (col == "logLik") out[order(-out[[col]]), ] else out[order(out[[col]]), ]
  rownames(out) <- NULL
  round_df(out, 4)
}
