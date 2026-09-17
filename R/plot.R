.slt_theme <- function() {
  ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}

.plot_cdf_data <- function(fit, npoints = 300) {
  Tobs <- sort(c(fit$data$Y, fit$data$Z))
  grid <- seq(1e-6, max(Tobs) * 1.05, length.out = npoints)
  fitted <- pslt(grid, fit$family, fit$theta0, fit$theta1, fit$tau1)
  list(Tobs = Tobs, grid = grid, fitted = fitted)
}

#' Fitted vs empirical CDF plot
#'
#' Plots the fitted population CDF ([pslt()]) against the empirical CDF
#' of the observed data.
#'
#' @details
#' Assumes every survivor was transferred at `tau1` (the classical
#' simple-step-stress design, `prop_transfer = 1` in [rslt()]) -- see the
#' note in [slt_gof_ks()] for why: the pooled `Y`/`Z` sample only has
#' [pslt()] as its correct marginal CDF in that case. This applies to all
#' four diagnostic plots in this file.
#'
#' @param fit A converged `"slt_fit"` object, fit to data from a
#'   simple-step-stress ([rslt()] with `prop_transfer = 1`) design.
#' @return A `ggplot` object.
#' @export
#' @examples
#' set.seed(1)
#' d <- rslt(150, "weibull", c(0.05, 1.5), c(0.2, 1.5), tau1 = 5, prop_transfer = 1)
#' fit <- fit_slt(d, "weibull", "shape1", tau1 = 5)
#' plot_cdf(fit)
plot_cdf <- function(fit) {
  if (!isTRUE(fit$converged)) stop("'fit' did not converge.", call. = FALSE)
  dd <- .plot_cdf_data(fit)
  df_fit <- data.frame(t = dd$grid, F = dd$fitted)
  ggplot2::ggplot() +
    ggplot2::stat_ecdf(data = data.frame(t = dd$Tobs), ggplot2::aes(x = .data$t),
                        geom = "point", pad = FALSE, alpha = 0.6) +
    ggplot2::geom_line(data = df_fit, ggplot2::aes(x = .data$t, y = .data$F), colour = "steelblue") +
    ggplot2::geom_vline(xintercept = fit$tau1, linetype = "dashed", colour = "grey40") +
    ggplot2::labs(x = "Time", y = "F(t)",
                  title = sprintf("Fitted vs. empirical CDF (%s, %s)", fit$family_name, fit$case)) +
    .slt_theme()
}

#' Fitted vs empirical survival curve
#' @inheritParams plot_cdf
#' @return A `ggplot` object.
#' @export
#' @examples
#' set.seed(1)
#' d <- rslt(150, "weibull", c(0.05, 1.5), c(0.2, 1.5), tau1 = 5, prop_transfer = 1)
#' fit <- fit_slt(d, "weibull", "shape1", tau1 = 5)
#' plot_survival(fit)
plot_survival <- function(fit) {
  if (!isTRUE(fit$converged)) stop("'fit' did not converge.", call. = FALSE)
  dd <- .plot_cdf_data(fit)
  df_fit <- data.frame(t = dd$grid, S = 1 - dd$fitted)
  ggplot2::ggplot() +
    ggplot2::stat_ecdf(data = data.frame(t = dd$Tobs), ggplot2::aes(x = .data$t, y = 1 - ggplot2::after_stat(.data$y)),
                        geom = "point", pad = FALSE, alpha = 0.6) +
    ggplot2::geom_line(data = df_fit, ggplot2::aes(x = .data$t, y = .data$S), colour = "firebrick") +
    ggplot2::geom_vline(xintercept = fit$tau1, linetype = "dashed", colour = "grey40") +
    ggplot2::labs(x = "Time", y = "S(t) = 1 - F(t)",
                  title = sprintf("Fitted vs. empirical survival (%s, %s)", fit$family_name, fit$case)) +
    .slt_theme()
}

#' P-P (probability-probability) plot
#' @inheritParams plot_cdf
#' @return A `ggplot` object.
#' @export
#' @examples
#' set.seed(1)
#' d <- rslt(150, "weibull", c(0.05, 1.5), c(0.2, 1.5), tau1 = 5, prop_transfer = 1)
#' fit <- fit_slt(d, "weibull", "shape1", tau1 = 5)
#' plot_pp(fit)
plot_pp <- function(fit) {
  if (!isTRUE(fit$converged)) stop("'fit' did not converge.", call. = FALSE)
  Tobs <- sort(c(fit$data$Y, fit$data$Z))
  n <- length(Tobs)
  emp <- stats::ppoints(n)
  fitted_p <- pslt(Tobs, fit$family, fit$theta0, fit$theta1, fit$tau1)
  df <- data.frame(empirical = emp, fitted = fitted_p)
  ggplot2::ggplot(df, ggplot2::aes(x = .data$fitted, y = .data$empirical)) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "steelblue", linetype = "dashed") +
    ggplot2::labs(x = "Fitted F(t)", y = "Empirical probability",
                  title = sprintf("P-P plot (%s, %s)", fit$family_name, fit$case)) +
    .slt_theme()
}

#' Q-Q (quantile-quantile) plot
#' @inheritParams plot_cdf
#' @return A `ggplot` object.
#' @export
#' @examples
#' set.seed(1)
#' d <- rslt(150, "weibull", c(0.05, 1.5), c(0.2, 1.5), tau1 = 5, prop_transfer = 1)
#' fit <- fit_slt(d, "weibull", "shape1", tau1 = 5)
#' plot_qq(fit)
plot_qq <- function(fit) {
  if (!isTRUE(fit$converged)) stop("'fit' did not converge.", call. = FALSE)
  Tobs <- sort(c(fit$data$Y, fit$data$Z))
  n <- length(Tobs)
  emp <- stats::ppoints(n)
  fitted_q <- qslt(emp, fit$family, fit$theta0, fit$theta1, fit$tau1)
  df <- data.frame(empirical = Tobs, fitted = fitted_q)
  rng <- range(c(df$empirical, df$fitted), finite = TRUE)
  ggplot2::ggplot(df, ggplot2::aes(x = .data$fitted, y = .data$empirical)) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "steelblue", linetype = "dashed") +
    ggplot2::coord_cartesian(xlim = rng, ylim = rng) +
    ggplot2::labs(x = "Fitted quantile", y = "Empirical (observed) quantile",
                  title = sprintf("Q-Q plot (%s, %s)", fit$family_name, fit$case)) +
    .slt_theme()
}

#' Diagnostic plots for a fitted two-stage life-testing model
#'
#' See [plot_cdf()] for an important caveat: these assume the classical
#' simple-step-stress design (`prop_transfer = 1` in [rslt()]).
#'
#' @param x A converged `"slt_fit"` object.
#' @param which One of `"cdf"`, `"survival"`, `"pp"`, `"qq"`, or `"all"`
#'   (returns a named list of all four `ggplot` objects).
#' @param ... Currently unused.
#' @return A `ggplot` object, or (if `which = "all"`) a named list of four
#'   `ggplot` objects.
#' @export
#' @examples
#' set.seed(1)
#' d <- rslt(150, "weibull", c(0.05, 1.5), c(0.2, 1.5), tau1 = 5, prop_transfer = 1)
#' fit <- fit_slt(d, "weibull", "shape1", tau1 = 5)
#' plot(fit, which = "cdf")
plot.slt_fit <- function(x, which = c("cdf", "survival", "pp", "qq", "all"), ...) {
  which <- match.arg(which)
  if (which == "all") {
    return(list(cdf = plot_cdf(x), survival = plot_survival(x),
                pp = plot_pp(x), qq = plot_qq(x)))
  }
  switch(which,
    cdf = plot_cdf(x), survival = plot_survival(x),
    pp = plot_pp(x), qq = plot_qq(x)
  )
}
