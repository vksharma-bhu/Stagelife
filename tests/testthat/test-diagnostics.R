plot_fit_example <- function() {
  set.seed(41)
  theta0 <- c(0.08, 1.4); theta1 <- c(0.25, 1)
  tau1 <- slt_tau("weibull", theta0, 0.5)
  d <- rslt(200, "weibull", theta0, theta1, tau1, prop_transfer = 1)
  fit_slt(d, "weibull", "shape1", tau1 = tau1)
}

test_that("all diagnostic plot functions return ggplot objects, and plot.slt_fit dispatches correctly", {
  fit <- plot_fit_example()
  plots <- list(cdf = plot_cdf(fit), survival = plot_survival(fit), pp = plot_pp(fit), qq = plot_qq(fit))
  expect_true(all(vapply(plots, inherits, logical(1), "ggplot")))

  all_plots <- plot(fit, which = "all")
  expect_named(all_plots, c("cdf", "survival", "pp", "qq"))
  expect_true(all(vapply(all_plots, inherits, logical(1), "ggplot")))
})

test_that("plot functions error clearly on a non-converged fit", {
  fake <- structure(list(converged = FALSE), class = "slt_fit")
  expect_error(plot_cdf(fake), "did not converge")
})
