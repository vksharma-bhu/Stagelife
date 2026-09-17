setup_fits <- function(seed = 31, n = 400) {
  set.seed(seed)
  theta0 <- c(0.08, 1); theta1 <- c(0.25, 1)  # shape11-true data
  tau1 <- slt_tau("weibull", theta0, 0.5)
  d <- rslt(n, "weibull", theta0, theta1, tau1, prop_transfer = 1)
  list(
    data = d, tau1 = tau1,
    shape1  = fit_slt(d, "weibull", "shape1", tau1 = tau1),
    shape11 = fit_slt(d, "weibull", "shape11", tau1 = tau1),
    common  = fit_slt(d, "weibull", "common", tau1 = tau1)
  )
}

test_that("slt_gof_ks returns a valid D statistic and both p-values in [0,1]", {
  skip_on_cran()
  fits <- setup_fits()
  set.seed(200)
  g <- slt_gof_ks(fits$shape1, B = 60)
  expect_s3_class(g, "slt_gof")
  expect_true(g$D >= 0 && g$D <= 1 && g$p_naive >= 0 && g$p_naive <= 1 &&
                g$p_bootstrap >= 0 && g$p_bootstrap <= 1 && g$B_converged > 0)
})

test_that("slt_lrt identifies the two nested pairs, rejects non-nested pairs, and requires matching family", {
  fits <- setup_fits()
  lrt1 <- slt_lrt(fits$shape1, fits$shape11)
  lrt2 <- slt_lrt(fits$common, fits$shape11)
  expect_true(lrt1$df == 1L && lrt2$df == 1L &&
                lrt1$statistic >= 0 && lrt1$p.value >= 0 && lrt1$p.value <= 1)
  expect_error(slt_lrt(fits$shape1, fits$common), "supports")

  d2 <- rslt(200, "lomax", c(0.1, 1), c(0.3, 1), tau1 = fits$tau1, prop_transfer = 1)
  other_fam_fit <- fit_slt(d2, "lomax", "shape11", tau1 = fits$tau1)
  expect_error(slt_lrt(fits$shape1, other_fam_fit), "same family")
})

test_that("slt_model_comparison ranks by AIC by default and names unnamed lists", {
  fits <- setup_fits()
  tab <- slt_model_comparison(list(shape1 = fits$shape1, shape11 = fits$shape11, common = fits$common))
  expect_true(nrow(tab) == 3 && all(diff(tab$AIC) >= -1e-8) && all(tab$converged))

  tab2 <- slt_model_comparison(list(fits$shape1, fits$shape11))
  expect_true(all(grepl("^model", tab2$model)))
})

test_that("slt_model_comparison works with a single model", {
  fits <- setup_fits()
  tab <- slt_model_comparison(list(only_model = fits$shape1))
  expect_equal(nrow(tab), 1)
  expect_true(tab$converged[1])
})
