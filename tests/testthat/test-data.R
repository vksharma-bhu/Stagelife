test_that("all four bundled datasets have the expected structure", {
  ok <- vapply(c("slt_han_kundu_2014", "slt_bobotas_kateri_2015", "slt_wang_fei_2003", "slt_zhu_2010"),
               function(nm) {
    dat <- get(nm, envir = asNamespace("stagelife"))
    all(c("y0", "y1", "tau1") %in% names(dat)) &&
      is.numeric(dat$y0) && length(dat$y0) > 0 &&
      is.numeric(dat$y1) && length(dat$y1) > 0 &&
      all(dat$y0 <= dat$tau1) && all(dat$y1 > dat$tau1)
  }, logical(1))
  expect_true(all(ok), info = paste("Failed for:", paste(names(ok)[!ok], collapse = ", ")))
})

test_that("a real dataset can be fit end-to-end", {
  fit <- fit_slt(list(Y = slt_han_kundu_2014$y0, Z = slt_han_kundu_2014$y1),
                    "weibull", "shape1", tau1 = slt_han_kundu_2014$tau1)
  expect_true(fit$converged && is.finite(fit$loglik))
})
