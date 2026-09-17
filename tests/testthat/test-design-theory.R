# Checks the package against Table 5 of Laumen & Cramer (2019, Lemma 4.1),
# matched exactly, plus rslt()'s design-generation behaviour.


test_that("slt_prob_stage1() exactly reproduces Table 5 of Laumen & Cramer (2019)", {
  r1 <- slt_prob_stage1(12, 0.5, c(1, 1), design = "typeP", pi1 = 0.25)
  expect_equal(r1$n_star, 8L)
  expect_equal(r1$prob, 0.9864756680, tolerance = 1e-9)

  r2 <- slt_prob_stage1(12, 0.5, c(1, 1), design = "typeP", pi1 = 0.50)
  expect_equal(r2$n_star, 10L)
  expect_equal(r2$prob, 0.9997315155, tolerance = 1e-9)

  r3 <- slt_prob_stage1(12, 0.5, c(1, 1), design = "typeP", pi1 = 0.75)
  expect_equal(r3$n_star, 10L)
  expect_equal(r3$prob, 0.9997315155, tolerance = 1e-9)

  for (R0 in c(3, 6, 9)) {
    r4 <- slt_prob_stage1(12, 0.5, c(1, 1), design = "typeM", R0_1 = R0)
    expect_equal(r4$n_star, 11L)
    expect_equal(r4$prob, 0.9999862301, tolerance = 1e-9)
  }
})

test_that("slt_prob_stage1() validates its inputs and handles boundary cases", {
  expect_error(slt_prob_stage1(12, 0.5, c(1, 1), design = "typeP"), "pi1")
  expect_error(slt_prob_stage1(12, 0.5, c(1, 1), design = "typeM"), "R0_1")
  expect_equal(slt_prob_stage1(12, 0.5, c(1, 1), design = "typeM", R0_1 = 0)$prob, 0)
  # pi1 = 1 (transfer every survivor) must match Type-M's n* = n - 1 exactly,
  # since both then require only "at least one survivor" (D1 <= n - 1)
  r_full <- slt_prob_stage1(12, 0.5, c(1, 1), design = "typeP", pi1 = 1)
  r_typeM <- slt_prob_stage1(12, 0.5, c(1, 1), design = "typeM", R0_1 = 12)
  expect_equal(r_full$n_star, 11L)
  expect_equal(r_full$prob, r_typeM$prob)
})

test_that("rslt() Type-M respects R0_1, capped by available survivors", {
  set.seed(5)
  theta0 <- c(0.08, 1.4); theta1 <- c(0.25, 1.4)
  tau1 <- slt_tau("weibull", theta0, 0.5)
  d <- rslt(500, "weibull", theta0, theta1, tau1 = tau1, design = "typeM", R0_1 = 100)
  survivors <- 500 - d$n_early
  expect_equal(d$n_transferred, min(survivors, 100))

  # R0_1 larger than any plausible survivor count: transfers ALL survivors
  d2 <- rslt(50, "weibull", theta0, theta1, tau1 = tau1, design = "typeM", R0_1 = 10000)
  expect_equal(d2$n_transferred, 50 - d2$n_early)
})

test_that("rslt() validates design-specific inputs", {
  theta0 <- c(0.08, 1.4); theta1 <- c(0.25, 1.4)
  expect_error(rslt(10, "weibull", theta0, theta1, tau1 = 1, design = "typeM"), "R0_1")
  expect_error(rslt(10, "weibull", theta0, theta1, tau1 = 1, design = "typeP", prop_transfer = 1.5),
               "prop_transfer")
})

test_that("rslt()'s default design = 'typeP' keeps the old call signature working", {
  set.seed(6)
  theta0 <- c(0.08, 1.4); theta1 <- c(0.25, 1.4)
  tau1 <- slt_tau("weibull", theta0, 0.5)
  d <- rslt(200, "weibull", theta0, theta1, tau1 = tau1, prop_transfer = 0.5)
  expect_true(is.list(d) && all(c("Y", "Z", "n_early", "n_transferred") %in% names(d)))
})
