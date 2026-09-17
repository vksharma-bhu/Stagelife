test_that("rslt(prop_transfer = 1) returns sorted, positive, correctly partitioned failure times", {
  set.seed(1)
  theta0 <- c(0.05, 1.4); theta1 <- c(0.2, 1.4); tau1 <- 5
  d <- rslt(300, "weibull", theta0, theta1, tau1, prop_transfer = 1)
  expect_true(all(d$Y > 0) && all(d$Y <= tau1) && all(d$Z > tau1) &&
                !is.unsorted(d$Y) && !is.unsorted(d$Z) &&
                length(d$Y) + length(d$Z) == 300)
})

test_that("rslt validates its inputs", {
  expect_error(rslt(10, "weibull", c(-1, 1), c(1, 1), tau1 = 5, prop_transfer = 1), "positive")
  expect_error(rslt(10, "weibull", c(1, 1), c(1, 1), tau1 = 0, prop_transfer = 1), "tau1")
})

test_that("rslt respects prop_transfer, including the prop_transfer = 0 edge case", {
  set.seed(2)
  theta0 <- c(0.05, 1.4); theta1 <- c(0.2, 1.4)
  tau1 <- slt_tau("weibull", theta0, 0.5)
  d <- rslt(1000, "weibull", theta0, theta1, tau1 = tau1, prop_transfer = 0.4)
  survivors <- 1000 - d$n_early
  expect_true(length(d$Y) + length(d$Z) == 1000 &&
                d$n_transferred == length(d$Z) &&
                d$n_transferred == floor(0.4 * survivors) &&
                all(d$Z > tau1))

  d0 <- rslt(200, "weibull", theta0, theta1, tau1 = tau1, prop_transfer = 0)
  expect_true(length(d0$Z) == 0L && length(d0$Y) == 200L)
})

test_that("rslt validates prop_transfer range", {
  theta0 <- c(0.05, 1.4); theta1 <- c(0.2, 1.4)
  expect_error(rslt(10, "weibull", theta0, theta1, tau1 = 1, prop_transfer = 1.5), "prop_transfer")
})

test_that("simulated data are consistent with the generating model's cdf (KS sanity check)", {
  set.seed(4)
  theta0 <- c(0.05, 1.4); theta1 <- c(0.2, 1.4); tau1 <- 5
  d <- rslt(5000, "weibull", theta0, theta1, tau1, prop_transfer = 1)
  Tobs <- c(d$Y, d$Z)
  ks <- suppressWarnings(stats::ks.test(Tobs, function(t) pslt(t, "weibull", theta0, theta1, tau1)))
  expect_gt(ks$p.value, 0.01)
})
