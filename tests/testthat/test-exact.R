test_that("exact simulated example", {
  rm(list = ls())
  n <- 1e3
  x <- sensemakr:::resid_maker(n, rep(1, n))
  w <- sensemakr:::resid_maker(n, x)
  u <- sensemakr:::resid_maker(n, cbind(x, w))
  z <- x + w + sensemakr:::resid_maker(n, cbind(x, w, u))
  d <- z + x+ w + u + sensemakr:::resid_maker(n, cbind(x,w,u,z))
  y <- 0*d + x+ w + u + sensemakr:::resid_maker(n, cbind(x,w,u,z))
  data <- data.frame(y,d,z,x)
  mr.sense <- mr_sensemakr(outcome = "y",
                           exposure = "d",
                           instrument = "z",
                           covariates = "x",
                           benchmark_covariates = list(x = "x"),
                           alpha = 0.05,
                           data = data)
  expect_equivalent(mr.sense$outcome$bounds$adjusted_t_outcome, 0)
})
