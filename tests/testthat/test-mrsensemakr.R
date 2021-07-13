test_that("multiplication works", {
  # loads package
  library(mrsensemakr)
  library(testthat)

  ## simulated data example
  data("sim_data")

  ## create vectors indicating variable names in the data
  outcome    <- "out.trait" # name of outcome trait
  exposure   <- "exp.trait" # name of exposure trait
  instrument <- "prs" # genetic instrument (e.g, polygenic risk score)
  age.sex    <- c("age", "sex") # age and sex variables (if applicable)
  alc.smok   <- c("alcohol", "smoking") # putative pleoitropic vars.
  pcs        <- paste0("pc", 1:20) # first 20 principal components pc1 ... pc20

  ## runs MR sensitivity analysis
  mr.sense.05 <- mr_sensemakr(outcome = outcome,
                           exposure = exposure,
                           instrument = instrument,
                           covariates = c(age.sex, alc.smok, pcs),
                           data = sim_data,
                           alpha = 0.05,
                           benchmark_covariates = list(alc.smok = alc.smok,
                                                       pcs = pcs))
  ## print results
  mr.sense.01 <- mr_sensemakr(outcome = outcome,
                              exposure = exposure,
                              instrument = instrument,
                              covariates = c(age.sex, alc.smok, pcs),
                              data = sim_data,
                              alpha = 0.01,
                              benchmark_covariates = list(alc.smok = alc.smok,
                                                          pcs = pcs))
  expect_equal(mr.sense.01$outcome$sensitivity$rv, 0, ignore_attr=T)
  expect_equal(mr.sense.05$outcome$sensitivity$rv, 0.0176069377, ignore_attr=T)
})
