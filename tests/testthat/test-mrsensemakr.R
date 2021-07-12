test_that("handles NA", {
  data("sim_data")
  n <- nrow(sim_data)
  na.out.trait <- sample(n, 20)
  na.exp.trait <- sample(n, 20)
  sim_data$out.trait[na.out.trait] <- NA
  sim_data$exp.trait[na.exp.trait] <- NA
  expect_warning(mr.out <- mr_sensemakr(outcome = "out.trait", exposure = "exp.trait", instrument = "prs", data = sim_data))
})



test_that("negative value plot", {
  data("sim_data")
  ## create vectors indicating variable names in the data
  outcome    <- "out.trait" # name of outcome trait
  exposure   <- "exp.trait" # name of exposure trait
  instrument <- "prs" # genetic instrument (e.g, polygenic risk score)
  age.sex    <- c("age", "sex") # age and sex variables (if applicable)
  alc.smok   <- c("alcohol", "smoking") # putative pleoitropic vars.
  pcs        <- paste0("pc", 1:20) # first 20 principal components pc1 ... pc20

  sim_data$prs2 <- - sim_data$prs
  ## runs MR sensitivity analysis
  mr.sense <- mr_sensemakr(outcome = outcome,
                           exposure = exposure,
                           instrument = "prs2",
                           covariates = c(age.sex, alc.smok, pcs),
                           data = sim_data,
                           benchmark_covariates = list(alc.smok = alc.smok,
                                                       pcs = pcs), alpha = 0.1)
  ## print results
  mr.sense
  ## sensitivity contour plots
  plot(mr.sense,
       benchmark_covariates = list(alc.smok = alc.smok, pcs = pcs),
       k = list(alc.smok = 25, pcs = 35))

})
