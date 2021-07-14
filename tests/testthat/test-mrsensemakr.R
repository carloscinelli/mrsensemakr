test_that("alpha", {
  # cleans workspace
  rm(list = ls())

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

test_that("sign and alpha for plots",
          {
            # cleans workspace
            rm(list = ls())

            # loads pacakges
            library(mrsensemakr)
            library(testthat)

            ## simulated data example
            data("sim_data")

            ## changes sign
            sim_data$prs <- sim_data$prs*(-1)

            ## create vectors indicating variable names in the data
            outcome    <- "out.trait" # name of outcome trait
            exposure   <- "exp.trait" # name of exposure trait
            instrument <- "prs" # genetic instrument (e.g, polygenic risk score)
            age.sex    <- c("age", "sex") # age and sex variables (if applicable)
            alc.smok   <- c("alcohol", "smoking") # putative pleoitropic vars.
            pcs        <- paste0("pc", 1:20) # first 20 principal components pc1 ... pc20

            ## runs MR sensitivity analysis
            mr.sense.neg <- mr_sensemakr(outcome = outcome,
                                        exposure = exposure,
                                        instrument = instrument,
                                        covariates = c(age.sex, alc.smok, pcs),
                                        data = sim_data,
                                        alpha = 0.05,
                                        benchmark_covariates = list(alc.smok = alc.smok,
                                                                    pcs = pcs))
            print <- capture_output(print(mr.sense.neg))

            check_plot <- plot(mr.sense.neg,
                               benchmark_covariates = list(alc.smok = alc.smok, pcs = pcs),
                               k = list(alc.smok = 1, pcs = 1))
            expect_equal(check_plot$info$t.thr, -1.96240507)

            check_plot2 <- plot(mr.sense.neg,
                               benchmark_covariates = list(alc.smok = alc.smok, pcs = pcs),
                               k = 1)
            expect_equal(check_plot$contours, check_plot2$contours)

            check_plot3 <- plot(mr.sense.neg, alpha = 0.1,
                                benchmark_covariates = list(alc.smok = alc.smok, pcs = pcs),
                                k = 1)
            expect_equal(check_plot3$info$t.thr, -1.6464212)

            check_plot4 <- plot(mr.sense.neg)
            expect_equal(check_plot4$graphics$lim.x, mr.sense.neg$outcome$sensitivity$rv*1.5, ignore_attr = T)

            check_plot5 <- plot(mr.sense.neg, type = "exposure")
            expect_equal(check_plot5$graphics$lim.x, mr.sense.neg$exposure$sensitivity$rv*1.5, ignore_attr = T)

            expect_error(plot(mr.sense.neg, alpha = -0.1))

          }
          )

test_that("NA",{
  # cleans workspace
  rm(list = ls())

  # loads pacakges
  library(mrsensemakr)
  library(testthat)

  ## simulated data example
  data("sim_data")

  ## creates NAs
  na.prs <- sample(nrow(sim_data), 10)
  na.out <- sample(nrow(sim_data), 15)
  na.exp <- sample(nrow(sim_data), 20)
  sim_data$prs[na.prs] <- NA
  sim_data$exp.trait[na.exp] <- NA
  sim_data$out.trait[na.out] <- NA

  ## create vectors indicating variable names in the data
  outcome    <- "out.trait" # name of outcome trait
  exposure   <- "exp.trait" # name of exposure trait
  instrument <- "prs" # genetic instrument (e.g, polygenic risk score)
  age.sex    <- c("age", "sex") # age and sex variables (if applicable)
  alc.smok   <- c("alcohol", "smoking") # putative pleoitropic vars.
  pcs        <- paste0("pc", 1:20) # first 20 principal components pc1 ... pc20

  expect_warning(mr.sense.na <- mr_sensemakr(outcome = outcome,
                               exposure = exposure,
                               instrument = instrument,
                               covariates = c(age.sex, alc.smok, pcs),
                               data = sim_data,
                               alpha = 0.05,
                               benchmark_covariates = list(alc.smok = alc.smok,
                                                           pcs = pcs)))

})

test_that("errors", {
  rm(list = ls())

  # loads pacakges
  library(mrsensemakr)
  library(testthat)

  ## simulated data example
  data("sim_data")

  ## changes sign
  sim_data$prs <- sim_data$prs*(-1)

  ## create vectors indicating variable names in the data
  outcome    <- "out.trait" # name of outcome trait
  exposure   <- "exp.trait" # name of exposure trait
  instrument <- "prs" # genetic instrument (e.g, polygenic risk score)
  age.sex    <- c("age", "sex") # age and sex variables (if applicable)
  alc.smok   <- c("alcohol", "smoking") # putative pleoitropic vars.
  pcs        <- paste0("pc", 1:20) # first 20 principal components pc1 ... pc20

  ## runs MR sensitivity analysis
  expect_error(mr.sense.neg <- mr_sensemakr(outcome = "a",
                               exposure = exposure,
                               instrument = instrument,
                               covariates = c(age.sex, alc.smok, pcs),
                               data = sim_data,
                               alpha = 0.05,
                               benchmark_covariates = list(alc.smok = alc.smok,
                                                           pcs = pcs)))

  expect_error(mr.sense.neg <- mr_sensemakr(outcome = outcome,
                                            exposure = "b",
                                            instrument = instrument,
                                            covariates = c(age.sex, alc.smok, pcs),
                                            data = sim_data,
                                            alpha = 0.05,
                                            benchmark_covariates = list(alc.smok = alc.smok,
                                                                        pcs = pcs)))

  expect_error(mr.sense.neg <- mr_sensemakr(outcome = outcome,
                                            exposure = exposure,
                                            instrument = "c",
                                            covariates = c(age.sex, alc.smok, pcs),
                                            data = sim_data,
                                            alpha = 0.05,
                                            benchmark_covariates = list(alc.smok = alc.smok,
                                                                        pcs = pcs)))

  expect_error(mr.sense.neg <- mr_sensemakr(outcome = outcome,
                                            exposure = exposure,
                                            instrument = instrument,
                                            covariates = "d",
                                            data = sim_data,
                                            alpha = 0.05,
                                            benchmark_covariates = list(alc.smok = alc.smok,
                                                                        pcs = pcs)))

  expect_error(mr.sense.neg <- mr_sensemakr(outcome = outcome,
                                            exposure = exposure,
                                            instrument = instrument,
                                            covariates =  c(age.sex, alc.smok, pcs),
                                            data = mtcars,
                                            alpha = 0.05,
                                            benchmark_covariates = list(alc.smok = alc.smok,
                                                                        pcs = pcs)))


  expect_error(mr.sense.neg <- mr_sensemakr(outcome = outcome,
                                            exposure = exposure,
                                            instrument = instrument,
                                            covariates =  c(age.sex, alc.smok, pcs),
                                            data = sim_data,
                                            alpha = 0.05,
                                            benchmark_covariates = "e"))

  expect_error(mr.sense.neg <- mr_sensemakr(outcome = outcome,
                                            exposure = exposure,
                                            instrument = instrument,
                                            data = sim_data,
                                            alpha = 0.05,
                                            benchmark_covariates = age.sex))

  expect_error(mr.sense.neg <- mr_sensemakr(outcome = outcome,
                                            exposure = exposure,
                                            instrument = instrument,
                                            covariates = alc.smok,
                                            data = sim_data,
                                            alpha = 0.05,
                                            benchmark_covariates = age.sex))

  expect_error(mr.sense.neg <- mr_sensemakr(outcome = outcome,
                                            exposure = NULL,
                                            instrument = instrument,
                                            data = sim_data,
                                            alpha = 0.05))

  expect_error(plot(mr.sense.neg, alpha = -0.1))
  expect_error(plot(mr.sense.neg, alpha = 1.1))

  file.remove("Rplots.pdf")

})
