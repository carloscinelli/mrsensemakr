make_formula <- function(y,x){
  x   <- paste(x, collapse = " + ")
  yx  <- paste(y, x, sep = " ~ ")
  as.formula(yx)
}



sense_trait <- function(model,
                        instrument,
                        trait,
                        alpha){
  r2 <- sensemakr::partial_r2(model, covariates = instrument)
  rv <- sensemakr::robustness_value(model, covariates = instrument, q = 1, alpha = alpha)
  out <- list(instrument = instrument,
              trait      = trait,
              type.trait = deparse(substitute(trait)),
              partial.r2 = r2,
              rv = rv)
  class(out) <- "sense_trait"
  return(out)
}


##'@export
print.sense_trait <- function(x, digits = 4, ...){
  cat("Sensitivity genetic instrument (",
      x$instrument,
      ") -> ",
      x$type.trait,
      " (",
      x$trait, ")\n",
      "  Partial R2: ", round(x$partial.r2*100, digits), "%\n",
      "  RV (alpha = ", attr(x$rv, "alpha"),
      "): ",
      round(x$rv*100, digits), "%\n",
      sep ="")
}

##' MR-sensemakr
##'
##'
##'@param outcome A character vector with the name of the outcome trait.
##'@param exposure A character vector with the name of the exposure trait.
##'@param instrument A character vector with the name of the genetic instrument.
##'@param covariates A character vector with the name of the control covariates, such as age, sex, genomic principal components, batch effect dummies and putative pleiotropic pathways.
##'@param data An object of the class data.frame containing the variables used in the analysis.
##'@param benchmark_covariates Covariates for benchmarking. Must be a subset of the \code{covariates} argument. The user has two options: (i) character vector of the names of covariates that will be used to bound the plausible strength of the unobserved confounders. Each variable will be considered separately; (ii) a named list with character vector names of covariates that will be used, as a group, to bound the plausible strength of the unobserved confounders. The names of the list will be used for the benchmark labels.
##'
##'@param k numeric vector. Parameterizes how many times stronger residual biases are related to the treatment and the outcome in comparison to the observed benchmark covariates.
##'@param alpha significance level
##'
##'@examples
## loads package
##'library(mrsensemakr)
##'
##'## simulated data example
##'data("sim_data")
##'
##'## create vectors indicating variable names in the data
##'outcome    <- "out.trait" # name of outcome trait
##'exposure   <- "exp.trait" # name of exposure trait
##'instrument <- "prs" # genetic instrument (e.g, polygenic risk score)
##'age.sex    <- c("age", "sex") # age and sex variables (if applicable)
##'alc.smok   <- c("alcohol", "smoking") # putative pleoitropic vars.
##'pcs        <- paste0("pc", 1:20) # first 20 principal components pc1 ... pc20
##'
##'## runs MR sensitivity analysis
##'mr.sense <- mr_sensemakr(outcome = outcome,
##'                         exposure = exposure,
##'                         instrument = instrument,
##'                         covariates = c(age.sex, alc.smok, pcs),
##'                         data = sim_data,
##'                         benchmark_covariates = list(alc.smok = alc.smok,
##'                                                     pcs = pcs))
##'## print results
##'mr.sense
##'## sensitivity contour plots
##'plot(mr.sense,
##'     benchmark_covariates = list(alc.smok = alc.smok, pcs = pcs),
##'     k = list(alc.smok = 25, pcs = 35))
##'@importFrom stats na.omit
##'@export
mr_sensemakr <- function(outcome,
                         exposure,
                         instrument,
                         covariates = NULL,
                         data,
                         benchmark_covariates = NULL,
                         k = 1,
                         alpha = 0.05){

  # coerce to data.frame
  data <- as.data.frame(data)

  # remove NAs
  any.na <- any(is.na(data))
  if(any.na){
    data <- na.omit(data)
    warning("Missing values (NA) were found and were omitted from the analysis.")
    NAs <- "Missing data removed from analysis."
  } else {
    NAs <- "No missing data found."
  }

  # check names on data.frame
  all.names <- unique(unlist(c(outcome, exposure,instrument, covariates, benchmark_covariates)))
  data.names <- names(data)
  not.ok <- which(!all.names %in% data.names)
  if (length(not.ok) > 0){
    stop(paste("variables", paste(all.names[not.ok], collapse = ", "), "were not found"))
  }

  # check benchmark and covariates
  benchmark_covariates_check <- unique(unlist(benchmark_covariates))
  if (is.null(covariates) & !is.null(benchmark_covariates_check)){
    stop("'benchmark_covariates' must be a subset of 'coviates:' argument `covariates` not provided.")
  }
  bench.not.ok <- which(!benchmark_covariates_check %in% covariates)
  if( length(bench.not.ok) > 0){
    stop("'benchmark_covariates' must be a subset of 'coviates. Variables not found: ",
         paste(benchmark_covariates_check[bench.not.ok], collapse = ", "))
  }

  # output list
  out <- list()

  # info
  out$info     <- list(outcome    = outcome,
                       exposure   = exposure,
                       instrument = instrument,
                       covariates = covariates,
                       alpha= alpha,
                       NAs = NAs)

  # check if either exposure or outcome was provided
  if (is.null(outcome) | is.null(exposure)){
    stop("Both the outcome trait and exposure trait must be provided")
  }

  # first stage
  fs.form      <- make_formula(y = exposure, x = c(instrument, covariates))
  first.stage  <- lm(fs.form, data = data)

  # first stage sensitivity
  fs.sense     <- sense_trait(model      = first.stage,
                              instrument = instrument,
                              trait      = exposure,
                              alpha      = alpha)

  out$exposure <- list(model       = first.stage,
                       sensitivity = fs.sense)

  if(!is.null(benchmark_covariates)){

    benchmark_covariates <- lapply(benchmark_covariates,
                                   clean_benchmarks,
                                   model = first.stage)

    fs.bounds <- sensemakr::ovb_bounds(model = first.stage,
                                       treatment = instrument,
                                       benchmark_covariates = benchmark_covariates,
                                       kd = k,
                                       alpha = alpha)
    fs.bounds <- fs.bounds[c(1,2,3,7)]
    names(fs.bounds)[2:3] <- c("r2zw.x", "r2dw.zx")
    names(fs.bounds)[4] <- c("adjusted_t_exposure")
    out$exposure$bounds <- fs.bounds
  }


  # reduced form

  rf.form      <- make_formula(y = outcome,  x = c(instrument, covariates))
  reduced.form <- lm(rf.form, data = data)


  # reduced form sensitivity
  rf.sense     <- sense_trait(model      = reduced.form,
                              instrument = instrument,
                              trait      = outcome,
                              alpha      = alpha)

  out$outcome  <- list(model       = reduced.form,
                       sensitivity = rf.sense)

  if(!is.null(benchmark_covariates)){

    benchmark_covariates <- lapply(benchmark_covariates,
                                   clean_benchmarks,
                                   model = reduced.form)

    rf.bounds <- sensemakr::ovb_bounds(model = reduced.form,
                                       treatment = instrument,
                                       benchmark_covariates = benchmark_covariates,
                                       kd = k, alpha = alpha)
    rf.bounds <- rf.bounds[c(1,2,3,7)]
    names(rf.bounds)[2:3] <- c("r2zw.x", "r2yw.zx")
    names(rf.bounds)[4] <- c("adjusted_t_outcome")
    out$outcome$bounds <- rf.bounds
  }



  # traditional MR

  trad.mr <- mr_estimates(fs = first.stage,
                          rf = reduced.form,
                          exposure = exposure,
                          outcome  = outcome,
                          instrument = instrument,
                          alpha = alpha)
  out$mr <- trad.mr

  class(out) <- "mr_sensemakr"
  return(out)
}

##'@export
print.mr_sensemakr <- function(x, digits = 2, ...){
  cat("Sensitivity Analysis for Mendelian Randomization (MR)\n", sep ="")
  cat(" Exposure: ", x$info$exposure, "\n",
      " Outcome: ", x$info$outcome, "\n",
      " Genetic Instrument: ", x$info$instrument,"\n",
      " Missing Data: ", x$info$NAs,"\n",
      sep ="")
  cat("\n")
  print(x$mr)
  cat("\n")
  print(x$exposure$sensitivity, digits = digits)
  cat("\n")
  print(x$outcome$sensitivity, digits = digits)
  cat("\n")
  if(!is.null(x$exposure$bounds) & !is.null(x$outcome$bounds)){
    bounds <- cbind(x$exposure$bounds[,1:3, drop = F],
                    x$outcome$bounds[,3,drop = F],
                    x$exposure$bounds[,4,drop = F],
                    x$outcome$bounds[,4,drop = F])
    cat("Bounds on the maximum explanatory power of omitted variables W, if it were as strong as:\n")
    bounds[,2:4] <-  lapply(bounds[,2:4], function(x) paste0(round(x*100, digits = digits), "%"))
    print(bounds, row.names = F)
    cat("\n")
  }
}


clean_benchmarks <- function(bench, model){
  classes <- sapply(model$model[bench], class)
  change  <- which(classes %in% c("factor", "character"))
  for(i in change){
    bench <- c(bench, unique(paste0(bench[i], model$model[[bench[i]]])))
  }
  intersect(bench,names(coef(model)))
}


multiple_bounds <- function(model,
                            covariate,
                            benchmark_covariates,
                            k,
                            alpha = 0.05){

  check_is_list   <- is.list(benchmark_covariates)
  bench_names     <- names(benchmark_covariates)
  check_null_name <- is.null(bench_names)

  if(!check_is_list | check_null_name){
    stop("benchmark_covariates should be a named list")
  }

  # if k is numeric vector, replicate it for each benchmark list

  if(is.numeric(k)){
    num <- k
    k   <- list()
    for(i in bench_names){
      k[[i]] <- num
    }
  }

  k_names           <- names(k)
  check_names_match <- all(bench_names == k_names)

  if(!check_names_match){
    stop("names of 'k' should match names of 'benchmark_covariates'")
  }


  bounds <- list()
  for(i in bench_names){
    bounds[[i]] <- sensemakr::ovb_bounds(model = model,
                                         treatment = covariate,
                                         benchmark_covariates = benchmark_covariates[i],
                                         kd = k[[i]],
                                         alpha = alpha)
  }
  bounds <- do.call("rbind", bounds)
  row.names(bounds) <- NULL
  return(bounds)
}

