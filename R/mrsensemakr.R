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
  out <- list()

  # first stage
  fs.form <- make_formula(y = exposure, x = c(instrument, covariates))
  first.stage <- lm(fs.form, data = data)


  # reduced form
  rf.form <- make_formula(y = outcome,  x = c(instrument, covariates))
  reduced.form <- lm(rf.form, data = data)

  info <- list(outcome = outcome,
               exposure = exposure,
               instrument = instrument,
               covariates = covariates)

  # traditional MR
  trad.mr <- mr_estimates(fs = first.stage,
                          rf = reduced.form,
                          exposure = exposure,
                          outcome  = outcome,
                          instrument = instrument,
                          alpha = alpha)
  out$mr <- trad.mr


  # first stage sensitivity
  fs.sense <- sense_trait(model = first.stage,
                          instrument = instrument,
                          trait = exposure,
                          alpha = alpha)

  out$exposure <- list(model = first.stage,
                       sensitivity = fs.sense)

  # reduced form sensitivity
  rf.sense <- sense_trait(model = reduced.form,
                          instrument = instrument,
                          trait = outcome,
                          alpha = alpha)

  out$outcome <- list(model = reduced.form,
                       sensitivity = rf.sense)

  # bounds
  if(!is.null(benchmark_covariates)){

    benchmark_covariates <- lapply(benchmark_covariates,
                                   clean_benchmarks,
                                   data = data,
                                   model = first.stage)

    fs.bounds <- sensemakr::ovb_partial_r2_bound(model = first.stage,
                                                 treatment = instrument,
                                                 benchmark_covariates = benchmark_covariates,
                                                 kd = k)
    names(fs.bounds)[2:3] <- c("r2zw.x", "r2dw.zx")
    out$exposure$bounds <- fs.bounds

    rf.bounds <- sensemakr::ovb_partial_r2_bound(model = reduced.form,
                                                 treatment = instrument,
                                                 benchmark_covariates = benchmark_covariates,
                                                 kd = k)
    names(rf.bounds)[2:3] <- c("r2zw.x", "r2yw.zx")
    out$outcome$bounds <- rf.bounds

  }

  class(out) <- "mr_sensemakr"
  return(out)
}

##'@export
print.mr_sensemakr <- function(x, digits = 2, ...){
  cat("Sensitivity Analysis for Mendelian Randomization (MR)\n", sep ="")
  cat(" Exposure: ", x$mr$exposure, "\n",
      " Outcome: ", x$mr$outcome, "\n",
      " Genetic instrument: ", x$mr$instrument,"\n", sep ="")
  cat("\n")
  print(x$mr)
  cat("\n")
  print(x$exposure$sensitivity, digits = digits)
  cat("\n")
  print(x$outcome$sensitivity, digits = digits)
  cat("\n")
  if(!is.null(x$exposure$bounds) & !is.null(x$outcome$bounds)){
    bounds <- cbind(x$exposure$bounds, x$outcome$bounds[,3,drop = F])
    cat("Bounds on the maximum strength of omitted variables W\n")
    bounds[,2:4] <-  lapply(bounds[,2:4], function(x) paste0(round(x*100, digits = digits), "%"))
    print(bounds, row.names = F)
    cat("\n")
  }
}


clean_benchmarks <- function(bench, data, model){
  classes <- sapply(data[bench], class)
  change  <- which(classes %in% c("factor", "character"))
  for(i in change){
    bench <- c(bench, unique(paste0(bench[i], data[[bench[i]]])))
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
  return(bounds)
}

