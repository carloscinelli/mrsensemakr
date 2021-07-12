# given first stage and reduced form estimates, produces the ILS estimates


ils <- function(fs, rf, instrument, ...){
  UseMethod("ils")
}


#'@importFrom sensemakr model_helper
#'@importFrom stats as.formula coef cor lm pt qt resid
ils.lm <- function(fs, rf, instrument, ...){
  summ.fs <- model_helper(model = fs, covariates = instrument)
  summ.rf <- model_helper(model = rf, covariates = instrument)
  rho <- rho(fs = fs, rf = rf)
  estimates <- ils_estimates(fs.coef = summ.fs$estimate,
                             rf.coef = summ.rf$estimate,
                             fs.se = summ.fs$se,
                             rf.se = summ.rf$se,
                             rho = rho)
  return(estimates)
}


ils_estimates <- function(fs.coef, rf.coef, fs.se, rf.se, rho){

  iv.estimate <- ils_coef(fs.coef = fs.coef, rf.coef = rf.coef)

  iv.se       <- ils_se(fs.coef = fs.coef,
                        rf.coef = rf.coef,
                        fs.se = fs.se,
                        rf.se = rf.se,
                        rho = rho)

  data.frame(estimate = iv.estimate,
             se = iv.se)
}


ils_coef <- function(fs.coef, rf.coef){
  tau <- rf.coef/fs.coef
  return(tau)
}


ils_var <- function(fs.coef, rf.coef, fs.se, rf.se, rho){
  tau <- ils_coef(fs.coef = fs.coef, rf.coef = rf.coef)
  fs.var <- fs.se^2
  rf.var <- rf.se^2
  var.tau <- (1/fs.coef^2)*(rf.var + (tau^2)*fs.var - 2*tau*rho*sqrt(fs.var*rf.var))
  return(var.tau)
}


ils_se <-  function(fs.coef, rf.coef, fs.se, rf.se, rho){
  se.tau <- sqrt(ils_var(fs.coef = fs.coef,
                         rf.coef = rf.coef,
                         fs.se = fs.se,
                         rf.se = rf.se,
                         rho = rho))
  return(se.tau)
}



# Get summary statistics from RF and FS regressions -----------------------


rho <- function(fs, rf) cor(resid(fs), resid(rf))

## These functions get all the data we need for making inference and performing sensitivity analysis

iv_model_helper <- function(...){
  UseMethod("iv_model_helper")
}


iv_model_helper.lm <- function(fs, rf, instrument, ...){

  summ.fs <- sensemakr::model_helper(model = fs, covariates = instrument)
  summ.rf <- sensemakr::model_helper(model = rf, covariates = instrument)
  rho     <- rho(fs, rf)

  if (summ.fs$dof != summ.rf$dof) stop("Degrees of freedom of first-stage and reduced form regressions differ.")

  iv.data <- list(fs.coef   = summ.fs$estimate,
                  fs.se     = summ.fs$se,
                  rf.coef   = summ.rf$estimate,
                  rf.se     = summ.rf$se,
                  rho       = rho,
                  dof       = summ.rf$dof)
  return(iv.data)
}



mr_estimates <- function(fs,
                         rf,
                         exposure = exposure,
                         outcome  = outcome,
                         instrument,
                         alpha){
  trad.mr <- ils(fs = fs, rf = rf, instrument = instrument)
  dof <- rf$df.residual
  t.value <- trad.mr$estimate/trad.mr$se
  crit.value <- qt(p = 1- alpha/2, df = dof)
  ci.low  <- trad.mr$estimate - crit.value*trad.mr$se
  ci.up   <- trad.mr$estimate + crit.value*trad.mr$se
  p.value <- 2*(pt(-abs(t.value), df = dof, lower.tail = T))
  out <- list(
    exposure = exposure,
    outcome  = outcome,
    instrument = instrument,
    estimate = trad.mr$estimate,
    se       = trad.mr$se,
    t.value  = t.value,
    ci.low   = ci.low,
    ci.up    = ci.up,
    p.value  = p.value,
    dof      = dof,
    alpha    = alpha)
  class(out) <- "trad.mr"
  return(out)
}

print.trad.mr <- function(x, digits  = 3,...){
  cat("Traditional MR results (2SLS)\n")

  if(is.null(x)){
    "Both exposure and outcome data are needed for full MR estimates\n"
  } else {
  cat("  MR Estimate (", round((1-x$alpha)*100,digits), "% CI): ",
      round(x$estimate, digits = digits),
      " (",
      round(x$ci.low, digits = digits), " - ",
      round(x$ci.up, digits = digits),
      ")\n",
      "  P-value: ",
      ifelse(x$p.value < 2e-16, " < 2x10^16", x$p.value),
      "\n",
      sep = "")
  }
}
