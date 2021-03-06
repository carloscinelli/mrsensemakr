
##'@export
plot.mr_sensemakr <- function(x,
                              type = c("outcome", "exposure"),
                              benchmark_covariates = NULL,
                              k = 1,
                              alpha = NULL,
                              nlevels = 7,
                              lim.x = NULL,
                              lim.y = NULL,
                              ...
                              ){
  contour.out <- list()

  type <- match.arg(type)

  model <- switch(type,
                  outcome = x$outcome$model,
                  exposure = x$exposure$model)

  t.value <- coef(summary(model))[x$info$instrument, "t value"]

  if(is.null(alpha)){
    alpha <- x$info$alpha
  }
  if(!(alpha >= 0 & alpha <= 1)){
    stop("alpha must be between 0 and 1")
  }
  dof     <- model$df.residual
  t.thr   <- abs(qt(alpha/2, df = dof - 1))*sign(t.value)

  contour.out$info <- list(type = type, k = k, alpha = alpha, t.thr = t.thr)


  rv <- x[[type]]$sensitivity$rv

    if(!is.null(benchmark_covariates)){

      benchmark_covariates <- lapply(benchmark_covariates,
                                     clean_benchmarks,
                                     model = model)

      bounds <- multiple_bounds(model = model,
                                covariate = x$mr$instrument,
                                benchmark_covariates = benchmark_covariates,
                                k = k,
                                alpha = alpha)

      contour.out$bounds <- bounds

      max.r2dz.x  <- max(c(bounds$r2dz.x*1.5,  1.5*rv))
      max.r2yz.dx <- max(c(bounds$r2yz.dx*1.5, 1.5*rv))

      if(is.null(lim.x)){
        lim.x <- max.r2dz.x
      }

      if(is.null(lim.y)){
        lim.y <- max.r2yz.dx
      }
    }

  if(is.null(lim.x)){
    lim.x <- max(c(rv*1.5, 0.001))
  }

  if(is.null(lim.y)){
    lim.y <- max(c(rv*1.5, 0.001))
  }

  xlab <- expression(paste("Partial ", R^2, " of unobservables with gen. instrument"))

  if(type == "outcome"){
    ylab <- expression(paste("Partial ", R^2, " of unobservables with outcome trait"))
  }else{
    ylab <- expression(paste("Partial ", R^2, " of unobservables with exposure trait"))
  }

  contour.out$graphics <- list(lim.x = lim.x, lim.y = lim.y, nlevels = nlevels, xlab = xlab, ylab = ylab)

  sensemakr::ovb_contour_plot(model = model,
                              treatment = x$mr$instrument,
                              sensitivity.of = "t-value",
                              t.threshold = t.thr,
                              xlab = xlab,
                              ylab = ylab,
                              lim = lim.x,
                              lim.y = lim.y,
                              nlevels = nlevels) -> contours

  names(contours) <- c("r2zw.x", "r2yw.zx", "value")
  contour.out$contours <- contours

  if(!is.null(benchmark_covariates)){
    sensemakr::add_bound_to_contour(bounds, treatment = x$mr$instrument)
  }

  invisible(contour.out)

}
