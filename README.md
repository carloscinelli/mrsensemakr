
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MR-SENSEMAKR

<!-- badges: start -->

[![Codecov test
coverage](https://codecov.io/gh/carloscinelli/mrsensemakr/branch/master/graph/badge.svg)](https://codecov.io/gh/carloscinelli/mrsensemakr?branch=master)
[![Travis build
status](https://travis-ci.com/carloscinelli/mrsensemakr.svg?branch=master)](https://travis-ci.com/carloscinelli/mrsensemakr)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/carloscinelli/mrsensemakr?branch=master&svg=true)](https://ci.appveyor.com/project/carloscinelli/mrsensemakr)
<!-- badges: end -->

The R package `mrsensemakr` implements sensitivity analysis tools for
Mendelian Randomization, as discussed in [Cinelli et al (2020). Robust
Mendelian randomization in the presence of residual population
stratification, batch effects and horizontal
pleiotropy.](https://www.biorxiv.org/content/10.1101/2020.10.21.347773v1)

## Development version

<!-- commit -->

To install the development version on GitHub make sure you have the R
package `devtools` installed. Also make sure to have the latest version
of `sensemakr` [(link)](https://github.com/carloscinelli/sensemakr)
installed.

``` r
# install.packages("devtools")
devtools::install_github("carloscinelli/sensemakr")
devtools::install_github("carloscinelli/mrsensemakr")
```

CRAN version coming soon.

## Basic usage

``` r
## loads package
library(mrsensemakr)

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
mr.sense <- mr_sensemakr(outcome = outcome,
                         exposure = exposure,
                         instrument = instrument,
                         covariates = c(age.sex, alc.smok, pcs), 
                         data = sim_data, 
                         benchmark_covariates = list(alc.smok = alc.smok,
                                                     pcs = pcs))
## print results
mr.sense
#> Sensitivity Analysis for Mendelian Randomization (MR)
#>  Exposure: exp.trait
#>  Outcome: out.trait
#>  Genetic Instrument: prs
#>  Missing Data: No missing data found.
#> 
#> Traditional MR results (2SLS)
#>   MR Estimate (95% CI): 0.807 (0.111 - 1.503)
#>   P-value: 0.02317659
#> 
#> Sensitivity genetic instrument (prs) -> exposure (exp.trait)
#>   Partial R2: 1%
#>   RV (alpha = 0.05): 3.71%
#> 
#> Sensitivity genetic instrument (prs) -> outcome (out.trait)
#>   Partial R2: 0.65%
#>   RV (alpha = 0.05): 1.76%
#> 
#> Bounds on the maximum explanatory power of omitted variables W, if it were as strong as:
#>  bound_label r2zw.x r2dw.zx r2yw.zx adjusted_t_exposure adjusted_t_outcome
#>  1x alc.smok  0.41%   0.05%   0.76%            3.091924            2.34717
#>       1x pcs  2.27%   3.14%   2.51%            2.309587            1.76642


## sensitivity contour plots
plot(mr.sense, 
     benchmark_covariates = list(alc.smok = alc.smok,
                                 pcs = pcs),
     k = list(alc.smok = 1, pcs = 1))
```

<img src="man/figures/README-example-1.png" width="100%" />

## Simulations

Code to reproduce the simulations of the paper can be found in the
vignettes folder.
