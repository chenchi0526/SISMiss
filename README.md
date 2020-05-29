SIsMiss
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

`SIsMiss` is an R package for variable **s**election and statistical
**i**nference under **s**hadow variable assumption for a linear
regression with **miss**ing subjects in response. `SIsMiss` is
applicable to linear regression with regularization and without
regularization. Robust inference will be obtained for both scenario of
missing at random (MAR) and missing not at random (MNAR). The estimates
and standard error for the unknown regression coefficients will be
returned, along with optional confidence intervals.

The underlying estimating method is based on conditional likelihood
discussed in paper (…).

## Installation

You can install SIsMiss from github with:

``` r
# install.packages("devtools")
devtools::install_github("chenchi0526/SIsMiss")
```

## Example

For simplicity, consider a linear regression with no missing subjects.

``` r
rm(list = ls())
library(SIsMiss)
n <- 50
p <- 8
beta <- c(3, 0, 1.5, 0, 2, rep(0, p-5))
gamma <- 3
xm <- matrix(rnorm(n*p), ncol = p, nrow = n)
z <- rnorm(n, 0, 1)
y <- xm %*% beta + gamma*z + rnorm(n)
```

### Unregularized linear regression

For unregularized linear regression, the standard error can be estimated
via asymptotic theory or perturbation method.

When estimating standard error via asymptotic theory, the symmetric
asymptotic confidence interval will be returned.

``` r
SIsMiss(y, z, u, regularize = FALSE, cov.names = NULL,
        se.method = "asymp", CI.alpha = 0.05,
        M = NULL, seed_num = NULL)
```

When estimating standard error via perturbation, the lower bound and
upper bound for confidence interval are the \(\alpha/2\)-th quantile and
\(1-\alpha/2\)-th quantile for the samples of perturbated estimates.

``` r
SIsMiss(y, z, u, regularize = FALSE, cov.names = NULL,
        se.method = "perturb", CI.alpha = 0.05,
        M = 200, seed_num = 123)
```

### Regularized linear regression

For regularized linear regression, the adaptive LASSO penalty is
considered where the tuning parameter is determined by BIC. The standard
error of coefficients is estimated via perturbation method only.

``` r
SIsMiss(y, z, u, regularize = TRUE, cov.names = NULL,
        se.method = "perturb", CI.alpha = 0.05,
        M = 200, seed_num = 123)
```
