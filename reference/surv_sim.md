# Survival data generation using either weibull or inverse CDF method

Survival data generation using either weibull or inverse CDF method

## Usage

``` r
surv_sim(
  n_patients,
  n_vars,
  vars_cate,
  scale_weibull,
  shape_weibull,
  beta,
  rateC,
  seed = 926,
  method = c("InverseCDF", "Weibull")
)
```

## Arguments

- n_patients:

  number of patients

- n_vars:

  number of covariates

- vars_cate:

  vector of "continuous", "binary"

- scale_weibull:

  scale parameter for baseline hazard of weibull distribution

- shape_weibull:

  shape parameter of weibull, shape_weibull \< 1 decreasing hazard,
  shape_weibull \> 1 increasing hazard, shape_weibull=1 constant hazard

- beta:

  Fixed effect parameter, could control noise variable by set betai = 0

- rateC:

  censoring rate,

- seed:

  seed

- method:

  "InverseCDF" or "Weibull"

## Value

data frame with columns: Id, HazardWOBaseline, Time, Censor, X, Status,
x1, x2, ...

## Examples

``` r
library(survival)
set.seed(926)
nsim <- 100
beta_true <- c(0, 0.8, -0.8, -1.5, 0, 3, 0, -1, 0.8, 1.5, 3, 0, -0.8, 0, -3)
betaHat <- matrix(nrow = nsim, ncol = 15)
for (k in 1:nsim) {
  df_surv <- surv_sim(
    n_patients = 1000,
    scale_weibull = 1,
    shape_weibull = 1,
    beta = beta_true,
    rateC = 0.3,
    n_vars = 15,
    vars_cate = c(rep("binary", 7), rep("continous", 8)),
    method = "InverseCDF",
    seed = 926 + k
  )
  fit <- coxph(Surv(X, Status) ~ ., data = df_surv[, -c(1:4)])
  betaHat[k, ] <- fit$coef
}
colMeans(betaHat) - beta_true
#>  [1] -0.0034724470  0.0021831449 -0.0074545855 -0.0108587385  0.0062379842
#>  [6]  0.0308821523  0.0025897539 -0.0011369898  0.0153488615  0.0042831387
#> [11]  0.0339854245 -0.0018421230 -0.0087771831  0.0009899273 -0.0399372835
```
