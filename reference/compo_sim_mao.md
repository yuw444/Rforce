# composite event generator for simulation

composite event generator for simulation, this simulator is constructed
based on Mao. L(2015), it assumes constant baseline hazard.

## Usage

``` r
compo_sim_mao(
  n_patients = 1000,
  n_vars = 10,
  vars_cate = c(rep("binary", 6), rep("continuous", 4)),
  true_beta = c(0, 0, 0, 0.6, 0, 0, 0.8, 0, 0.7, 0),
  non_linear_hazard = FALSE,
  non_linear_function = NULL,
  sigma_scale_gamma = 0.25,
  seed = 926,
  s = 2000
)
```

## Arguments

- n_patients:

  number of patients

- n_vars:

  number of covariates

- vars_cate:

  vector of "continuous", "binary"

- true_beta:

  True parameters

- non_linear_hazard:

  whether to use non-linear hazard function

- non_linear_function:

  the non-linear function to generate hazard

- sigma_scale_gamma:

  variance of frailty term

- seed:

  seed for random number generation

- s:

  number of events one patient can up to have, default is 2000

## Value

a list of simulated data and parameters

## Examples

``` r
# example code
library(doParallel)
registerDoParallel(cores = 16)
rst <- foreach(i = 1:48) %dopar%{
  data_list <- compo_sim_mao(
    n_patients = 200,
    seed = i
  )
  library(dplyr)
  df_train <- manual_censoring(data_list[[1]], 0.8)

  estimate_list <- wcompo_est(
    data = df_train,
    weight = c(1, 1)
  )

  estimate_list$beta
}
#> Error in {    data_list <- compo_sim_mao(n_patients = 200, seed = i)    library(dplyr)    df_train <- manual_censoring(data_list[[1]], 0.8)    estimate_list <- wcompo_est(data = df_train, weight = c(1,         1))    estimate_list$beta}: task 1 failed - "could not find function "manual_censoring""
df_rst <- t(do.call("cbind", rst))
#> Error: object 'rst' not found
colMeans(df_rst)
#> Error: object 'df_rst' not found
boxplot(df_rst)
#> Error: object 'df_rst' not found
```
