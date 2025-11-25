# composite event generator for simulation

composite event generator for simulation by non-homogeneous Poisson
process

## Usage

``` r
compo_sim(
  n_patients = 1000,
  n_vars = 10,
  vars_cate = c(rep("binary", 6), rep("continuous", 4)),
  true_beta = c(0, 0, 0, 0.6, 0, 0, 0.8, 0, 0.7, 0),
  lambda = 0.01,
  a_shape_weibull = 1,
  sigma_scale_weibull = 1,
  sigma_scale_gamma = 0.05,
  non_linear_hazard = FALSE,
  non_linear_function = NULL,
  seed = 926,
  verbose = FALSE
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

- lambda:

  baseline death hazard for entire population

- a_shape_weibull:

  shape parameter of Weibull distribution

- sigma_scale_weibull:

  scale parameter of Weibull distribution

- sigma_scale_gamma:

  scale parameter of gamma distribution

- non_linear_hazard:

  whether to use non-linear hazard function

- non_linear_function:

  the non-linear function to generate hazard

- seed:

  seed for random number generation

- verbose:

  whether to print out summary of the number of recurrent events per
  patients

## Value

a list of simulated data and parameters

## Examples

``` r
# example
library(doParallel)
#> Loading required package: foreach
#> Loading required package: iterators
#> Loading required package: parallel
registerDoParallel(cores = 16)
rst <- foreach(i = 1:48) %dopar%{
  data_list <- compo_sim(
    n_patients = 200,
    seed = i,
    verbose = FALSE
  )
  dim(data_list[[1]])
  library(dplyr)
  df_train <- data_list[[1]] %>%
    dplyr::mutate(X = Time) %>%
    dplyr::select(-c("Time"))
  estimate_list <- wcompo_est(
    data = df_train,
    weight = c(1, 1)
  )
  estimate_list$beta
}
#> Error in {    data_list <- compo_sim(n_patients = 200, seed = i, verbose = FALSE)    dim(data_list[[1]])    library(dplyr)    df_train <- data_list[[1]] %>% dplyr::mutate(X = Time) %>%         dplyr::select(-c("Time"))    estimate_list <- wcompo_est(data = df_train, weight = c(1,         1))    estimate_list$beta}: task 37 failed - "zero non-NA points"
df_rst <- t(do.call("cbind", rst))
#> Error: object 'rst' not found
colMeans(df_rst)
#> Error: object 'df_rst' not found
boxplot(df_rst)
#> Error: object 'df_rst' not found
```
