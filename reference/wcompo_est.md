# Fast implementation of Lu Mao's Wcompo method by Kim So Young

This function is a fast implementation of Lu Mao's Wcompo method by Kim
So Young

## Usage

``` r
wcompo_est(data, weight)
```

## Arguments

- data:

  data frame including Id, X, Status, and covariates

- weight:

  weight for death and recurrent events

## Value

list, including beta estimate, y, and t

## Examples

``` r
data_list <- compo_sim(
  n_patients = 500,
  seed = 926,
  verbose = FALSE
)
library(dplyr)
#> 
#> Attaching package: ‘dplyr’
#> The following objects are masked from ‘package:stats’:
#> 
#>     filter, lag
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, setequal, union
df_train <- data_list[[1]] %>%
  dplyr::mutate(X = Time) %>%
  dplyr::select(-c("Time"))

estimate_list <- wcompo_est(
  data = df_train,
  weight = c(1, 1)
)
estimate_list$beta
#>              [,1]
#>  [1,] -0.02843479
#>  [2,]  0.03137461
#>  [3,] -0.06931644
#>  [4,]  0.55227028
#>  [5,]  0.01203906
#>  [6,] -0.02577392
#>  [7,]  0.68532311
#>  [8,] -0.06726258
#>  [9,]  0.89298230
#> [10,]  0.14333866
```
