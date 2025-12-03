# introduction-to-Rforce

``` r
library(Rforce)
library(doParallel)
library(dplyr)
registerDoParallel(cores = 10)
```

``` r
n_sims <- if (interactive()) 20 else 5
n_patients <- if (interactive()) 500 else 200

rst <- foreach(i = 1:n_sims) %dopar%{
  data_list <- compo_sim(
    n_patients = n_patients,
    seed = i,
    true_beta = c(0, 0, 0, 0.6, 0, 0, 0.8, 0, 0.7, 0),
    verbose = FALSE
  )
  df_train <- data_list[[1]] %>%
    dplyr::mutate(X = Time) %>%
    dplyr::select(-c("Time"))

  estimate_list <- wcompo_est(
    data = df_train,
    weight = c(1, 1)
  )
  estimate_list$beta
}
```

``` r
df_rst <- t(do.call("cbind", rst))
colMeans(df_rst)
#>  [1]  0.02425514 -0.01624627  0.07268047  0.61358282  0.04531289 -0.01801206
#>  [7]  0.77200064  0.11502782  0.63473937  0.14195386
boxplot(df_rst)
```

![The beta estimates from
simulations](introduction-to-Rforce_files/figure-html/unnamed-chunk-3-1.png)
