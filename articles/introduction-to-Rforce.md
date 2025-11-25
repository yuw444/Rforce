# introduction-to-Rforce

``` r
library(Rforce)
library(doParallel)
#> Loading required package: foreach
#> Loading required package: iterators
#> Loading required package: parallel
registerDoParallel(cores = 10)
rst <- foreach(i = 1:50) %dopar%{
  data_list <- compo_sim(
    n_patients = 500,
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
df_rst <- t(do.call("cbind", rst))
colMeans(df_rst)
#>  [1]  9.647290e-04 -4.256601e-03 -7.894755e-03  5.941230e-01 -3.709555e-05
#>  [6]  1.744629e-02  7.873844e-01 -9.023901e-03  6.645270e-01 -1.914386e-02
boxplot(df_rst)
```

![](introduction-to-Rforce_files/figure-html/setup-1.png)
