test_that("compo_sim", {

  library(doParallel)
  registerDoParallel(cores = 16)
  rst <- foreach(i = 1:100) %dopar%{
    data_list <- compo_sim(
      n_patients = 200,
      seed = i,
      verbose = FALSE
    )
    library(dplyr)
    df_train <- random_censoring(
      data_list$dataset,
      0.9
    )
    sum(df_train$Status == 1)/200
    estimate_list <- wcompo_est(
       data = df_train,
       weight = c(1, 1)
    )
    estimate_list$beta
  }
  df_rst <- t(do.call("cbind", rst))
  temp <- colMeans(df_rst)
  expect_equal(temp[4], 0.6, tolerance = 0.05)
  expect_equal(temp[7], 0.8, tolerance = 0.05)
  expect_equal(temp[9], 0.7, tolerance = 0.05)
})


test_that("compo_sim_mao", {

  library(doParallel)
  registerDoParallel(cores = 16)
  rst <- foreach(i = 1:100) %dopar%{
    data_list <- compo_sim_mao(
      n_patients = 200,
      seed = i
    )
    dim(data_list[[1]])
    library(dplyr)
    df_train <- manual_censoring(data_list[[1]], 0.8)

    estimate_list <- wcompo_est(
      data = df_train,
      weight = c(1, 1)
    )

    estimate_list$beta
  }
  df_rst <- t(do.call("cbind", rst))
  temp <- colMeans(df_rst)
  boxplot(df_rst)
  expect_equal(temp[4], 0.6, tolerance = 0.05)
  expect_equal(temp[7], 0.8, tolerance = 0.05)
  expect_equal(temp[9], 0.7, tolerance = 0.05)
})


