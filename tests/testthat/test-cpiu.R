test_that("cpiu", {

  data_list <- compo_sim(
    n_patients = 200,
    verbose = FALSE
  )

  df_censor1 <- random_censoring(
    data_list[[1]],
    0.5
  )

  df_censor2 <- random_censoring(
    data_list[[1]],
    0.5
  )

  identical(df_censor1, df_censor2)

  temp <- empirical_Y(
    data_list,
    x = rep(1, 10),
    time_points = c(1, 2, 3),
    n_sims = 3,
    n_size = 2000
  )

  data_list_mao <- compo_sim_mao(
    n_patients = 200
  )

  temp_mao <- empirical_Y(
    data_list_mao,
    weibull_baseline = FALSE,
    x = rep(1, 10),
    time_points = c(1, 2, 3),
    n_sims = 3,
    n_size = 200
  )

  temp <- true_Y_numerical_form(
    3,
    a_shape_weibull = data_list$a_shape_weibull,
    sigma_scale_weibull = data_list$sigma_scale_weibull,
    sigma_scale_gamma = data_list$sigma_scale_gamma,
    lambda = data_list$lambda,
    lambdaZ = data_list$hazard_function(rep(1,10))[1,1]
  )

  temp2 <- true_Y_numerical_form(
    1,
    constant_baseline_hazard = TRUE,
    sigma_scale_gamma = data_list_mao$sigma_scale_gamma,
    lambda = data_list_mao$lambda,
    lambdaZ = data_list_mao$hazard_function(rep(1,10))[1,1]
  )

  df_train_admin <- admin_censoring(
    data_list$dataset,
    0.99
  )

  expect_warning(df_train <- random_censoring(
    data_list$dataset,
    0.9,
    verbose = TRUE
  ))

  expect_equal(sum(df_train$Status ==1 )/200, 0.86)

  df_train <- random_censoring(
    data_list$dataset,
    0.73,
    verbose = TRUE
  )

  expect_equal(sum(df_train$Status ==1 )/200, 0.795)

  hist(df_train$X)

  length_cpius <- diff(quantile(c(0,df_train$X, max(df_train$X) + 0.01), probs = seq(0, 1, length.out = 10 + 1)))

  expect_gte(sum(length_cpius), max(df_train$X))

  df_cpiu_wide <- patients_to_cpius(
    df_train,
    length_cpius,
    c(0, 1, 1),
    pseudo_risk = TRUE,
    wide_format = TRUE
  )

  expect_gt(as.numeric(df_cpiu_wide$auxiliaryFeatures[2, 5+2]), 0)

  df_cpiu_wide_nort <- patients_to_cpius(
    df_train,
    length_cpius,
    c(0, 1, 1),
    pseudo_risk = FALSE,
    wide_format = TRUE
  )

  expect_equal(as.numeric(df_cpiu_wide_nowt$auxiliaryFeatures[2, 5+2]),0)

})
