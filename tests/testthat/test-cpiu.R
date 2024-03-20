test_that("cpiu", {
  data_list <- compo_sim(
    n_patients = 200,
    verbose = FALSE
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

  expect_equal(sum(df_train$Status ==1 )/200, 0.73)

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

  df_cpiu_wide_nowt <- patients_to_cpius(
    df_train,
    length_cpius,
    c(0, 1, 1),
    pseudo_risk = FALSE,
    wide_format = TRUE
  )

  expect_equal(as.numeric(df_cpiu_wide_nowt$auxiliaryFeatures[2, 5+2]),0)

})
