test_that("cpiu", {
  data_list <- compo_sim(
    n_patients = 200,
    verbose = FALSE
  )

  df_train <- manual_censoring(
    data_list$dataset,
    0.99
  )

  df_train <- random_censoring(
    data_list$dataset,
    0.9,
    verbose = TRUE
  )

  sum(df_train$Status ==1 )/200

  hist(df_train$X)

  length_cpius <- diff(quantile(c(0,df_train$X, max(df_train$X) + 0.01), probs = seq(0, 1, length.out = 10 + 1)))

  df_cpiu_wide <- patients_to_cpius(
    df_train,
    length_cpius,
    c(0, 1, 1),
    pseudo_risk = TRUE,
    wide_format = TRUE
  )

  df_cpiu_wide_nowt <- patients_to_cpius(
    df_train,
    length_cpius,
    c(0, 1, 1),
    pseudo_risk = FALSE,
    wide_format = TRUE
  )


})
