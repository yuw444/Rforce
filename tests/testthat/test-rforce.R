test_that("Rforce fits with formula API and returns correct structure", {
  data_list <- compo_sim(n_patients = 100, seed = 42, verbose = FALSE)
  df_train <- random_censoring(data_list$dataset, 0.8)

  covariate_cols <- setdiff(colnames(df_train), c("Id", "X", "Status"))
  formula <- as.formula(paste(
    "Surv(Id, X, Status) ~",
    paste(covariate_cols, collapse = " + ")
  ))

  forest <- Rforce(
    data = df_train,
    formula = formula,
    n_intervals = 5,
    n_trees = 10,
    mtry = 3,
    n_splits = 3,
    max_depth = 3,
    min_node_size = 5,
    min_gain = 0,
    split_rule = "Rforce-QLR",
    seed = 926
  )

  expect_s3_class(forest, "Rforce")
  expect_equal(forest$nTrees, 10)
  expect_equal(nrow(forest$predicted), nrow(forest$cpius$designMatrix_Y))
  expect_equal(ncol(forest$predicted), 5)
})

test_that("predict.Rforce returns matrix of correct dimensions", {
  data_list <- compo_sim(n_patients = 100, seed = 42, verbose = FALSE)
  df_train <- random_censoring(data_list$dataset, 0.8)

  covariate_cols <- setdiff(colnames(df_train), c("Id", "X", "Status"))
  formula <- as.formula(paste(
    "Surv(Id, X, Status) ~",
    paste(covariate_cols, collapse = " + ")
  ))

  forest <- Rforce(
    data = df_train,
    formula = formula,
    n_intervals = 5,
    n_trees = 10,
    mtry = 3,
    n_splits = 3,
    max_depth = 3,
    min_node_size = 5,
    min_gain = 0,
    split_rule = "Rforce-QLR",
    seed = 926
  )

  # in-bag predictions (no new data)
  pred_inbag <- predict(forest)
  expect_true(is.matrix(pred_inbag) || is.data.frame(pred_inbag))

  # predict on design matrix
  design_mat <- as.matrix(forest$cpius$designMatrix_Y[
    1:5,
    !startsWith(colnames(forest$cpius$designMatrix_Y), "nEvents")
  ])
  pred_new <- predict(forest, design_mat)
  expect_equal(nrow(pred_new), 5)
  expect_equal(ncol(pred_new), 5)
  expect_true(all(pred_new >= 0))
})

test_that("vimp.Rforce returns a named numeric vector", {
  data_list <- compo_sim(n_patients = 100, seed = 42, verbose = FALSE)
  df_train <- random_censoring(data_list$dataset, 0.8)

  covariate_cols <- setdiff(colnames(df_train), c("Id", "X", "Status"))
  formula <- as.formula(paste(
    "Surv(Id, X, Status) ~",
    paste(covariate_cols, collapse = " + ")
  ))

  forest <- Rforce(
    data = df_train,
    formula = formula,
    n_intervals = 5,
    n_trees = 10,
    mtry = 3,
    n_splits = 3,
    max_depth = 3,
    min_node_size = 5,
    min_gain = 0,
    split_rule = "Rforce-QLR",
    seed = 926
  )

  imp <- vimp(forest)
  expect_true(is.numeric(imp))
  expect_true(!is.null(names(imp)))
  expect_equal(length(imp), length(covariate_cols))
  expect_true(all(imp >= 0 & imp <= 1))
})

test_that("saveRforce and loadRforce round-trip preserves predictions", {
  data_list <- compo_sim(n_patients = 100, seed = 42, verbose = FALSE)
  df_train <- random_censoring(data_list$dataset, 0.8)

  covariate_cols <- setdiff(colnames(df_train), c("Id", "X", "Status"))
  formula <- as.formula(paste(
    "Surv(Id, X, Status) ~",
    paste(covariate_cols, collapse = " + ")
  ))

  forest <- Rforce(
    data = df_train,
    formula = formula,
    n_intervals = 5,
    n_trees = 10,
    mtry = 3,
    n_splits = 3,
    max_depth = 3,
    min_node_size = 5,
    min_gain = 0,
    split_rule = "Rforce-QLR",
    seed = 926
  )

  out_dir <- tempfile(pattern = "rforce_test_")
  dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)
  saveRforce(forest, out_dir)

  forest2 <- loadRforce(out_dir)
  expect_s3_class(forest2, "Rforce")
  expect_equal(forest2$nTrees, forest$nTrees)

  design_mat <- as.matrix(forest$cpius$designMatrix_Y[
    1:3,
    !startsWith(colnames(forest$cpius$designMatrix_Y), "nEvents")
  ])
  pred1 <- predict(forest, design_mat)
  pred2 <- predict(forest2, design_mat)
  expect_equal(pred1, pred2, tolerance = 1e-3)
})
