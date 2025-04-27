# test_that("rforce works", {
  # Developing help function
  library(dplyr)
  library(sjmisc)
  library(Rforce)

  data <- readRDS(file = "/home/yu89975/r-dev/Rforce/data/test_data.rds") %>% filter(Id < 1000)
  units_of_cpius <- diff(c(0, quantile(data$X, 1 / 10 * 1:10)))
  table(data$Status)

  data_to_dummy <- data %>%
    dplyr::select(-c(X, Status, Id))

  lst <- sapply(data_to_dummy, function(x) {
    if (is.factor(x)) {
      sjmisc::to_dummy(x)[, -1]
    } else {
      x
    }
  })

  data_to_convert <- cbind.data.frame(
    do.call("cbind.data.frame", lst),
    data[, c("Id", "X", "Status")]
  )

  variable_Ids <- colnames(do.call("cbind.data.frame", lst))

  variable_Ids <- gsub("\\.x_\\d+", "", variable_Ids[])

  unique_vars <- unique(variable_Ids)

  variable_Ids <- match(variable_Ids, unique_vars) - 1
  
  lst_cpiu_wide <- patients_to_cpius(
    data_to_convert = data_to_convert,
    units_of_cpiu = units_of_cpius,
    weights_by_status = c(0, 1, 1, 1, 1),
    pseudo_risk = TRUE,
    wide_format = TRUE
  )

  design_matrix_Y <- as.matrix(lst_cpiu_wide$designMatrix_Y)
  auxiliary_features <- as.matrix(lst_cpiu_wide$auxiliaryFeatures)

  design_matrix_Y[is.na(design_matrix_Y)] <- 999999.0

  temp <- Rforce(
    design_matrix_Y = design_matrix_Y,
    auxiliary_features = auxiliary_features,
    variable_Ids = variable_Ids,
    units_of_cpius = units_of_cpius,
    split_rule = "Rforce-QLR",
    n_trees = 4,
    mtry = 3,
    n_splits = 2,
    seed = 926
  )

# })
