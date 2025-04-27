#' Cummax
#' @useDynLib Rforce, .registration=TRUE
#' @export
#'
RC_Cummax <- function(x) {
  .Call("R_Cummax", x)
}

#' ColsPermute
#' @useDynLib Rforce, .registration=TRUE
#' @export
#' @param x A matrix, it is stored in column-major order when passed to C
RC_ColsPermute <- function(x, colsToPermute, seed) {
  x_dim <- dim(x)
  x_new <- matrix(as.numeric(x), nrow = x_dim[1], ncol = x_dim[2])
  temp <- .Call("R_ColsPermute", x_new, as.integer(colsToPermute - 1), as.integer(seed))
  return(temp)
}

#' Sum
#' @useDynLib Rforce, .registration=TRUE
#' @export
#' @param x A vector of numbers
#' @return The sum of the numbers in x
#'
RC_Sum <- function(x, nthreads = 4) {
  .Call("R_Sum", as.numeric(x), as.integer(nthreads))
}


#' matrix add
#' @useDynLib Rforce, .registration=TRUE
#' @export
#' @param x A matrix
#' @param y A matrix
#' @return The sum of the two matrices
#'
RC_MatrixAdd <- function(x, y) {
  x_new <- matrix(as.numeric(x), nrow = dim(x)[1], ncol = dim(x)[2])
  y_new <- matrix(as.numeric(y), nrow = dim(y)[1], ncol = dim(y)[2])
  .Call("R_MatrixAdd", x_new, y_new)
}

#' return a list to R
#' @useDynLib Rforce, .registration=TRUE
#' @export
#' @return A list
#'
RC_ReturnList <- function() {
  .Call("R_ListOfVectors")
}


#' Main Rforce function
#' @useDynLib Rforce, .registration=TRUE
#' @export
#' @importFrom sjmisc to_dummy
#' @param design_matrix_Y data.frame that includes design matrix and response(number of events at each CPIU)
#' @param auxiliary_features data.frame that includes Id, X, Status, pseudo-risktime-at-duration
#'

Rforce <- function(
    design_matrix_Y,
    auxiliary_features,
    variable_Ids,
    units_of_cpius,
    split_rule = c("Rforce-QLR", "Rforce-GEE", "Rforce-GEEIntr", "RF-SLAM"),
    n_trees = 500,
    max_depth = 8,
    min_node_size = 10,
    min_gain = 0,
    mtry = NA,
    n_splits = NA,
    seed = 926) {
  ## Developing help function
  # data <- readRDS(file = "/home/yu89975/r-dev/Rforce/data/test_data.rds")
  # units_of_cpius <- diff(c(0, quantile(data$X[data$Status != 0], 1 / 4 * 1:4)))
  # table(data$Status)

  # data_to_dummy <- data %>%
  #   dplyr::select(-c(X, Status, Id))

  # lst <- sapply(data_to_dummy, function(x) {
  #   if (is.factor(x)) {
  #     sjmisc::to_dummy(x)[, -1]
  #   } else {
  #     x
  #   }
  # })

  # data_to_convert <- cbind.data.frame(
  #   do.call("cbind.data.frame", lst),
  #   data[, c("Id", "X", "Status")]
  # )

  # variable_Ids <- colnames(do.call("cbind.data.frame", lst))

  # variable_Ids <- gsub("\\.x_\\d+", "", variable_Ids[])

  # unique_vars <- unique(variable_Ids)

  # variable_Ids <- match(variable_Ids, unique_vars) - 1
  
  # lst_cpiu_wide <- patients_to_cpius(
  #   data_to_convert = data,
  #   units_of_cpiu = units_of_cpius,
  #   weights_by_status = c(0, 1, 1, 1, 1),
  #   pseudo_risk = TRUE,
  #   wide_format = TRUE
  # )

  # design_matrix_Y <- lst_cpiu_wide$designMatrix_Y
  # auxiliary_features <- lst_cpiu_wide$auxiliaryFeatures

  ## choose splitting rule
  split_rule_index <- switch(split_rule,
    "Rforce-QLR" = 0,
    "Rforce-GEE" = 1,
    "Rforce-GEEIntr" = 1,
    "RF-SLAM" = 2,
    -1
  )

  if (split_rule_index == -1) {
    stop('Invalid split_rule, must choose from c("Rforce-QLR", "Rforce-GEE", "Rforce-GEEIntr", "RF-SLAM")')
  }

  gee_interaction <- switch(split_rule,
    "Rforce-GEEIntr" = 1,
    0
  )

  .Call(
    "R_Rforce",
    as.integer(split_rule_index),
    as.integer(gee_interaction),
    matrix(as.numeric(design_matrix_Y), nrow = dim(design_matrix_Y)[1], ncol = dim(design_matrix_Y)[2]),
    matrix(as.numeric(auxiliary_features), nrow = dim(auxiliary_features)[1], ncol = dim(auxiliary_features)[2]),
    as.integer(variable_Ids),
    as.numeric(units_of_cpius),
    as.integer(n_trees),
    as.integer(max_depth),
    as.integer(min_node_size),
    as.numeric(min_gain),
    as.integer(mtry),
    as.integer(n_splits),
    as.integer(seed)
  )

}
