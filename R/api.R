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

  temp_C_return <- .Call(
    "R_Rforce",
    as.integer(split_rule_index),
    as.integer(gee_interaction),
    matrix(as.numeric(design_matrix_Y), nrow = dim(design_matrix_Y)[1], ncol = dim(design_matrix_Y)[2]),
    matrix(as.numeric(auxiliary_features), nrow = dim(auxiliary_features)[1], ncol = dim(auxiliary_features)[2]),
    as.integer(variable_Ids),
    as.integer(length(unique(variable_Ids))),
    as.numeric(units_of_cpius),
    as.integer(n_trees),
    as.integer(max_depth),
    as.integer(min_node_size),
    as.numeric(min_gain),
    as.integer(mtry),
    as.integer(n_splits),
    as.integer(seed)
  )

  # Extract struct values from the returned list
  rforce_obj <- list(
    bagMatrix = temp_C_return[["bagMatrix"]],
    treePhi = temp_C_return[["treePhi"]],
    forestMatrix = temp_C_return[["forestMatrix"]],
    vimpStat = temp_C_return[["vimpStat"]],
    predicted = temp_C_return[["predicted"]],
    oobPredicted = temp_C_return[["oobPredicted"]],
    # You can add more fields here if your C code returns them
    nrowsDesign = nrow(design_matrix_Y),
    ncolsDesign = ncol(design_matrix_Y),
    varIDs = variable_Ids,
    nVars = length(unique(variable_Ids)),
    nTrees = n_trees,
    unitsOfCPIU = units_of_cpius,
    nUnits = length(units_of_cpius),
    lenOutput = ncol(temp_C_return[["predicted"]]),
    maxDepth = max_depth,
    minNodeSize = min_node_size,
    minGain = min_gain,
    mtry = mtry,
    nsplits = n_splits,
    seed = seed,
    `_external_forest_C_Ptr` = temp_C_return[["_external_forest_C_Ptr"]]
  )

  class(rforce_obj) <- "Rforce"
  return(rforce_obj)
}
