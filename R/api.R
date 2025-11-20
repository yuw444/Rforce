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
#' @importFrom fastDummies dummy_cols
#' @param design_matrix_Y data.frame that includes design matrix and response(number of events at each CPIU)
#' @param auxiliary_features data.frame that includes Id, X, Status, pseudo-risktime-at-duration
#'

# ...existing code...
Rforce <- function(
    data = NULL,
    formula = NULL, # e.g. Surv(X, Status) ~ covariate1 + covariate2
    id = "Id",
    n_intervals = 10,
    design_matrix_Y = NULL,
    auxiliary_features = NULL,
    variable_Ids = NULL,
    units_of_cpius = NULL,
    event_types = NULL,
    event_weights = NULL, # e.g. c(0,1,1,1,1)
    split_rule = c("Rforce-QLR", "Rforce-GEE", "Rforce-GEEIntr", "RF-SLAM"),
    n_trees = NULL,
    max_depth = NULL,
    min_node_size = NULL,
    min_gain = NULL,
    mtry = NA,
    n_splits = NA,
    seed = 926) {
  ## If user provided raw data + formula, build design_matrix_Y and auxiliary_features here
  if (is.null(design_matrix_Y) || is.null(auxiliary_features)) {
    if (is.null(data) || is.null(formula)) {
      stop("If design_matrix_Y or auxiliary_features are not provided you must supply 'data' and a 'formula' (e.g. Surv(time, status) ~ covariates).")
    }

    # Extract time and status variable names from left hand side of the formula
    lhs_vars <- all.vars(formula[[2]])
    if (length(lhs_vars) != 3) {
      stop("Formula left hand side must be a Surv(id, x, status) with patient id, event time, and status variables, e.g. Surv(id, time, status) ~ covariate_1 + covariate_1.")
    }
    idName <- lhs_vars[1]
    timeName <- lhs_vars[2]
    statusName <- lhs_vars[3]

    if (!all(c(idName, timeName, statusName) %in% colnames(data))) {
      stop(sprintf("Data must include columns '%s', '%s' and '%s'.", idName, timeName, statusName))
    }

    # compute units_of_cpius if not supplied (default: deciles of follow-up time)
    if (is.null(units_of_cpius)) {
      probs <- seq(0, 1, length.out = n_intervals + 1)
      break_points <- stats::quantile(data[[timeName]], probs = probs, na.rm = TRUE)
      break_points[1] <- 0.0 # ensure starting at 0
      break_points[n_intervals + 1] <- max(data[[timeName]], na.rm = TRUE) * 1.001 # ensure ending at max time
      units_of_cpius <- diff(break_points) # compute widths of intervals
    }
    # remove zeros if any (defensive)
    units_of_cpius <- units_of_cpius[units_of_cpius > 0]
    if (length(units_of_cpius) == 0) stop("Computed units_of_cpius is empty; check your time variable.")

    # prepare covariate set (drop Id, time, status)
    cov_names <- setdiff(colnames(data), c(idName, timeName, statusName))
    data_cov <- data[, cov_names, drop = FALSE]

    # dummy-encode factors (keep same approach as tests: remove reference column)
    covariate_df <- fastDummies::dummy_cols(
      data_cov,
      remove_first_dummy = TRUE, # drop reference level per factor
      ignore_na = TRUE,
      remove_selected_columns = TRUE # drop original factor/char columns
    )
    # ensure no rownames and keep as data.frame
    covariate_df <- as.data.frame(covariate_df, stringsAsFactors = FALSE)

    # create variable_Ids mapping if not provided
    if (is.null(variable_Ids)) {
      variable_Ids_names <- colnames(covariate_df)
      # remove internal suffixes that to_dummy may add like ".x_1"
      variable_Ids_names <- gsub("\\_\\d+", "", variable_Ids_names)
      unique_vars <- unique(variable_Ids_names)
      variable_Ids_int <- match(variable_Ids_names, unique_vars) - 1
    } else {
      variable_Ids_int <- as.integer(variable_Ids)
    }

    # build full table expected by patients_to_cpius: covariates + Id + time + status
    data_to_convert <- cbind.data.frame(
      data[, c(idName, timeName, statusName), drop = FALSE],
      covariate_df,
      stringsAsFactors = FALSE
    )

    # make sure status include 0 and 1, or "0", "1"
    if (!all(c(0, 1) %in% unique(data[[statusName]])) &&
      !all(c("0", "1") %in% unique(data[[statusName]]))) {
      stop(sprintf("Status variable '%s' must include at least two event types: 0 (censoring) and 1 (terminal).", statusName))
    }
    # record unique event types and convert to 0-based integer
    event_types <- sort(unique(data[[statusName]]))
    data[[statusName]] <- as.integer(as.factor(data[[statusName]])) - 1

    # event_weights default
    if (is.null(event_weights)) {
      event_weights <- c(0, rep(1, length(event_types) - 1))
    }

    # check event_weights length
    if (length(event_weights) != length(event_types)) {
      stop("Length of event_weights must match number of event types in the data.")
    }

    ## choose splitting rule
    split_rule_index <- switch(split_rule,
      "Rforce-QLR" = 0,
      "Rforce-GEE" = 1,
      "Rforce-GEEIntr" = 1,
      "RF-SLAM" = 2,
      -1
    )

    if (is.null(design_matrix_Y) || is.null(auxiliary_features)) {
      message("Converting data to design matrix and auxiliary features...")
      # Call patients_to_cpius to produce wide-format design matrix and auxiliary features
      lst_cpiu_wide <- patients_to_cpius(
        data_to_convert = data_to_convert,
        units_of_cpiu = units_of_cpius,
        weights_by_status = event_weights,
        pseudo_risk = ifelse(split_rule_index == 2, FALSE, TRUE),
        wide_format = TRUE
      )

      # Extract design matrix and auxiliary features
      design_matrix_Y <- as.matrix(lst_cpiu_wide$designMatrix_Y)
      auxiliary_features <- as.matrix(lst_cpiu_wide$auxiliaryFeatures)
    }

    # NA is handled in C backend
    design_matrix_Y[1:10, 1:10] <- NA

    # store computed variable ids
    variable_Ids <- variable_Ids_int
    # keep units_of_cpius possibly updated above
  } # end internal conversion

  if (split_rule_index == -1) {
    stop('Invalid split_rule, must choose from c("Rforce-QLR", "Rforce-GEE", "Rforce-GEEIntr", "RF-SLAM")')
  }

  gee_interaction <- switch(split_rule,
    "Rforce-GEEIntr" = 1,
    0
  )

  # sink("temp.txt", split = TRUE)
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
  # sink()

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
    `_external_forest_C_Ptr` = temp_C_return[["_external_forest_C_Ptr"]],
    # newly stored fields for convenience / reproducibility
    design_matrix_Y = design_matrix_Y,
    auxiliary_features = auxiliary_features,
    variable_Ids = variable_Ids,
    units_of_cpiu = units_of_cpius,
    data = if (exists("data")) data else NULL,
    formula = if (exists("formula")) formula else NULL,
    event_types = if (exists("event_types")) event_types else NULL,
    event_weights = if (exists("event_weights")) event_weights else NULL
  )

  class(rforce_obj) <- "Rforce"
  return(rforce_obj)
}
