#' Predict function for Rforce
#' @useDynLib Rforce, .registration=TRUE
#' @export
#' @param forest An Rforce forest object
#' @param designMatrix A matrix of design features for prediction
#' @return A list containing the predicted values
#'
predict.Rforce <- function(forest, designMatrix) {
  if (!inherits(forest, "Rforce")) {
    stop("forest must be an Rforce object.")
  }

  if (!is.data.frame(designMatrix) && !is.matrix(designMatrix)) {
    stop("designMatrix must be a data frame or matrix.")
  }

  .Call(
    "R_ForestPredict",
    forest$`_external_forest_C_Ptr`,
    as.matrix(designMatrix)
  )
}

#' Variable importance function for Rforce object
#' @export
#' @param object An Rforce object
#' @return A matrix containing variable importance statistics for each variable in each tree of the forest, 
#' row correspond to trees and columns correspond to variables.

vimp.Rforce <- function(forest) {
  if (!inherits(forest, "Rforce")) {
    stop("forest must be an Rforce object.")
  }

  return(forest$vimpStat)
  
}


#' saveRDS function for Rforce
#' @useDynLib Rforce, .registration=TRUE
#' @export
#' @param object An Rforce object
#' @param file A character string giving the name of the file where the object should be saved
saveRforce <- function(forest, file, ...) {
  if (!inherits(forest, "Rforce")) {
    stop("object must be an Rforce object.")
  }

  if (!is.character(file) || length(file) != 1) {
    stop("file must be a single character string.")
  }

  .Call(
    "R_SaveRforce",
    forest$`_external_forest_C_Ptr`,
    as.character(file)
  )
  forest1 <- forest
  forest1$`_external_forest_C_Ptr` <- NULL
  saveRDS(forest1, paste0(file, "/Rforce.rds"), ...)
  rm(forest1)
}

#' loadRforce function for Rforce
#' @useDynLib Rforce, .registration=TRUE
#' @export
#' @param path A character string giving the name of the path from which the object should be loaded
loadRforce <- function(path) {
  if (!is.character(path) || length(path) != 1) {
    stop("file must be a single character string.")
  }

  forest <- readRDS(paste0(path, "/Rforce.rds"))
  if (!inherits(forest, "Rforce")) {
    stop("The loaded object is not an Rforce object.")
  }

  forest_path <- list.dirs(path, full.names = TRUE, recursive = FALSE)

  if(length(forest_path) != 1) {
    stop("The specified path does not contain a valid Rforce directory structure.")
  }

  forest$`_external_forest_C_Ptr` <- .Call(
    "R_LoadRforce",
    as.character(forest_path)
  )
  if (is.null(forest$`_external_forest_C_Ptr`)) {
    stop("Failed to load the Rforce object from the specified path.")
  }
  return(forest)
}
