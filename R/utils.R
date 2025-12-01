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

vimp <- function(object, ...) {
  UseMethod("vimp")
}

#' Variable importance function for Rforce object
#' @export
#' @param forest An Rforce object
#' @return A matrix containing variable importance statistics for each variable in each tree of the forest,
#' row correspond to trees and columns correspond to variables.
vimp.Rforce <- function(forest) {
  if (!inherits(forest, "Rforce")) {
    stop("forest must be an Rforce object.")
  }

  vec <- colMeans(forest$vimpStat)
  vec <- (vec - min(vec)) / (max(vec) - min(vec))
  names(vec) <- forest$cpius$variableUsed
  return(vec)
}


#' Save an Rforce Object to Disk
#'
#' \code{saveRforce()} saves a trained \code{Rforce} object to disk, including both
#' the R-side metadata (e.g., model settings, variable importance, parameters) and
#' the full C-based forest structure. The external C pointer is handled separately
#' and excluded from the saved \code{.rds} file to ensure portability and session
#' safety.
#'
#' @export
#'
#' @param forest An object of class \code{"Rforce"} representing the fitted random
#'   survival forest.
#' @param file A character string specifying the directory path where the object
#'   should be saved. The directory will be created by the underlying C routine
#'   if it does not already exist.
#' @param ... Additional arguments passed to \code{saveRDS()}, such as
#'   \code{compress} or \code{version}.
#'
#' @details
#' The function performs the following actions:
#' \enumerate{
#'   \item Calls the native routine \code{R_SaveRforce} to serialize and write the
#'         forest structure (all trees) to the specified directory.
#'   \item Removes the external C pointer (\code{_external_forest_C_Ptr}) from the
#'         R object to prevent session-specific pointer issues.
#'   \item Saves the remaining R-side object to \code{Rforce.rds} using
#'         \code{saveRDS()}.
#' }
#'
#' The resulting directory structure is compatible with \code{loadRforce()}, which
#' restores both the R object and its native forest pointer.
#'
#' @return
#' Invisibly returns \code{NULL}. A message is not printed, but files are written
#' to the specified directory.
#'
#' @seealso
#' \code{\link{loadRforce}} to reload a saved Rforce object.
#'
#' @examples
#' \dontrun{
#' # Fit an Rforce model and then save it
#' saveRforce(rf_model, file = "saved_model")
#' }
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

#' Load an Rforce Object from Disk
#'
#' \code{loadRforce()} restores a previously saved \code{Rforce} object from a
#' specified directory. It reconstructs both the R-side object (stored as
#' \code{Rforce.rds}) and the associated C++ forest structure via the
#' \code{R_LoadRforce} routine, ensuring full functionality.
#'
#' @useDynLib Rforce, .registration=TRUE
#' @export
#'
#' @param path A character string specifying the directory containing the saved
#'   \code{Rforce} object. The directory must include \code{Rforce.rds} and a
#'   single forest data subdirectory created by \code{saveRforce()}.
#'
#' @details
#' The function performs the following tasks:
#' \enumerate{
#'   \item Reads the \code{Rforce.rds} file to reconstruct the R-side object.
#'   \item Validates that the restored object inherits class \code{"Rforce"}.
#'   \item Locates the associated forest data directory.
#'   \item Calls the native routine \code{R_LoadRforce} to rebuild the C-allocated
#'         forest structure, restoring the external pointer at
#'         \code{_external_forest_C_Ptr}.
#' }
#'
#' If the directory structure is invalid or the C pointer fails to load, an
#' error is raised.
#'
#' @return A fully functional \code{Rforce} object with a reconstructed external
#'   forest pointer, ready for prediction, visualization, or further analysis.
#'
#' @examples
#' \dontrun{
#' # Load a previously saved Rforce model
#' rf <- loadRforce("path/to/saved/model")
#' rf
#' }

loadRforce <- function(path) {
  if (!is.character(path) || length(path) != 1) {
    stop("file must be a single character string.")
  }

  forest <- readRDS(paste0(path, "/Rforce.rds"))
  if (!inherits(forest, "Rforce")) {
    stop("The loaded object is not an Rforce object.")
  }

  forest_path <- list.dirs(path, full.names = TRUE, recursive = FALSE)

  if (length(forest_path) != 1) {
    stop(
      "The specified path does not contain a valid Rforce directory structure."
    )
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

#' Visualize a Single Tree from an Rforce Forest
#'
#' This function extracts a single decision tree from an `Rforce` forest,
#' generates its Graphviz DOT representation via the underlying C backend,
#' and renders it using \pkg{DiagrammeR}. The DOT file is stored temporarily
#' and automatically removed after the R session.
#'
#' @param forest An object of class \code{Rforce} containing the trained forest
#'   and a valid external C pointer in \code{`_external_forest_C_Ptr`}.
#' @param treeIndex An integer specifying the index of the tree to visualize.
#'   Must be between \code{1} and \code{forest\$numTrees}.
#'
#' @details
#' Internally, this function:
#' \enumerate{
#'   \item Validates input.
#'   \item Creates a temporary \code{.dot} file.
#'   \item Calls the C function \code{R_PrintTree} via \code{.Call} to write
#'         the DOT graph representation of the tree.
#'   \item Displays the tree structure using \code{DiagrammeR::grViz()}.
#' }
#'
#' The DOT file persists only for the current R session and is stored in a
#' temporary directory. It is automatically cleaned up when the session ends.
#'
#' @return
#' A DiagrammeR graph object representing the visualized decision tree.
#' Additionally, a message is printed showing the location of the temporary
#' DOT file.
#'
#' @importFrom DiagrammeR grViz
#'
#' @examples
#' \dontrun{
#' # Assuming `rf_model` is an Rforce object with at least one tree:
#' printTree(rf_model, treeIndex = 1)
#' }
#'
#' @export

printTree <- function(forest, treeIndex) {
  if (!inherits(forest, "Rforce")) {
    stop("forest must be an Rforce object.")
  }

  if (!is.numeric(treeIndex) || length(treeIndex) != 1 ||
      treeIndex < 1 || treeIndex > forest$nTrees) {
    stop("treeIndex must be a valid tree index within the forest.")
  }

  file_name <- tempfile(fileext = ".dot")

  .Call(
    "R_PrintTree",
    forest$`_external_forest_C_Ptr`,
    as.integer(treeIndex),
    as.character(file_name)
  )

  message(paste("Tree structure saved to", file_name))

  DiagrammeR::grViz(readLines(file_name))

}