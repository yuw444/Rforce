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