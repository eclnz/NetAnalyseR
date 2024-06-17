#' @title Validate a Square Numeric Matrix
#'
#' @description This function checks if the input is a valid square numeric matrix without NA values,
#' ensuring it is symmetric. If the matrix is not symmetric, it modifies the matrix
#' to make it symmetric by mirroring the lower triangular part onto the upper triangular part.
#'
#' @param W A matrix to validate.
#' @return Returns a valid symmetric square numeric matrix.
#' If the matrix fails any validation checks, the function stops with an appropriate error message.
#'
validate_matrix <- function(W) {
  # Check if input is a matrix
  if (!is.matrix(W)) {stop("Input must be a matrix.")}

  # Check if all entries in the matrix are numeric
  if (!all(is.numeric(W))) {stop("All entries in matrix W must be numeric.")}

  # Check for NA values in the matrix
  if (any(is.na(W))) {stop("Matrix contains NA values, which are not allowed.")}

  # Check for empty matrix
  if (all(W[upper.tri(W, diag = TRUE)]==0)) {stop("Matrix must have connection values greater than zero.")}

  # Check if the matrix is square
  if (nrow(W) != ncol(W)) {stop("Matrix W must be square.")}

  # Check if the matrix has more than one node (dimension greater than 1x1)
  if(any(dim(W)==1)) {stop("Matrix must have more than one node")}

  # Ensure the matrix is symmetric. If not, modify it to be symmetric
  if (any(t(W) != W)) {
    # For a non-symmetric matrix, copy the lower triangular part to the upper triangular part
    W[upper.tri(W)] <- t(W)[upper.tri(W)]
    # Warn the user about the modification
    warning("Matrix is not symmetric. Setting upper half to match lower half.")
  }

  # Return the validated (and possibly modified) matrix
  return(W)
}
