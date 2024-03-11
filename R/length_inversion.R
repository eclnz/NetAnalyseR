#' Convert Matrix Weights to Lengths
#'
#' Converts the weights in a numeric matrix to lengths based on the principle that length is
#' inversely proportional to the connection strength. This function first validates the input matrix,
#' sets its diagonal elements to zero to avoid self-connections, and then converts all non-zero
#' weights to their reciprocal values, effectively transforming weights into lengths.
#'
#' @param W A numeric matrix whose weights represent connection strengths.
#' @param wcm Unused in this version of the function, included for compatibility with previous versions.
#'
#' @return A matrix where each non-zero weight has been converted to its reciprocal value, representing the length.
#'
#' @examples
#' W <- matrix(c(0, 2, 0, 0.5, 2, 0, 3, 0, 0, 0.5, 3, 0, 0, 0, 0, 0), 4, 4)
#' lengths_W <- weight_conversion(W, 'lengths')
#'
#' Validates the input matrix and ensures no self-connections by setting diagonal elements to zero.
#' @export
#'
length_inversion <- function(W){
  W <- validate_matrix(W) # Check if the matrix is valid
  diag(W) <- 0 # Set diagonal to zero

  # Converts weights to lengths for all non-zero elements.
  E <- which(W != 0, arr.ind = TRUE) # Identify non-zero elements
  W[E] <- 1 / W[E] # Convert weights to lengths

return(W)
}

