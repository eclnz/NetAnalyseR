#' @title Convert Matrix Weights to Lengths
#'
#' @description Converts the weights in a numeric matrix to lengths based on the principle that length is
#' inversely proportional to the connection strength. This function first validates the input matrix,
#' sets its diagonal elements to zero to avoid self-connections, and then converts all non-zero
#' weights to their reciprocal values, effectively transforming weights into lengths.
#' @param W A numeric matrix whose weights represent connection strengths.
#' @param validate Whether to validate the input matrix.
#' @return A matrix where each non-zero weight has been converted to its reciprocal value, representing the length.
#' @examples
#' W <- matrix(c(0, 2, 1, 4, 2, 0, 3, 5, 1, 3, 0, 6, 4, 5, 6, 0), nrow = 4, byrow = TRUE)
#' lengths_W <- length_inversion(W)
#'
#' @export
#'
length_inversion <- function(W, validate = TRUE){
  if(validate){W <- validate_matrix(W)} # Check if the matrix is valid
  diag(W) <- 0 # Set diagonal to zero

  # Converts weights to lengths for all non-zero elements.
  E <- which(W != 0, arr.ind = TRUE) # Identify non-zero elements
  W[E] <- 1 / W[E] # Convert weights to lengths

return(W)
}

