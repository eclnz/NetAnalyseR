

#' Calculate Normalized Characteristic Path Length
#'
#' @title Normalized Characteristic Path Length Calculation
#' @description Computes the normalized characteristic path length for a given network, which can be provided either as a matrix or within a list alongside randomized matrices.
#' This function first checks if the input is a matrix or a list. If no array of randomized matrices is provided (or not included in the list), one is generated.
#' The function then calculates the average characteristic path length for the random array and normalizes the characteristic path length of the input network by this average.
#' The list elements must not be named.
#' @param mat_list Either a square matrix representing the network for which the characteristic path length is to be calculated, or a list containing the network matrix and an optional 3D array of randomized matrices for comparison.
#' @param rand_array An optional 3D array of randomized matrices for comparison if the input is a single matrix. If NULL and the input is a matrix,
#' a set of randomized matrices is generated. This parameter is ignored if the input is a list.
#' @return The normalized characteristic path length, computed as the ratio of the characteristic path length of the input network to the average characteristic path length of the randomized networks.
#' @examples
#' # Example for a single matrix input
#' W <- matrix(c(0, 2, 1, 0, 2, 0, 3, 5, 1, 3, 0, 6, 0, 5, 6, 0), nrow = 4, byrow = TRUE)
#' norm_char_path_length <- normalised_characteristic_path_length(W)
#' W_rand <- generateRewiredMatrices(W, 100)
#' W_rand <- abind::abind(W_rand, along = 3)
#' norm_char_path_length <- normalised_characteristic_path_length(list(W, W_rand))
#' # Can also be called with the following:
#' norm_char_path_length <- normalised_characteristic_path_length(W, W_rand)
#' @export

normalised_characteristic_path_length <- function(mat_list, rand_array = NULL) {
  if(is.matrix(mat_list)){
    validate_matrix(mat_list)
    if(network_density(mat_list)==1){
      warning("Network density is equal to 1. Network cannot be rewired while maintaining degree distribution\n
              This may make results incorrect.")
    }
    # If the global clustering coefficient is 0 then the normalised clustering coefficient will be 0 and math errors will be created by 0 division.
    if(characteristic_path_length(mat_list)==0){
      return(0)
    }
    # Generate randomized matrices if not provided
    if (is.null(rand_array)) {
      rand_array <- generateRewiredMatrices(mat_list, 100) %>%
        abind::abind(along = 3)
    }
    # Validate that rand_array is indeed an array
    if (!is.array(rand_array)) {
      stop("Random array specified is not in the form of an array")
    }
    # Calculate the average clustering coefficient for the random array
    rand_cpl <- mean(apply(rand_array, MARGIN = 3, FUN = characteristic_path_length))
    # Calculate the clustering coefficient for the input matrix
    cpl <- characteristic_path_length(mat_list)
    # Normalize the clustering coefficient by the random array's average
    norm_cpl <- cpl / rand_cpl
    return(norm_cpl)
  }
  if(is.list(mat_list)){
    if(network_density(mat_list[[1]])==1){
      warning("Network density is equal to 1. Network cannot be rewired while maintaining degree distribution\n
              This may make results incorrect.")
    }
    # If the global clustering coefficient is 0 then the normalised clustering coefficient will be 0 and mat_list[[1]]h errors will be created by 0 division.
    if(characteristic_path_length(mat_list[[1]])==0){
      return(0)
    }
    # Validate that rand_array is indeed an array
    if (!is.array(mat_list[[2]])) {
      stop("Random array specified is not in the form of an array")
    }
    # Calculate the average clustering coefficient for the random array
    rand_cpl <- mean(apply(mat_list[[2]], MARGIN = 3, FUN = characteristic_path_length))
    # Calculate the clustering coefficient for the input matrix
    cpl <- characteristic_path_length(mat_list[[1]])
    # Normalize the clustering coefficient by the random array's average
    norm_cpl <- cpl / rand_cpl
    return(norm_cpl)
  }
}
