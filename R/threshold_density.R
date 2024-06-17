#' @title Threshold Density Function
#' @description Applies a threshold to a matrix to achieve a specified network density.
#'
#' @param mat A numeric matrix representing the network.
#' @param target_density A numeric value indicating the target density.
#'
#' @return A list containing the thresholded matrix (`matrix`), the applied threshold (`threshold`), and the resulting density (`density`).
#'
#' @examples
#' mat <- mat_gen(100)
#' target_density <- 0.2
#' result <- threshold_density(mat, target_density)
#' image(result)
#'
#' @export
threshold_density <- function(mat, target_density) {
  # Validate input matrix
  validate_matrix(mat)

  # Apply threshold to matrix
  threshold_mat <- function(mat, threshold) {
    return(as.numeric(mat > threshold) * mat)
  }

  # Calculate network density
  network_density <- function(mat) {
    num_edges <- sum(mat > 0)
    num_possible_edges <- length(mat) - nrow(mat)  # for undirected networks
    return(num_edges / num_possible_edges)
  }

  # Initial setup
  min_val <- min(mat)
  max_val <- max(mat)
  current_threshold <- (min_val + max_val) / 2
  step_size <- (max_val - min_val) / 4
  tolerance <- 0.001

  # Binary search for optimal threshold
  while (step_size > tolerance) {
    mat_thr <- threshold_mat(mat, current_threshold)
    current_density <- network_density(mat_thr)

    density_difference <- current_density - target_density

    if (abs(density_difference) <= tolerance) {
      break
    }

    if (density_difference > 0) {
      current_threshold <- current_threshold + step_size
    } else {
      current_threshold <- current_threshold - step_size
    }

    step_size <- step_size / 2
  }
  return(mat_thr)
  # return(list(matrix = mat_thr, threshold = current_threshold, density = current_density))
}
