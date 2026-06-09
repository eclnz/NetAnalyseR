#' @title Calculate Missing Weights in a Network
#'
#' @description Calculates the missing weights in a network by first validating the
#' input matrix, computing both inter-node and intra-node weights, summing these
#' weights to find the total (rounded to one significant figure), and finally
#' determining the missing weights by subtracting the sum of inter and intra-node
#' weights from the total.
#' @param W A matrix representing the network, where rows and columns correspond to
#' nodes, and the matrix elements represent the weight of the connections between
#' these nodes.
#' @return The calculated missing weights in the network.
#' @examples
#' # Example matrix representing weighted connections between nodes
#' W <- matrix(c(0, 2, 1, 4, 2, 0, 3, 5, 1, 3, 0, 6, 4, 5, 6, 0), nrow = 4, byrow = TRUE)
#' missing_weights(W)
#' @export
#'
missing_weights_ <- function(W) {
  inter <- inter_node_(W)
  intra <- intra_node_(W)
  total <- signif(inter + intra, 1)
  total - (inter + intra)
}

missing_weights <- function(W) {
  missing_weights_(validate_matrix(W))
}
