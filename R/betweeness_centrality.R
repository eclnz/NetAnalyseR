#' Calculate Betweenness Centrality for Weighted Networks
#'
#' @description This function computes the betweenness centrality for each node in a weighted network represented by a matrix.
#' Betweenness centrality measures the extent to which a node lies on paths between other nodes.
#' Nodes with higher betweenness centrality play a more significant role in the communication/interaction within the network.
#' The function assumes that a zero in the matrix represents no connection and infinities represent disconnected paths.
#' @param W A square numeric matrix representing the weighted network, where higher weights imply stronger connections.
#' @param validate Whether to validate the input matrix.
#' @return A numeric vector containing the betweenness centrality for each node in the network.
#' @examples
#' W<- matrix(c(0, 2, 1, 0, 2, 0, 3, 5, 1, 3, 0, 6, 0, 5, 6, 0), nrow = 4, byrow = TRUE)
#' betweenness_wei(W)
#' @export
#'
betweenness_wei <- function(W, validate = TRUE) {
  if(validate){validate_matrix(W)}
  diag(W) <- 0
  L <- length_inversion(W, FALSE)

  # Get shortest distances and number of paths using predefined functions
  D <- shortest_distance(L, FALSE)
  NP <- shortest_distance((L>0)+0, FALSE)
  # Call C++ function
  BC <-calculateBetweennessCentrality(D,NP)
  return(BC)
}
