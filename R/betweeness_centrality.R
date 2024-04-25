#' Calculate Betweenness Centrality for Weighted Networks
#'
#' @description This function computes the betweenness centrality for each node in a weighted network represented by a matrix.
#' Betweenness centrality measures the extent to which a node lies on paths between other nodes.
#' Nodes with higher betweenness centrality play a more significant role in the communication/interaction within the network.
#' The function assumes that a zero in the matrix represents no connection and infinities represent disconnected paths.
#' @param G A square numeric matrix representing the weighted network, where higher weights imply stronger connections.
#' @return A numeric vector containing the betweenness centrality for each node in the network.
#' @examples
#' G<- matrix(c(0, 2, 1, 0, 2, 0, 3, 5, 1, 3, 0, 6, 0, 5, 6, 0), nrow = 4, byrow = TRUE)
#' betweenness_wei(G)
#' @export
#'
betweenness_wei <- function(G) {
  n <- nrow(G)  # Number of nodes
  BC <- numeric(n)  # Initialize betweenness centrality vector
  diag(G) <- Inf
  G[G == 0] <- Inf  # Assuming 0 means no connection

  # Get shortest distances and number of paths using predefined functions
  D <- shortest_distance(G)
  NP <- shortest_distance((G>1&G<Inf)+0)

  # Calculate betweenness centrality
  for (s in 1:n) {
    for (t in 1:n) {
      if (s != t) {
        for (v in 1:n) {
          if (s != v && t != v) {
            if (D[s, v] + D[v, t] == D[s, t]) {  # Check if v is on the shortest path from s to t
              BC[v] <- BC[v] + (NP[s, v] * NP[v, t]) / NP[s, t]
            }
          }
        }
      }
    }
  }

  return(BC)
}
