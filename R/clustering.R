#' @title Weighted Semi-Synchronous Label Propagation
#' @description Perform weighted semi-synchronous label propagation on a given
#' weighted adjacency matrix. Automatically determines the number of clusters to
#' generate.
#' @param weighted_adj_matrix A square matrix representing the weighted adjacency matrix of the graph.
#' @return A vector of labels indicating the community assignment for each node.
#' @examples
#' # Create a sample weighted adjacency matrix
#' W <- matrix(c(0, 2, 1, 0, 2, 0, 3, 5, 1, 3, 0, 6, 0, 5, 6, 0), nrow = 4, byrow = TRUE)
#' # Perform weighted semi-synchronous label propagation
#' labels <- label_prop_clust(W)
#' print(labels)
#' @references Cordasco, G, and L Gargano. “Community Detection via Semi-Synchronous Label Propagation Algorithms.” 2010 IEEE International Workshop on: Business Applications of Social Network Analysis (BASNA). IEEE, 2010. 1–8. Web.
#' @export
label_prop_clust <- function(weighted_adj_matrix) {

  # Number of nodes
  num_nodes <- nrow(weighted_adj_matrix)

  # Initialize labels
  node_labels <- 1:num_nodes  # Each node has a unique label

  # Function to update labels for one semi-synchronous phase
  update_labels <- function() {
    updated_labels <- node_labels

    # Randomly determine update order to simulate semi-synchronicity
    update_order <- sample(num_nodes)

    for (node in update_order) {

      # Neighbors of the node and their corresponding weights
      neighbors <- which(weighted_adj_matrix[node, ] > 0)

      if (length(neighbors) > 0) {

        # Get labels and weights of neighbors
        neighbor_labels <- node_labels[neighbors]
        neighbor_weights <- weighted_adj_matrix[node, neighbors]

        # Calculate the weighted frequency of each label
        weighted_label_counts <- tapply(neighbor_weights, neighbor_labels, sum)

        # Find the label with the maximum weighted count
        most_common_label <- names(which.max(weighted_label_counts))
        updated_labels[node] <- as.integer(most_common_label)
      }
    }
    return(updated_labels)
  }

  iteration <- 0
  max_iterations <- 100  # Limit iterations to prevent infinite loops

  while(TRUE) {
    iteration <- iteration + 1
    new_labels <- update_labels()

    if (identical(new_labels, node_labels) || iteration >= max_iterations) {
      break
    }

    node_labels <- new_labels
  }

  return(node_labels)
}

#' @title K-Means Distance-Based Clustering
#' @description Perform K-means clustering based on the shortest distance matrix derived from the input matrix.
#' @param mat An adjacency matrix from which the shortest distance matrix will be computed.
#' @param k An integer specifying the number of clusters. Default is set to one-eighth of the number of rows in the matrix.
#' @return A vector of cluster labels for each row in the input matrix.
#' @examples
#' # Create a sample matrix
#' W <- matrix(c(0, 2, 1, 0, 2, 0, 3, 5, 1, 3, 0, 6, 0, 5, 6, 0), nrow = 4, byrow = TRUE)
#' # Perform K-means clustering
#' cluster_labels <- K_dis_clust(W)
#' print(cluster_labels)
#' @export
#' @importFrom stats kmeans
K_dis_clust <- function(mat, k = floor(dim(mat)[1]/2)) {
  if(k<2){
    return(rep(1,dim(mat)[1]))
  }
  # Compute the length inversion matrix from the input matrix
  len <- length_inversion(mat)
  # Compute the shortest distance matrix from the length inversion matrix
  dis <- shortest_distance(len)
  # Apply K-means clustering on the shortest distance matrix
  kmeans_result <- kmeans(dis, centers = k)
  # Return the cluster labels
  return(kmeans_result$cluster)
}

