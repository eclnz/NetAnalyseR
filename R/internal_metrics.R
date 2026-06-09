# Internal (unexported) metric helpers — no validation, trailing _ suffix.
# These are called by public wrappers after validate_matrix() has already run,
# and by compute_global_metrics / compute_nodal_metrics via the dispatch maps.
#
# C++ symbols used here (registered in package namespace, always available):
#   localClusteringCoefficientWei, localEfficiencyWei,
#   floydWarshallRcpp, dijkstraAllPairs, calculateBetweennessCentrality

# --- Building blocks ---

length_inversion_ <- function(W) {
  diag(W) <- 0
  E <- which(W != 0, arr.ind = TRUE)
  W[E] <- 1 / W[E]
  W
}

# Floyd-Warshall O(n^3) outperforms Dijkstra O(n * E * log n) when the graph
# is dense (most node pairs are directly connected); empirically the crossover
# is around density 0.55.
shortest_distance_ <- function(L) {
  diag(L) <- 0
  if (network_density_(L) > 0.55) floydWarshallRcpp(L) else dijkstraAllPairs(L)
}

network_density_ <- function(W) {
  diag(W) <- 0
  A <- W > 0
  mean(A[lower.tri(A)])
}

# --- Global metrics ---

characteristic_path_length_ <- function(W) {
  diag(W) <- 0
  D <- shortest_distance_(length_inversion_(W))
  D[is.infinite(D)] <- NA
  mean(D[lower.tri(D)])
}

global_clustering_coefficient_wei_ <- function(W) {
  mean(localClusteringCoefficientWei(W))
}

global_efficiency_wei_ <- function(W) {
  diag(W) <- 0
  n <- nrow(W)
  inv_D <- 1 / shortest_distance_(length_inversion_(W))
  sum(inv_D[upper.tri(inv_D)], na.rm = TRUE) / (n * (n - 1) / 2)
}

inter_node_ <- function(W) {
  diag(W) <- 0
  sum(W[lower.tri(W)])
}

intra_node_ <- function(W) {
  sum(diag(W))
}

missing_weights_ <- function(W) {
  inter <- inter_node_(W)
  intra <- intra_node_(W)
  total <- signif(inter + intra, 1)
  total - (inter + intra)
}

normalised_clustering_coefficient_ <- function(W, rand_array) {
  c_obs <- global_clustering_coefficient_wei_(W)
  if (c_obs == 0) return(0)
  rand_c <- mean(apply(rand_array, 3, global_clustering_coefficient_wei_))
  c_obs / rand_c
}

normalised_characteristic_path_length_ <- function(W, rand_array) {
  cpl <- characteristic_path_length_(W)
  if (cpl == 0) return(0)
  rand_cpl <- mean(apply(rand_array, 3, characteristic_path_length_))
  cpl / rand_cpl
}

# --- Nodal metrics ---

node_strength_ <- function(W) {
  diag(W) <- 0
  colSums(W)
}

local_clustering_coefficient_wei_ <- function(W) {
  localClusteringCoefficientWei(W)
}

local_efficiency_wei_ <- function(W) {
  localEfficiencyWei(W)
}

self_connectivity_ <- function(W) {
  diag(W)
}

betweenness_wei_ <- function(W) {
  diag(W) <- 0
  L <- length_inversion_(W)
  D <- shortest_distance_(L)
  NP <- shortest_distance_((L > 0) + 0)
  calculateBetweennessCentrality(D, NP)
}

normalise_inter_node_ <- function(W, target) {
  W * (target / inter_node_(W))
}
