#' @title Normalise a Network Matrix to a Target Total Inter-Node Weight
#'
#' @description Scales all edge weights in a weighted connection matrix so that the
#' total inter-node weight (sum of the lower triangle) equals the specified target.
#'
#' @param mat A square, numeric matrix representing the weighted connection matrix.
#' @param target A single numeric value giving the desired total inter-node weight.
#' @return A matrix of the same dimensions as \code{mat} with all weights rescaled.
#' @export
normalise_inter_node_ <- function(W, target) {
  W * (target / inter_node_(W))
}

normalise_inter_node <- function(mat, target) {
  normalise_inter_node_(validate_matrix(mat), target)
}
