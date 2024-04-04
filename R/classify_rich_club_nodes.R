#' @title Classify Rich Club Nodes
#'
#' @description This function identifies the top 10 rich club nodes in a network
#' based on their edge strengths, using edges produced by the function `process_matrices()`
#' and classified using the function `classify_connections()`. It specifically
#' focusing on edges between distinct nodes marked as "Network-Connected".
#'
#' @param data_frame A data frame containing the network data with columns:
#' "subject", "row", "col", "edge_strength", "self_con".
#' @return A vector of the top 10 rich club nodes identified within the network.
#' @examples
#' data_frame <- data.frame(subject = c(1, 1, 2, 2),
#'                          row = c("A", "B", "A", "C"),
#'                          col = c("B", "A", "C", "A"),
#'                          edge_strength = c(2, 2, 3, 3),
#'                          self_connectivity = c("Network-Connected",
#'                          "Network-Connected", "Network-Connected",
#'                          "Network-Connected"))
#' classify_rich_club_nodes(data_frame)
#'
#' @export
#'
#' @importFrom dplyr filter group_by summarize arrange slice count ungroup
#' @importFrom tidyr gather
#'
classify_rich_club_nodes <- function(data_frame) {

  # Validate the presence of required columns in the data frame
  required_columns <- c("subject", "row", "col", "edge_strength", "self_connectivity")
  missing_columns <- setdiff(required_columns, names(data_frame))
  if (length(missing_columns) > 0) {
    stop("The data frame is missing the following required columns: ", paste(missing_columns, collapse = ", "))
  }

  # Calculate rich club nodes focusing on network-connected entries
  rich_club_nodes <- data_frame %>%
    dplyr::filter(self_connectivity == "Network-Connected") %>%
    tidyr::pivot_longer(cols = c(row, col), names_to = "node_type", values_to = "node") %>%
    dplyr::group_by(subject, node) %>%
    dplyr::summarize(strength = sum(edge_strength), .groups = "drop") %>%
    dplyr::group_by(subject) %>%
    dplyr::arrange(dplyr::desc(strength)) %>%
    dplyr::slice(1:10) %>%
    dplyr::ungroup() %>%
    dplyr::count(node) %>%
    dplyr::slice_max(order_by = n, n = 10) %>%
    dplyr::ungroup()

  # Return the top 10 rich club nodes
  return(rich_club_nodes$node)
}
