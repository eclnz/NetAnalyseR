#' @title Classify Connections Based on Cortical and Subcortical Allocation
#'
#' @description This function classifies the connections in the edge data frame
#' produced by the function `process_matrices()` according to whether they are
#' cortical-cortical, subcortical-cortical, subcortical-subcortical, or unallocated.
#' It requires the data frame to contain specific columns ('row' and 'col') and
#' uses two vectors indicating cortical and subcortical allocations to classify the
#' connections.
#'
#' @param data_frame A data frame containing the connections to be classified.
#' It must include 'row' and 'col' columns representing the connections between nodes.
#' @param cortical_allocation A vector of nodes allocated as cortical.
#' @param subcortical_allocation A vector of nodes allocated as subcortical.
#' @return A modified version of the input data frame with an additional column named
#' 'connection_type' indicating the classification of each connection (Cortical-Cortical,
#' Subcortical-Cortical, Subcortical-Subcortical, or Unallocated).
#' @examples
#' cortical_nodes <- c(1, 2, 3)
#' subcortical_nodes <- c(4, 5)
#' connections <- data.frame(row = c(1, 2, 4, 5, 3), col = c(2, 3, 5, 1, 4))
#' classified_connections <- classify_connections(connections, cortical_nodes, subcortical_nodes)
#'
#' @export
#'
#' @importFrom dplyr mutate case_when
#'
classify_connections <- function(data_frame, cortical_allocation, subcortical_allocation) {

  # Validate presence of required columns in the data frame
  if (!("row" %in% names(data_frame)) || !("col" %in% names(data_frame))) {
    stop("The data_frame must contain 'row' and 'col' columns.")
  }

  # Identify and warn about unallocated nodes using base R functions
  all_nodes <- unique(c(data_frame$row, data_frame$col))
  unallocated_nodes <- setdiff(all_nodes, c(cortical_allocation, subcortical_allocation))
  if (length(unallocated_nodes) > 0) {
    warning("The following nodes are not allocated: ", paste(unallocated_nodes, collapse = ", "))
  }

  # Classify connections using dplyr for more readable and concise syntax
  data_frame <- dplyr::mutate(data_frame, connection_type = dplyr::case_when(
    (row %in% cortical_allocation & col %in% cortical_allocation) ~ "Cortical-Cortical",
    ((row %in% cortical_allocation & col %in% subcortical_allocation) |
       (row %in% subcortical_allocation & col %in% cortical_allocation)) ~ "Cortical-Subcortical",
    (row %in% subcortical_allocation & col %in% subcortical_allocation) ~ "Subcortical-Subcortical",
    TRUE ~ "Unallocated"
  ))

  # Return the updated data frame with classified connections
  return(data_frame)
}
