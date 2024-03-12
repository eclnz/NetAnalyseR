#' @title Compute Nodal Metrics for Brain Network Data
#'
#' @description This function calculates specified nodal metrics for a given set of brain network matrices
#' for multiple subjects. It computes the specified metrics for each node in each subject's brain network.
#'
#' @param matrices_array An array where each slice represents a connectivity matrix for a subject.
#' @param nodal_metrics A vector of strings specifying the nodal metrics to be calculated.
#' @param subject_names An optional vector of subject identifiers. If not provided, subjects will be named
#' sequentially as Subject1, Subject2, etc.
#'
#' @return A long-format data frame containing the calculated nodal metrics for each node and subject.
#' Each row corresponds to a node-subject combination with columns for node number, subject identifier, and the
#' calculated metric values.
#'
#' @examples
#' W <- matrix(c(0, 2, 1, 4, 2, 0, 3, 5, 1, 3, 0, 6, 4, 5, 6, 0), nrow = 4, byrow = TRUE)
#' matrices_array <- array(W, dim = c(3, 3, 1))
#' nodal_metrics <- c("local_clustering_coefficient_wei", "local_efficiency_wei")
#' subject_names <- c("Subject1")
#' metrics_df <- compute_nodal_metrics(matrices_array, nodal_metrics, subject_names)
#'
#' @export
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr arrange
#' @importFrom magrittr %>%
#'
compute_nodal_metrics <- function(matrices_array, nodal_metrics, subject_names = NULL) {

  # Define valid nodal metrics
  valid_nodal_metrics <- c("local_clustering_coefficient_wei", "local_efficiency_wei", "node_strength", "self_connectivity")

  # Identify any specified metrics that are not valid
  invalid_metrics <- nodal_metrics[!nodal_metrics %in% valid_nodal_metrics]
  if (length(invalid_metrics) > 0) {
    warning("The following nodal metrics are invalid: \n - ",
            paste0(invalid_metrics, collapse = ", \n - "))
  }

  # Filter for valid metrics specified by the user
  valid_user_metrics <- nodal_metrics[nodal_metrics %in% valid_nodal_metrics]
  if (length(valid_user_metrics) == 0) {
    stop("No valid nodal metrics were specified. Please use any of the following nodal metrics:\n - ",
         paste0(valid_nodal_metrics, collapse = ", \n - "))
  }

  # Warn if both valid and invalid metrics were specified
  if (length(valid_user_metrics) > 0 & length(invalid_metrics) > 0) {
    warning("The following nodal metrics can be used:\n - ",
            paste0(valid_nodal_metrics, collapse = ", \n - "))
  }

  # Compute specified valid metrics
  nodal <- lapply(valid_user_metrics, function(metric_function) {
    apply(matrices_array, MARGIN = 3, FUN = get(metric_function))
  })

  # Create data frames for each metric and reshape to long format
  metric_dfs <- lapply(seq_along(nodal), function(i) {
    metric_df <- as.data.frame(nodal[[i]])
    metric_name <- valid_user_metrics[i]
    if (is.null(subject_names)) {
      names(metric_df) <- paste0(metric_name, "-Subject", seq_len(ncol(metric_df)))
    } else {
      names(metric_df) <- paste0(metric_name, "-", subject_names)
    }
    metric_df$node <- seq_len(nrow(metric_df))

    tidyr::pivot_longer(metric_df, cols = dplyr::starts_with(metric_name), names_to = "subject", values_to = metric_name,
                        names_prefix = paste0(metric_name, "-")) %>%
      dplyr::arrange(.data$subject)
  })

  # Combine metric data frames
  combined_df <- Reduce(function(x, y) merge(x, y, by = c("node", "subject"), all = TRUE), metric_dfs) %>%
    dplyr::arrange(.data$subject, .data$node)

  return(combined_df)
}
