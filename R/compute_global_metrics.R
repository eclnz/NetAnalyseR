#' @title Compute Global Metrics for Brain Network Data
#'
#' @description This function calculates specified global metrics for a given set of brain network matrices
#' for multiple subjects. It validates the requested global metrics against a list of supported metrics, warns
#' about invalid metrics, and computes the specified metrics for each subject's brain network.
#'
#' @param matrices_array An array where each slice represents a connectivity matrix for a subject.
#' @param global_metrics A vector of strings specifying the global metrics to be calculated.
#' @param subject_names An optional vector of subject identifiers. If not provided, subjects will be named
#' sequentially as Subject1, Subject2, etc.
#'
#' @return A data frame containing the calculated global metrics for each subject, with each row corresponding
#' to a subject and columns for the subject identifier and the calculated metric values.
#'
#' @examples
#' W <- matrix(c(0, 2, 1, 4, 2, 0, 3, 5, 1, 3, 0, 6, 4, 5, 6, 0), nrow = 4, byrow = TRUE)
#' matrices_array <- array(W, dim = c(3, 3, 1))
#' global_metrics <- c("characteristic_path_length", "global_efficiency_wei")
#' subject_names <- c("Subject1")
#' global_df <- compute_global_metrics(matrices_array, global_metrics, subject_names)
#'
#' @importFrom dplyr select
#' @importFrom magrittr %>%
#' @export

compute_global_metrics <- function(matrices_array, global_metrics, subject_names = NULL) {
  # Ensure that 'dplyr' and 'magrittr' are available for the select function and %>% operator
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("The 'dplyr' package is required but not installed. Please install it using install.packages('dplyr').")
  }
  if (!requireNamespace("magrittr", quietly = TRUE)) {
    stop("The 'magrittr' package is required for the '%>%' operator but not installed. Please install it using install.packages('magrittr').")
  }

  # Define valid global metrics
  valid_global_metrics <- c("characteristic_path_length", "global_clustering_coefficient_wei", "global_efficiency_wei", "inter_node", "intra_node", "missing_weights", "network_density")

  # Identify any specified metrics that are not valid
  invalid_metrics <- global_metrics[!global_metrics %in% valid_global_metrics]
  if (length(invalid_metrics) > 0) {
    warning("The following global metrics are invalid: \n - ",
            paste0(invalid_metrics, collapse = ", \n - "))
  }

  # Filter for valid metrics specified by the user
  valid_user_metrics <- global_metrics[global_metrics %in% valid_global_metrics]
  if (length(valid_user_metrics) == 0) {
    stop("No valid global metrics were specified. Please use any of the following global metrics:\n - ",
         paste0(valid_global_metrics, collapse = ", \n - "))
  }

  # Warn if both valid and invalid metrics were specified
  if (length(valid_user_metrics) > 0 & length(invalid_metrics) > 0) {
    warning("The following global metrics can be used:\n - ",
            paste0(valid_global_metrics, collapse = ", \n - "))
  }

  # Compute specified valid metrics
  global <- lapply(valid_user_metrics, function(metric_function) {
    apply(matrices_array, MARGIN = 3, FUN = get(metric_function))
  })

  # Create a data frame from the global metrics
  global_df <- data.frame(lapply(global, as.numeric))
  colnames(global_df) <- valid_user_metrics

  # Add subject identifiers to the data frame
  if(is.null(subject_names)) {
    num_subjects <- dim(matrices_array)[3]
    global_df$subject <- paste0("Subject", seq_len(num_subjects))
  } else {
    global_df$subject <- subject_names
  }

  # Reorder columns to place subject identifier first
  global_df <- dplyr::select(global_df, .data$subject, dplyr::everything())

  return(global_df)
}
