
network_deviation <- function(control_array, matrices_array, subject_names = NULL, threshold = 2.5) {
  # Calculate the mean and standard deviation for the control group
  mean_control <- apply(control_array, c(1, 2), mean)
  mean_control_bin_35 <- apply(control_array>0, c(1, 2), mean)>0.35
  sd_control <- apply(control_array, c(1, 2), sd)

  # Function to calculate deviation for each subject
  calculate_deviation <- function(subject_matrix) {
    # Calculate Z-scores
    z_matrix <- abs((subject_matrix - mean_control) / sd_control)
    # Zero diagonal
    diag(z_matrix)<- 0
    # Replace NaNs created by 0/0 and Infs created by x/0 with 0
    z_matrix[is.nan(z_matrix)|is.infinite(z_matrix)] <- 0
    # Identify edges in subject which are present in 35% of controls
    z_matrix <- ((subject_matrix>0)*mean_control_bin_35)*z_matrix
    # Identify edges with deviations beyond the threshold
    edge_outlier <- z_matrix > threshold
    # Sum deviations for each node
    edge_deviations <- colSums(edge_outlier)
    # Network deviation as the sum of subjects with deviations beyond threshold
    sum(edge_deviations > 1)
  }

  # Apply deviation calculation across all subjects
  network_deviations <- apply(matrices_array, MARGIN = 3, FUN = calculate_deviation)

  # Generate subject identifiers
  subject_identifiers <- if(is.null(subject_names)) {
    paste0("subject_", seq_along(network_deviations))
  } else {
    subject_names
  }

  # Return a data frame of subject identifiers and their network deviations
  data.frame(subject = subject_identifiers, network_deviation = network_deviations)
}

#' Compute network deviations for a given dataset and matrices array
#'
#' @title Compute Network Deviations
#' @description This function calculates network deviations for subjects based on the provided control and subject groups within a data frame.
#' It extracts control and subject matrices from a 3D array based on group classification in the data frame and then calculates network deviations
#' using a specified threshold.
#' @param data_frame A data frame containing at least two columns: `group` and `subject` for participant identifiers. The
#' @param matrices_array A 3D array where each slice corresponds to a matrix for a specific participant.
#' @param control_group A character specifying the level of `group` which should be considered the control group.
#' @param threshold A numeric value defining the Z-score threshold for considering deviations as significant.
#' @return A data frame that joins the input data_frame with computed network deviations for each case.
#' @examples
#' # data_frame is a data frame with columns 'group' and 'subject'.
#' # matrices_array is a 3D array where the third dimension corresponds to individuals in the data_frame.
#' # Assuming network_deviation is a predefined function that takes a threshold as an argument.
#' result <- compute_network_deviation(data_frame, matrices_array, 2.5)
#' @references Gugger, J. J., Sinha, N., Huang, Y., Walter, A. E., Lynch, C., Kalyani, P., Smyk, N., Sandsmark, D., Diaz-Arrastia, R., & Davis, K. A. (2023). Structural brain network deviations predict recovery after traumatic brain injury. NeuroImage clinical, 38, 103392-103392. https://doi.org/10.1016/j.nicl.2023.103392
#' @export

compute_network_deviation <- function(data_frame, matrices_array, control_group, threshold = 3) {
  if(!"group" %in% names(data_frame)) {
    stop("The dataframe is missing the 'group' column.", call. = FALSE)
  }
  if(!"subject" %in% names(data_frame)) {
    stop("The dataframe is missing the 'subject' column.", call. = FALSE)
  }
  if(!is.array(matrices_array)){
    stop("The matrices array is not in array format.", call. = FALSE)
  }
  if(!(is.character(control_group) & length(control_group)==1)){
    stop("The control group specified must be a character string with a length of 1.", call. = FALSE)
  }
  if(!any(length(data_frame$group==control_group))){
    stop("No observations of the control group ", control_group, " were found", call. = FALSE)
  }
  # Identify control and case indices based on the group column
  ctrl_indices <- which(data_frame$group == control_group)
  # Extract control and case arrays based on identified indices
  control_array <- matrices_array[, , ctrl_indices]
  # Calculate the network deviation of each network from the control group
  network_deviation_df <- network_deviation(control_array, matrices_array, data_frame$subject, threshold = threshold)
  # Join the original data_frame with the computed network deviation data frame on subject
  merged_df <- dplyr::left_join(data_frame, network_deviation_df, by = "subject")

  return(merged_df)
}

nodal_deviation <- function(control_array, matrices_array, subject_names = NULL, threshold) {
  # Calculate the mean and standard deviation for the control group
  mean_control <- apply(control_array, c(1, 2), mean)
  sd_control <- apply(control_array, c(1, 2), sd)
  mean_control_bin_35 <- apply(control_array>0, c(1, 2), mean)>0.35


  # Function to calculate deviation for each case
  calculate_deviation <- function(case_matrix) {
    # Calculate Z-scores
    z_matrix <- abs((case_matrix - mean_control) / sd_control)
    # Zero diagonal
    diag(z_matrix)<- 0
    # Replace NaNs created by 0/0 and Infs created by x/0 with 0
    z_matrix[is.nan(z_matrix)|is.infinite(z_matrix)] <- 0
    # Identify edges in case which are present in 35% of controls
    z_matrix <- ((case_matrix>0)*mean_control_bin_35)*z_matrix
    # Identify edges with deviations beyond the threshold
    edge_outlier <- z_matrix > threshold
    # Sum deviations for each node
    edge_deviations <- colSums(edge_outlier)
    return(edge_deviations)
  }

  # Apply nodal deviation calculation across all cases
  nodal_deviations <- apply(matrices_array, MARGIN = 3, FUN = calculate_deviation)
  # Prepare data frame for results
  nodal_deviations_df <- as.data.frame(t(nodal_deviations))
  # Assign column names based on node index
  colnames(nodal_deviations_df) <- paste0("node", seq_len(ncol(nodal_deviations_df)))

  # Generate subject identifiers
  subject_identifiers <- if(is.null(subject_names)) {
    paste0("subject_", seq_along(nodal_deviations_df$node1))
  } else {
    subject_names
  }

  # Add case identifiers to the data frame
  nodal_deviations_df$subject <- subject_identifiers

  # Pivot the data frame to a long format for easier analysis
  nodal_deviations_df <- tidyr::pivot_longer(nodal_deviations_df, cols = starts_with("node"),
                                             names_to = "node", values_to = "nodal_deviation") %>%
    dplyr::mutate(node = as.integer(gsub("node", "", node)))

  return(nodal_deviations_df)
}

#' Compute nodal network deviations and merge with original data frame
#'
#' @title Compute Nodal Network Deviations
#' @description This function computes nodal network deviations for all subjects in a dataset,
#' combining control and case data based on a specified threshold. It leverages the `nodal_deviation`
#' function to calculate deviations and merges these results back into the original data frame.
#' @param data_frame A data frame where each row corresponds to a single subject. Typically produced using `compute_global_metrics()`.
#' @param matrices_array A 3D array where each slice corresponds to a connectivity matrix for a subject.
#' @param threshold A numeric value defining the Z-score threshold for considering deviations significant.
#' @return A modified version of the original data frame that includes nodal deviation calculations for each subject.
#' @examples
#' # Assuming a prepared data_frame with 'subject' and 'group' columns, and matrices_array where each slice [,,i] corresponds to a subject's connectivity matrix.
#' # A threshold value is also specified for deviation calculation.
#' merged_results <- compute_nodal_network_deviation(data_frame, matrices_array, threshold = 2.5)
#' @export

compute_nodal_network_deviation <- function(data_frame, matrices_array, control_group, threshold= 3) {
  if(!"group" %in% names(data_frame)) {
    stop("The dataframe is missing the 'group' column.", call. = FALSE)
  }
  if(!"subject" %in% names(data_frame)) {
    stop("The dataframe is missing the 'subject' column.", call. = FALSE)
  }
  if(!is.array(matrices_array)){
    stop("The matrices array is not in array format.", call. = FALSE)
  }
  if(!(is.character(control_group) & length(control_group)==1)){
    stop("The control group specified must be a character string with a length of 1.", call. = FALSE)
  }
  if(!any(length(data_frame$group==control_group))){
    stop("No observations of the control group ", control_group, " were found", call. = FALSE)
  }
  # Dataframe of controls not in long format (one row per subject the same as the matrices array)
  short_df <- data_frame %>% distinct(subject, group, subject)
  # Identify control and case indices based on the group column
  ctrl_indices <- which(short_df$group == control_group)
  # Extract control and case arrays based on identified indices
  control_array <- matrices_array[, , ctrl_indices]
  # Calculate nodal deviations using the control group and threshold
  nodal_network_deviation_df <- nodal_deviation(control_array = matrices_array[, , ctrl_indices],
                                                matrices_array = matrices_array,
                                                subject_names = short_df$subject,
                                                threshold = threshold)
  # Merge the calculated nodal network deviations with the original data frame
  merged_df <- dplyr::left_join(data_frame, nodal_network_deviation_df, by = c("subject", "node"))

  return(merged_df)
}


