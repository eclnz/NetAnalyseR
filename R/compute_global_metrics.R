#' @title Compute Global Metrics for Brain Network Data
#'
#' @description This function calculates specified global metrics for a given set of brain network matrices
#' for multiple subjects. It validates the requested global metrics against a list of supported metrics, warns
#' about invalid metrics, and computes the specified metrics for each subject's brain network.
#'
#' @param matrices_array An array where each slice represents a connectivity matrix for a subject.
#' @param global_metrics A vector of strings specifying the global metrics to be calculated.
#' @param density An optional float specifying the density all networks should be pruned to.
#' @param target An optional float specifying the density all networks should be pruned to.
#' @param subject_names An optional vector of subject identifiers. If not provided, subjects will be named
#' sequentially as Subject1, Subject2, etc.
#'
#' @return A data frame containing the calculated global metrics for each subject, with each row corresponding
#' to a subject and columns for the subject identifier and the calculated metric values.
#'
#' @examples
#' W <- matrix(c(0, 2, 1, 0, 2, 0, 3, 5, 1, 3, 0, 6, 0, 5, 6, 0), nrow = 4, byrow = TRUE)
#' valid_global_metrics <- c("characteristic_path_length",
#'                           "global_clustering_coefficient_wei",
#'                           "global_efficiency_wei",
#'                           "inter_node",
#'                           "intra_node",
#'                           "network_density",
#'                           "normalised_clustering_coefficient",
#'                           "normalised_characteristic_path_length",
#'                           "small_worldness")
#' subject_names <- c("Subject1")
#' global_df <- compute_global_metrics(W, valid_global_metrics, subject_names)
#' @importFrom dplyr select
#' @importFrom abind abind
#' @export

compute_global_metrics <- function(matrices_array, global_metrics, density_val = NULL, target = NULL, subject_names = NULL) {
  # If the user submits a single matrix the dimensions of the matrix as an array must be set.
  if(is.na(dim(matrices_array)[3])){
    dim(matrices_array)[3] <- 1
  }

  if(!is.array(matrices_array)){
    stop("Matrices array is not in array format")
  }
  # If subject names is specified, stop if it is not character.
  if(!is.null(subject_names)){
    if(!is.character(subject_names)){
      stop("Subject names is not in character format")
    }
  }
  if(!is.character(global_metrics)){
    stop("Global metrics is not in character format")
  }
  # Define valid global metrics
  valid_global_metrics <- c("characteristic_path_length", "global_clustering_coefficient_wei", "global_efficiency_wei", "inter_node", "intra_node", "missing_weights", "network_density", "normalised_clustering_coefficient", "normalised_characteristic_path_length", "small_worldness")
  valid_user_metrics <- global_metrics[global_metrics %in% valid_global_metrics]
  # Identify any specified metrics that are not valid
  invalid_metrics <- global_metrics[!global_metrics %in% valid_global_metrics]
  random_metrics <- c("normalised_clustering_coefficient", "normalised_characteristic_path_length", "small_worldness")
  user_non_random_metrics <- valid_user_metrics[!valid_user_metrics %in% random_metrics]
  user_random_metrics <-valid_user_metrics[valid_user_metrics %in% random_metrics]
  if (length(invalid_metrics) > 0) {
    warning("The following global metrics are invalid: \n - ",
            paste0(invalid_metrics, collapse = ", \n - "))
  }
  if (length(valid_user_metrics) == 0) {
    stop("\nNo valid global metrics were specified. Please use any of the following global metrics:\n - ",
         paste0(valid_global_metrics, collapse = ", \n - "))
  }
  # Warn if both valid and invalid metrics were specified
  if (length(valid_user_metrics) > 0 & length(invalid_metrics) > 0) {
    warning("The following global metrics can be used:\n - ",
            paste0(valid_global_metrics, collapse = ", \n - "))
  }

  # Warn if both valid and invalid metrics were specified
  if (any(apply(matrices_array, MARGIN= 3, FUN = network_density) == 1)) {
    warning("Network density is 1 in some networks. Random Networks cannot be created and random metrics will not be computed")
    user_random_metrics <- character(0)

  }

  # # Calculate Processing Time
  # If the user has a small number of smaller matrices to compute, then set benchmark to 1.
  if(length(matrices_array)<125000){
    user_benchmark <- 1
  } else{
    # If the user does not want to calculate intensive metrics then set user benchmark to 1.
    if(!any(c("characteristic_path_length", "global_clustering_coefficient_wei", "global_efficiency_wei") %in% valid_user_metrics)){
      user_benchmark <- 1
    } else{
      user_benchmark <- benchmark_performance()
    }
  }

  # Normalise arrays by density
  if(!is.null(density_val)){
     if(density_val>=1 | density_val<= 0){
       warning('density must be between 0 and 1')
     }
     else{
       for(i in 1:dim(matrices_array)[3]){
         matrices_array[,,i] <- threshold_density(matrices_array[,,i], density_val)
       }
     }
  }

  # Normalize arrays by target
  if(!is.null(target)){
    for(i in 1:dim(matrices_array)[3]){
      matrices_array[,,i] <- normalise_inter_node(matrices_array[,,i], target)
    }
  }


  estimate_total_duration(matrices_array,user_benchmark,valid_user_metrics)

  if(length(user_non_random_metrics)>0){
    global <- lapply(user_non_random_metrics, function(metric_function) {
      apply(matrices_array, MARGIN = 3, FUN = get(metric_function))
      })
  }

  # Create a data frame from the global metrics
  global_df <- data.frame(lapply(global, as.numeric))
  colnames(global_df) <- user_non_random_metrics

  if(length(user_random_metrics)>0){

    rand_array <- apply(matrices_array, MARGIN = 3, generateRewiredMatrices) %>%
      lapply(function(sublist) {
        abind(sublist, along = 3)
      })

    combined_list <- lapply(seq_len(dim(matrices_array)[3]), function(i) {
      list(
        original_matrix = matrices_array[,,i],  # Extract the i-th matrix from the 3D array
        associated_array = rand_array[[i]]  # Get the corresponding processed array from the list
      )
    })

    if("normalised_clustering_coefficient" %in% user_random_metrics || "small_worldness" %in% user_random_metrics) {
      norm_clust <- lapply(combined_list, function(item) {
        normalised_clustering_coefficient(list(item$original_matrix, item$associated_array))
      }) %>%
        unlist()
      global_df$normalised_clustering_coefficient <- norm_clust
    }

    if("normalised_characteristic_path_length" %in% user_random_metrics || "small_worldness" %in% user_random_metrics) {
      norm_cpl <- lapply(combined_list, function(item) {
        normalised_characteristic_path_length(list(item$original_matrix, item$associated_array))
      }) %>%
        unlist()
      global_df$normalised_characteristic_path_length <- norm_cpl
    }

    if( "small_worldness" %in% user_random_metrics){
      global_df$small_world <- norm_clust/norm_cpl
    }
  }

  # Add subject identifiers to the data frame
  if(is.null(subject_names)) {
    num_subjects <- dim(matrices_array)[3]
    global_df$subject <- paste0("Subject", seq_len(num_subjects))
  } else {
    global_df$subject <- subject_names
  }

  # Reorder columns to place subject identifier first
  global_df <- dplyr::select(global_df, subject, dplyr::everything())

  return(global_df)
}
