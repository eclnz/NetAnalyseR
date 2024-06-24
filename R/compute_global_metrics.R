#' @title Compute Global Metrics for Brain Network Data
#'
#' @description This function calculates specified global metrics for a given set of brain network matrices
#' for multiple subjects. It validates the requested global metrics against a list of supported metrics, warns
#' about invalid metrics, and computes the specified metrics for each subject's brain network.
#'
#' @param matrices_array An array where each slice represents a connectivity matrix for a subject.
#' @param global_metrics A vector of strings specifying the global metrics to be calculated.
#' @param density_val An optional float specifying the density all networks should be pruned to.
#' @param target An optional float specifying the total network strength all networks should be normalized to
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
#' global_df <- compute_global_metrics(W, valid_global_metrics, NULL, NULL, subject_names)
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

  # Normalise arrays by density
  if (!is.null(density_val)) {
     if(density_val>=1 || density_val<= 0){
       warning('density must be between 0 and 1')
     }
     else{
      for (i in 1:dim(matrices_array)[3]) {
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

  # Validate matrices within array
  apply(matrices_array, MARGIN = 3, FUN = validate_matrix)

  # Initialize a list to store results for each metric
  global <- list()

  # Max character number
  max_nchar <- nchar("normalised_characteristic_path_length")

  # Process user non-random metrics
  if (length(user_non_random_metrics) > 0) {
    for (metric_idx in seq_along(user_non_random_metrics)) {
      metric_function <- user_non_random_metrics[[metric_idx]]
      metric_name <- as.character(metric_function)

      # Initialize a vector to store the results for the current metric
      metric_results <- vector("list", dim(matrices_array)[3])

      # Record the start time
      start_time <- Sys.time()

      # Loop over each slice of the matrix array
      for (slice_idx in 1:dim(matrices_array)[3]) {
        # av the metric function to the current slice
        metric_results[[slice_idx]] <- get(metric_function)(matrices_array[,,slice_idx], FALSE)

        # Update progress
        update_progress(slice_idx, dim(matrices_array)[3],start_time, metric_name, max_nchar)
      }

      # Store the results for the current metric
      global[[metric_name]] <- metric_results
      cat("\n")
    }
  }
  # Create a data frame from the global metrics
  global_df <- data.frame(lapply(global, as.numeric))
  colnames(global_df) <- user_non_random_metrics

  # Process user random metrics
  if (length(user_random_metrics) > 0) {
    # Initialize the rand_array
    rand_array <- vector("list", dim(matrices_array)[3])

    # Record the start time for rand_array generation
    start_time_rand <- Sys.time()

    # Generate rewired matrices with progress and ETA
    for (slice_idx in 1:dim(matrices_array)[3]) {
      rand_array[[slice_idx]] <- generateRewiredMatrices(matrices_array[,,slice_idx])

      update_progress(slice_idx, dim(matrices_array)[3],start_time_rand, "Rewired Networks Generation", max_nchar)
    }
    cat("\n")
    # Combine rand_array into a single 3D array
    rand_array <- lapply(rand_array, function(sublist) {
      abind(sublist, along = 3)
    })

    combined_list <- lapply(seq_len(dim(matrices_array)[3]), function(i) {
      list(
        original_matrix = matrices_array[,,i],  # Extract the i-th matrix from the 3D array
        associated_array = rand_array[[i]]  # Get the corresponding processed array from the list
      )
    })

    # Start time for random metrics
    start_time_metrics <- Sys.time()

    # Initialize a counter for slices processed
    slice_counter <- 0

    if ("normalised_clustering_coefficient" %in% user_random_metrics || "small_worldness" %in% user_random_metrics) {
      norm_clust <- lapply(combined_list, function(item) {
        slice_counter <<- slice_counter + 1
        # Apply the metric function to the current slice
        result <- normalised_clustering_coefficient(list(item$original_matrix, item$associated_array), validate = FALSE)

        # Update progress for normalised_clustering_coefficient
        update_progress(slice_counter, dim(matrices_array)[3], start_time_metrics, "normalised_clustering_coefficient", max_nchar)
        return(result)
      }) %>% unlist()
      global_df$normalised_clustering_coefficient <- norm_clust
    }
    cat("\n")

    if ("normalised_characteristic_path_length" %in% user_random_metrics || "small_worldness" %in% user_random_metrics) {
      slice_counter <- 0  # Reset counter for next metric
      norm_cpl <- lapply(combined_list, function(item) {
        slice_counter <<- slice_counter + 1
        # Apply the metric function to the current slice
        result <- normalised_characteristic_path_length(list(item$original_matrix, item$associated_array), validate = FALSE)

        # Update progress for normalised_characteristic_path_length
        update_progress(slice_counter, dim(matrices_array)[3], start_time_metrics, "normalised_characteristic_path_length", max_nchar)
        return(result)
      }) %>% unlist()
      cat("\n")
      global_df$normalised_characteristic_path_length <- norm_cpl

    }

    if ("small_worldness" %in% user_random_metrics) {
      global_df$small_world <- norm_clust / norm_cpl
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


