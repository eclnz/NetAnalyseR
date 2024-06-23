#' Allocate subjects to groups based on conditions in a grouping list
#'
#' @title Subject Group Allocation
#' @description Dynamically allocates subjects to groups based on matching conditions specified within a grouping list.
#' This function iterates over each subject identifier, checks against provided patterns in the grouping list, and assigns
#' the corresponding group label if a match is found. It enforces strict input types for the data frame and grouping list
#' and ensures the presence of a 'subject' column in the data frame.
#' @param data_frame A data frame that includes a 'subject' column for participant identifiers. Generally produced with `compute_global_metrics()` or `compute_nodal_metrics()`
#' @param grouping_list A list where each element is named after a group and contains strings or patterns to match against 'subject' identifiers.
#' @return A data frame with an additional 'group' column indicating the group allocation for each subject.
#' @examples
#' data_dir <- system.file("extdata", package = "NetAnalyseR")
#' subjects <- c("A", "B", "C", "D")
#' file_convention <- ".csv"
#' output <- process_matrices(data_dir, subjects, file_convention)
#' global_metrics <- c("characteristic_path_length", "global_efficiency_wei")
#' global_df <- compute_global_metrics(matrices_array = output$matrices,
#'                                     global_metrics = global_metrics,
#'                                     subject_names = output$subjects)
#' grouping_list <- list(
#'   control = c("A", "B"),
#'   treatment = c("C", "D")
#' )
#' allocated_df <- allocate_groups(global_df, grouping_list)
#' @export
#' @importFrom dplyr mutate
#' @importFrom purrr map_chr
#' @importFrom stringr str_detect

allocate_groups <- function(data_frame, grouping_list) {
  # Validate input data_frame format
  if (!is.data.frame(data_frame)) {
    stop("Dataframe is not in dataframe format")
  }

  # Validate input grouping_list format
  if (!is.list(grouping_list)) {
    stop("Grouping list is not in list format")
  }

  # Ensure 'subject' column exists in the data_frame
  if (!"subject" %in% names(data_frame)) {
    stop("Dataframe does not contain a 'subject' column")
  }

  data_frame %>%
    dplyr::mutate(group = purrr::map_chr(subject, function(subject_id) {
      # Initialize found_group as NA
      found_group <- NA_character_
      # Loop through each group in the grouping list
      for (group_name in names(grouping_list)) {
        # Check if subject_id matches any condition in the current group
        if (any(sapply(grouping_list[[group_name]], function(condition) stringr::str_detect(subject_id, condition)))) {
          found_group <- group_name
          break # Break the loop once a match is found
        }
      }
      found_group
    }))
}

#' @title Generate Symmetric Matrix with Random Values
#' @description Creates a symmetric matrix with absolute random normal values.
#' @param dim Integer specifying the dimension of the matrix.
#' @return A symmetric matrix with dimensions `dim` x `dim`.
#' @examples
#' mat <- mat_gen(5)
#' @export
mat_gen <- function(dim) {
  # Generate a matrix with absolute random normal values
  mat <- matrix(abs(rnorm(dim * dim, mean = 0, sd = 100)), nrow = dim, ncol = dim)

  # Make the matrix symmetric
  mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]

  return(mat)
}

#' @title Generate Array of Thresholded Symmetric Matrices
#' @description Creates an array of symmetric matrices with values above a specified threshold.
#' @param dim Integer specifying the dimension of each matrix.
#' @param length Integer specifying the number of matrices in the array.
#' @param threshold Numeric value to threshold matrix entries (default is 0).
#' @return An array of symmetric matrices with dimensions `dim` x `dim` x `length`.
#' @examples
#' mat_array <- mat_array_gen(5, 10, threshold = 50)
#' @export
mat_array_gen <- function(dim, length, threshold = 0) {
  # Initialize an array to hold matrices
  mat_array <- array(dim = c(dim, dim, length))

  # Generate and threshold each matrix
  for (i in 1:length) {
    mat_array[,,i] <- threshold_mat(mat_gen(dim), threshold)
  }

  return(mat_array)
}

#' @title Apply Threshold to Matrix
#' @description Applies a threshold to a matrix, setting elements below the threshold to zero.
#' @param mat A numeric matrix.
#' @param threshold Numeric value to threshold matrix entries.
#' @return A matrix with values below the threshold set to zero.
#' @examples
#' thresholded_mat <- threshold_mat(mat_gen(5), 50)
#' @export
threshold_mat <- function(mat, threshold) {
  mat[abs(mat) < threshold] <- 0
  return(mat)
}

#' @title Update Progress of Matrix Processing
#' @description Displays the progress of processing matrices, including elapsed time, remaining time, and estimated time of completion.
#' @param slice_counter Integer specifying the number of matrices processed so far.
#' @param total_slices Integer specifying the total number of matrices to process.
#' @param start_time POSIXct object indicating the start time of the processing.
#' @param metric_name Character string specifying the name of the metric being processed.
#' @return None. Prints progress information to the console.
#' @examples
#' \dontrun{
#' start_time <- Sys.time()
#' for (i in 1:100) {
#'   Sys.sleep(0.1)
#'   update_progress(i, 100, start_time, "Example Metric")
#' }
#' }
#' @export
update_progress <- function(slice_counter, total_slices, start_time, metric_name, max_metric_length) {
  current_time <- Sys.time()
  elapsed_time <- as.numeric(difftime(current_time, start_time, units = "secs"))
  rate <- slice_counter / elapsed_time
  remaining_slices <- total_slices - slice_counter
  estimated_remaining_time <- remaining_slices / rate
  finish_time <- current_time + estimated_remaining_time
  finish_time_formatted <- format(finish_time, "%Y-%m-%d %H:%M:%S")

  # Format remaining time as minutes and seconds
  remaining_time_formatted <- ifelse(estimated_remaining_time >= 60,
                                     sprintf("%d min %.0f sec", floor(estimated_remaining_time / 60), estimated_remaining_time %% 60),
                                     sprintf("%.2f sec", estimated_remaining_time))
  # Make sure max character length is greater than metric name
  if(max_metric_length<nchar(metric_name)){
    max_metric_length <- nchar(metric_name)
  }
  # Create padding for alignment
  padding <- paste0(rep(" ", max_metric_length - nchar(metric_name)), collapse = "")

  cat(sprintf(
    "\rMetric %s%s: Completed %d of %d slices | Elapsed: %.2f s | Remaining: %s | ETA: %s",
    metric_name, padding, slice_counter, total_slices, elapsed_time, remaining_time_formatted, finish_time_formatted))

  flush.console()
}

