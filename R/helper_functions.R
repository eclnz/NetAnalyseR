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
#' global_df <- compute_global_metrics(output$matrices, global_metrics, output$subjects)
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


#' @title Generate a Warning Once
#' @description This function generates a specified warning message only once during the R session.
#' @param warning_var_name A string representing the name of the variable used to track whether the warning has been issued.
#' @param warning_message A string containing the warning message to be displayed.
#' @return No return value; the function solely issues a warning and updates the global environment.
#' @keywords internal
generate_warning <- function(warning_var_name, warning_message) {
  if (!get(warning_var_name, envir = .GlobalEnv)) {
    warning(warning_message)
    assign(warning_var_name, TRUE, envir = .GlobalEnv)
  }
}
