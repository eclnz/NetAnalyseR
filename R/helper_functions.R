#' Allocate subjects to groups based on conditions in a grouping list
#'
#' @title Subject Group Allocation
#' @description Dynamically allocates subjects to groups based on matching conditions specified within a grouping list.
#' This function iterates over each subject identifier, checks against provided patterns in the grouping list, and assigns
#' the corresponding group label if a match is found. It enforces strict input types for the data frame and grouping list
#' and ensures the presence of a 'subject' column in the data frame.
#' @param data_frame A data frame that includes a 'subject' column for participant identifiers.
#' @param grouping_list A list where each element is named after a group and contains strings or patterns to match against 'subject' identifiers.
#' @return A data frame with an additional 'group' column indicating the group allocation for each subject.
#' @examples
#' # Assuming a data_frame with 'subject' column and a predefined grouping_list with patterns for each group:
#' grouping_list <- list(
#'   control = c("^ctrl", "control$"),
#'   treatment = c("treat", "experiment")
#' )
#' allocated_df <- allocate_groups(data_frame, grouping_list)
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
