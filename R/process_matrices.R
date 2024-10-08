#' @title Process and analyse matrices from files within a directory
#'
#' @description This function processes matrices stored in CSV files within a specified directory,
#'   filtering by subject presence and adhering to a specified file naming convention. It validates
#'   file existence, dimensions consistency, and generates an edge dataframe and an array of matrices.
#'
#' @param directory The path to the directory containing the matrix files.
#' @param file_convention The naming convention for the matrix files related to each subject.
#' @param subjects_specified A vector of subject identifiers to be processed.
#' @param output_edge A boolean whether to calculate and output a dataframe containing the edge strengths of all connections in all networks.
#'
#' @return A list containing three elements: 'edge_df' which is a dataframe of edge strengths between lower
#'   triangular matrix elements across all subjects, 'matrices' which is a 3D array of matrices,
#'   where the matrices of different subjects are stacked in the Z axis, and 'subjects' which is a vector
#'   containing the names of subjects who have matrices included in the dataset.
#'
#' @examples
#' data_dir <- system.file("extdata", package = "NetAnalyseR")
#' file_convention <- ".csv"
#' process_matrices(data_dir, file_convention)
#'
#' @importFrom dplyr mutate
#' @importFrom magrittr %>%
#' @importFrom abind abind
#' @importFrom utils read.csv
#'
#' @export
#'
process_matrices <- function(directory, file_convention, subjects_specified=NULL, output_edge = FALSE) {
  # Check if the specified directory exists
  if (!dir.exists(directory)) {
    stop("Directory does not exist: ", directory)
  }
  if(is.null(subjects_specified)){
    subjects_specified <- unique(list.files(directory, recursive=FALSE))
    subjects_specified <- subjects_specified[grepl(file_convention, subjects_specified)]
  } else{
    subjects_specified <- paste0(subjects_specified,file_convention)
  }
  if (!is.character(subjects_specified) & !is.numeric(subjects_specified)){
      stop("If specifying subjects, it must be supplied as a character or number: ", subjects_specified)
  }

  if (is.character(file_convention)!=TRUE){
    stop("File convention must be supplied as a character: ", file_convention)
  }
  if (length(file_convention)!=1){
    stop("File convention must have a length of 1: ", file_convention)
  }

  # TODO: When the file convention doesn't match anything it needs to warn the user.

  # Initialize lists to store results
  subjects_present <- c()
  edge_df_list <- list()
  matrices_list <- list()
  expected_dims <- NULL

  # Iterate over subjects to process their matrices
  for (subj in subjects_specified) {
    # Construct file path from directory, subject ID, and file convention
    file_path <- file.path(directory, subj)
    subj <- gsub(file_convention, '', subj)

    if (!file.exists(file_path)) {
      warning("Missing file for subject: ", subj)
      next
    }
    subjects_present <- c(subjects_present,subj)

    cat("Processing subject: ", subj, "\n")

    # Read the matrix from the CSV file
    matrix <- as.matrix(read.csv(file_path, header = FALSE))
    # Set or check dimensions against the first matrix processed
    if (is.null(expected_dims)) {
      expected_dims <- dim(matrix)
    } else if (!all(dim(matrix) == expected_dims)) {
      stop("Dimension mismatch for subject: ", subj)
    }

    # Store the matrix in the list
    matrices_list[[subj]] <- matrix

    # Identify lower triangle indices
    lower_triangle_indices <- which(lower.tri(matrix, diag = TRUE), arr.ind = TRUE)

    # Create a dataframe for edges in the lower triangle
    subject_df <- data.frame(
      subject = subj,
      row = lower_triangle_indices[, 1],
      col = lower_triangle_indices[, 2],
      edge_strength = matrix[lower_triangle_indices]
    )

    # Generate an edge identifier
    subject_df$edge <- interaction(subject_df$row, subject_df$col)

    # Store the dataframe in the list
    edge_df_list[[subj]] <- subject_df
  }

  # Ensure at least one file was processed
  if (length(subjects_present) == 0) {
    stop("No files were found for the subjects specified")
  }

  # Combine individual matrices into a 3D array
  matrix_df <- abind(matrices_list, along = 3)

  # Validate matrices within array
  apply(matrix_df, MARGIN = 3, FUN = validate_matrix)

  # Combine individual subject edge dataframes into one
  if (output_edge==TRUE){
    edge_df <- do.call(rbind, edge_df_list) %>%
      dplyr::mutate(self_connectivity = ifelse(row == col, "Self-Connected", "Network-Connected"))
    return(list(matrices = matrix_df, subjects = subjects_present, edge_df = edge_df))
    }

  return(list(matrices = matrix_df, subjects = subjects_present))
}


