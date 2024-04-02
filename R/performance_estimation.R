#' @keywords internal
#' @importFrom stats rnorm
mat_gen <- function(dim){
  mat <- matrix(abs(rnorm(dim*dim, 0, 100)),nrow = dim, ncol=dim,)
  mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]
  return(mat)
}

#' @keywords internal
mat_threshold <- function(mat, threshold){
  return(as.numeric(mat>threshold) * mat)
}

#' @keywords internal
mat_array_gen <- function(dim,length, threshold = 0){
  mat_array <- array(dim = c(dim,dim,length))
  for(i in 1:length){
    mat_array[,,i] <- mat_threshold(mat_gen(dim),threshold)
  }
  return(mat_array)
}

#' @title Quantifies user thread speed
#' @description This function benchmarks the performance of a network analysis function
#' by comparing observed execution times with estimated times derived from a predictive model.
#' It generates networks with specified size, sparsity, and repeats the process to average results.
#' @return Proportion difference between the mean observed execution time and the mean estimated execution time.
#' @keywords internal
benchmark_performance <- function(){
  size = 100
  repeats = 10
  sparsity = 0.2
  number = 10

  # model_path <- paste0(getwd(),"/inst/extdata/global_efficiency_execution_time_model.rds")
  model_path <- system.file("extdata", "global_efficiency_execution_time_model.rds", package = "NetAnalyseR")
  model <- readRDS(model_path)

  observed <- c()
  estimated <- c()
  for(i in 1:repeats){
    mat_array <- array(dim = c(size,size,number))
    for (j in 1:number) {
      mat_array[,,j] <- mat_threshold(mat_gen(size),100)
    }
    observed[i] <- time_calc(array_process, n = number, global_efficiency_wei, mat_array)

    predict_df <- data.frame(
      Size = apply(mat_array, MARGIN = if(is.matrix(mat_array)) 1 else 3, FUN = nrow),
      MeanDensity = apply(mat_array, MARGIN = if(is.matrix(mat_array)) 1 else 3, FUN = network_density))

    estimated[i] <- (predict(model, predict_df)^3)[1]
  }
  # Returns proportion difference between execution time of machine and expected execution time on M2 Macbook.
  proportion_difference <- (mean(observed)-mean(estimated))/ mean(observed)
  return(proportion_difference)
}

#' @keywords internal
time_calc <- function(function_, n, ...) {
  elapsed <- system.time({function_(...)})[3]
  elapsed_average <- elapsed/n
  return(elapsed_average)
}

#' @keywords internal
array_process <- function(functions, mat_array){
  lapply(list(functions), function(metric_function) {
    apply(mat_array, MARGIN = 3, FUN = metric_function)
  })
}

#' @keywords internal
#' @importFrom stats predict
prediction_estimate <- function(mat_array, model) {
  # Create a data frame for prediction
  predict_df <- data.frame(
    Size = apply(mat_array, MARGIN = if(is.matrix(mat_array)) 1 else 3, FUN = nrow),
    MeanDensity = apply(mat_array, MARGIN = if(is.matrix(mat_array)) 1 else 3, FUN = network_density)
  )

  # Predict the duration using the model
  estimated_duration <- predict(model, predict_df)^3

  # Return the first predicted value as the estimated duration
  return(estimated_duration[1])
}

#' @keywords internal
convert_seconds_to_hms <- function(seconds) {
  # Calculate hours, minutes, and seconds
  hours <- seconds %/% 3600
  minutes <- (seconds %% 3600) %/% 60
  secs <- seconds %% 60

  # Build the time string conditionally
  time_parts <- c()

  if (hours > 0) {time_parts <- c(time_parts, paste(hours, "hours"))}

  if (minutes > 0 || hours > 0) { # Include minutes if hours are present, even if minutes are 0
    time_parts <- c(time_parts, paste(minutes, "mins"))
  }

  # Seconds are always included
  time_parts <- c(time_parts, paste(secs, "secs"))

  # Concatenate the parts into a single string
  time_string <- paste(time_parts, collapse = " ")

  return(time_string)
}

#' @keywords internal
#' @importFrom stats rnorm
calculate_metric_duration <- function(metric_name, matrices_array, user_benchmark, metrics_info) {
  model_path <- system.file("extdata", metrics_info$model_file[metric_name], package = "NetAnalyseR")
  model <- readRDS(model_path)
  number <- dim(matrices_array)[3]
  duration <- prediction_estimate(matrices_array, model) * user_benchmark * number
  duration <- round(duration,2)
  duration_time <- convert_seconds_to_hms(as.numeric(duration))
  tab_spacings <- list(
    global_efficiency_wei = "\t\t\t",
    global_clustering_coefficient_wei = "\t",
    characteristic_path_length = "\t\t",
    # Add more metrics as needed
    default = "\t\t\t" # Default tab spacing
  )
  tab_spacing <- ifelse(is.null(tab_spacings[[metric_name]]), tab_spacings[["default"]], tab_spacings[[metric_name]])
  cat(paste("\n\t", metric_name, tab_spacing, duration_time))
  return(duration)
}

# Define metrics and their models in a centralized way

#' @keywords internal
estimate_total_duration <- function(matrices_array, user_benchmark, valid_user_metrics) {

  metrics_info <- list(
    valid_user_metrics = c("global_efficiency_wei", "global_clustering_coefficient_wei", "characteristic_path_length", "local_efficiency_wei"),
    model_file = c(global_efficiency_wei = "global_efficiency_execution_time_model.rds",
                   local_efficiency_wei = "local_efficiency_execution_time_model.rds",
                   local_clustering_coefficient_wei = "local_clustering_coefficient_execution_time_model.rds",
                   global_clustering_coefficient_wei = "global_clustering_execution_time_model.rds",
                   characteristic_path_length = "characteristic_path_length_execution_time_model.rds")
  )
  total_time <- 0
  cat(paste("\nTime To Compute"))
  for(metric_name in valid_user_metrics[valid_user_metrics %in% metrics_info$valid_user_metrics]){
    duration <- calculate_metric_duration(metric_name, matrices_array, user_benchmark, metrics_info)
    total_time <- total_time + duration
  }
  if(any(! valid_user_metrics %in% metrics_info$valid_user_metrics)){
    fast_metrics <- valid_user_metrics[!valid_user_metrics %in% metrics_info$valid_user_metrics]
    cat("\n\n","The calculation time of the following metrics should be negligible: ",paste0(fast_metrics, collapse = ", "), sep = "")
  }

  if(total_time > 60){
    current_time <- Sys.time()
    estimated_finish <- current_time + as.difftime(total_time, units="secs")
    cat(paste("\n\nGlobal metric calculations should finish at:", format(estimated_finish, "%Y-%m-%d %H:%M:%S")))
  }
  if(total_time > 60*60){
    cat(paste("\n\nTime to compute is greater than 1 hour, consider thresholding matrices to decrease density and processing time."))
  }
  return(total_time)
}

