#' @keywords internal
.onAttach <- function(libname, pkgname) {
  .package_env <- new.env(parent = emptyenv())

  # Store the performance difference as a rds file
  performance_difference <- benchmark_performance()
}

benchmark_performance <- function(){
  size = 100
  repeats = 10
  sparsity = 0.2
  number = 10

  # model_path <- paste0(getwd(),"/inst/extdata/global_efficiency_execution_time_model.rds")
  model_path <- system.file("extdata", "local_clustering_execution_time_model.rds", package = "NetAnalyseR")
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

  saveRDS(proportion_difference, file = system.file("extdata", "performance_difference.rds", package = "NetAnalyseR"))
}

time_calc <- function(function_, n, ...) {
  elapsed <- system.time({function_(...)})[3]
  elapsed_average <- elapsed/n
  return(elapsed_average)
}

array_process <- function(functions, mat_array){
  lapply(list(functions), function(metric_function) {
    apply(mat_array, MARGIN = 3, FUN = metric_function)
  })
}
