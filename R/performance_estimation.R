## Generates matrix of specified dimentions where edges are individually generated
## from a normal distribution centered on zero with only absolute values being returned
mat_gen <- function(dim){
  mat <- matrix(abs(rnorm(dim*dim, 0, 100)),nrow = dim, ncol=dim,)
  mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]
  return(mat)
}

## Estimates the execution time of analyzing a set of networks given their size
estimate_duration <- function(mat_array, linear_model){
  if(!is.array(mat_array) & !is.matrix(mat_array)){
    stop("Input must be a matrix array or a single matrix")
  }
  # Handling both single matrices and arrays uniformly
  num_matrices <- if (is.matrix(mat_array)) 1 else dim(mat_array)[3]
  size <- apply(mat_array, MARGIN = if(is.matrix(mat_array)) 1 else 3, FUN = nrow)
  density <- apply(mat_array, MARGIN = if(is.matrix(mat_array)) 1 else 3, FUN = network_density)

  df <- data.frame(Size = size, Density = density)
  pred <- predict(linear_model, df)^3  # Assuming you indeed want to cube the prediction

  if (is.matrix(mat_array)) {
    return(pred)
  } else {
    return(sum(pred))
  }

}

benchmark_performance <- function(){
  size = 100
  repeats = 50

  estimated_duration <- c()
  execution_time <- c()
  for(i in 1:repeats){
    mat_array <- array(dim = c(size,size,repeats))
    for (j in 1:repeats) {
      # Assuming mat_gen() generates an m x p matrix
      mat_array[,,j] <- threshold_mat(mat_gen(size), 40)
    }
    model_path <- paste0(getwd(),"/inst/extdata/local_clustering_execution_time_model.rds")
    model <- readRDS(model_path)
    estimated_duration[i] <- estimate_duration(mat_array,model)
    start_time <- Sys.time()
    execution_time <- system.time(apply(mat_array, MARGIN = 3,local_clustering_coefficient_wei))[3]
  }
  # Returns proportion difference between execution time of machine and expected execution time on M2 Macbook.
  proportion_difference <- (mean(execution_time)-mean(estimated_duration))/ mean(execution_time)
  saveRDS(proportion_difference, file = paste0(getwd(),"/inst/extdata/performance_difference.rds"))

}

benchmark_performance()
model_path <- paste0(getwd(),"/inst/extdata/performance_difference.rds")
readRDS(model_path)
