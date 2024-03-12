#forecast_duration
mat_gen <- function(dim){
  mat <- matrix(abs(rnorm(dim*dim, 0, 100)),nrow = dim, ncol=dim,)
  mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]
  return(mat)
}

mat <- mat_gen(100)
stream_seq <- seq(max(ceiling(mat)))
density_curve <- rep(NA,length(stream_seq))
sparsity_of_int <- seq(0.05, 0.9, by = 0.05)
sparsity_indices <- numeric(length(sparsity_of_int)) # To store indices when density falls below thresholds

for(i in seq_along(stream_seq)){
  adjusted_mat <- mat - stream_seq[i]
  density_curve[i] <- network_density(adjusted_mat) # Calculate density for values above 0 after adjustment
  # Check if density falls below any sparsity of interest thresholds
  for(j in seq_along(sparsity_of_int)) {
    if(is.na(density_curve[i]) && density_curve[i] < sparsity_of_int[j]) {
      sparsity_indices[j] <- i # Record the first index where density falls below the threshold
    }
  }
}
plot(density_curve,type = "l")
network_density(mat_gen(5))
