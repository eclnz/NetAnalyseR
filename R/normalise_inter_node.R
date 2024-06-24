normalise_inter_node <- function(mat, target, validate = TRUE){
  if(validate){validate_matrix(mat)}
  total <- inter_node(mat)
  return(mat * (target/total))
}
