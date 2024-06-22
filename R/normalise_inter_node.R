normalise_inter_node <- function(mat, target){
  total <- inter_node(mat)
  return(mat * (target/total))
}
