# Generate best penalties for auxiliary structures S 
# based on the random matrix theory mentioned in the manuscripts

lambda_S_gen = function(p_x,modules_index_S,n_sample){
  all_lambda = c()
  for (i in modules_index_S) {
    temp = sqrt(sum(n_sample[i])) + sqrt(p_x)
    all_lambda = c(all_lambda,temp)
  }
  return(all_lambda)
}

# Generate best penalties for covariate effects B 
# based on the random matrix theory mentioned in the manuscripts

lambda_B_gen = function(B_list){
  all_lambda = c()
  for (i in B_list) {
    matrix_size = dim(i)
    temp = sqrt(matrix_size[1]) + sqrt(matrix_size[2])
    all_lambda = c(all_lambda,temp)
  }
  return(all_lambda)
}
