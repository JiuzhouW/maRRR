
#' @title Generate penalty for S
#' @description  Generate best penalties for auxiliary structures S
#' based on random matrix theory from maRRR paper.
#' @author Jiuzhou Wang 
#' @param n_sample A numeric vector, number of samples for each cohort.
#' @param p_x A numeric, number of features of outcomes X.
#' @param modules_index_S list of numeric vectors, indicating the index of Si included
#'               e.g. c(1,2) will corresponds to c(1,1) in "modules_S";
#'               e.g. c(1,3,5) will corresponds to c(1,0,1,0,1) in "modules_S"
#' @return all_lambda: A vector of penalties for all auxiliary structures
#' @examples
#' lambdaSs = lambda_S_gen(p_x=100,
#' modules_index_S=list(c(1,2),c(1),c(2)),n_sample=c(100,100))

lambda_S_gen = function(p_x,modules_index_S,n_sample){
  all_lambda = c()
  for (i in modules_index_S) {
    temp = sqrt(sum(n_sample[i])) + sqrt(p_x)
    all_lambda = c(all_lambda,temp)
  }
  return(all_lambda)
}


#' @title Generate penalty for B
#' @description  Generate best penalties for covariate effects B 
#' based on random matrix theory from maRRR paper.
#' @author Jiuzhou Wang 
#' @param B_list A list of matrices, covariate effects of all modules. 
#' @return all_lambda: A vector of penalties for all auxiliary structures
#' @examples
#' lambdaBs = lambda_B_gen(B_list = dat[[3]])


lambda_B_gen = function(B_list){
  all_lambda = c()
  for (i in B_list) {
    matrix_size = dim(i)
    temp = sqrt(matrix_size[1]) + sqrt(matrix_size[2])
    all_lambda = c(all_lambda,temp)
  }
  return(all_lambda)
}
