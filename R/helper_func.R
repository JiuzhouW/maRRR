
#' @title Generate a random matrix or a zero matrix
#' @description A helper function. Generate V or V-shaped zero matrix.
#' @param nc A numeric, number of columns.
#' @param nr A numeric, number of rows.
#' @param nonzero A numeric, only 1 or 0 to indicate whether nonzero is needed or not.
#' @return V: A matrix.
V_generator = function(nc=10,nr=100,nonzero=1){
  if(nonzero==1){
    V = matrix(rnorm(nc*nr),ncol = nc)
  }else if(nonzero==0){
    V = matrix(rep(0,nc*nr),ncol = nc)
  }
  return(V)
}

#' @title Concatenate predictors from multiple cohorts
#' @description Generate multi-cohort Y matrix for each module from each cohort.
#' If a cohort is not included in the module, the corresponding entries will be zero.
#' @param Yi_list A list of matrices, covariate matrices for each cohort
#' @param mod A numeric vector, indicating which cohort has a covariate effect in each module,
#'               i th "1" in the j th vector means the i th cohort is included in the j th module
#'               e.g. c(1,1) means that it is a joint covariate effect of both Y1 and Y2;
#'                    c(0,1) means that it is an individual covariate effect of Y2
#' @return Y: A matrix.

Y_combiner = function(Yi_list,mod){
  Y = NULL
  for (i in 1:length(mod)) {
    if(mod[i]==1){
      Y = cbind(Y,Yi_list[[i]])
    }else if(mod[i]==0){
      msize = dim(Yi_list[[i]])
      Y = cbind(Y,matrix(rep(0,msize[1]*msize[2]),ncol = msize[2]))
    }
  }
  return(Y)
}

#' @title Concatenate predictors from multiple cohorts without zero-like matrix
#' @description Generate multi-cohort Y matrix for each module from each cohort.
#' If a cohort is not included in the module, there will be no corresponding entries.
#' @param Yi_list A list of matrices, covariate matrices for each cohort
#' @param mod A numeric vector, indicating which cohort has a covariate effect in each module,
#'               i th "1" in the j th vector means the i th cohort is included in the j th module
#'               e.g. c(1,1) means that it is a joint covariate effect of both Y1 and Y2;
#'                    c(0,1) means that it is an individual covariate effect of Y2
#' @return Y: A matrix.

Y_combiner_nonzero = function(Yi_list,mod){
  Y = NULL
  for (i in 1:length(mod)) {
    if(mod[i]==1){
      Y = cbind(Y,Yi_list[[i]])
    }
  }
  return(Y)
}



#' @title Concatenate all matrices in a module, including zeros to maintain the shape
#' @description Generate multi-cohort X matrix from each cohort.
#' If a cohort is not included in the module, the corresponding entries will be zero.
#' @param X A matrix, full concatenated outcome matrix without zeros. 
#' @param mod A numeric vector, indicating which cohort has a covariate effect in each module,
#'               i th "1" in the j th vector means the i th cohort is included in the j th module
#'               e.g. c(1,1) means that it is a joint covariate effect of both Y1 and Y2;
#'                    c(0,1) means that it is an individual covariate effect of Y2
#' @param n_sample A numeric vector, number of samples for each cohort.
#' @return temp: A matrix, concatenated outcome matrix for current module with zeros.


mod_refine = function(X,mod,n_sample){

  temp = NULL
  p_x = dim(X)[1]
  index_upper = 0
  index_lower = 1
  for(i in 1:length(mod)){
    index_upper = index_upper + n_sample[i]
    if(mod[i]==1){
      temp = cbind(temp, X[,index_lower:index_upper])
    }else if(mod[i]==0){
      temp = cbind(temp, matrix(rep(0,p_x*n_sample[i]),nrow = p_x) )
    }
    index_lower = index_lower + n_sample[i]
  }
  return(temp)
}

#' @title Replace all negative numbers with zeros in a vector
#' @param vector A vector
#' @return new_vec: A vector without negative numbers

thres = function(vector){
  new_vec = NULL
  for (i in vector) {
    new_vec = c(new_vec,max(0,i))
  }
  return(new_vec)
}

#' @title Calculate current
#' @description It can be computationally heavy if high dimension.

loss_calculator = function(X_tot,n_mod_B,n_mod_S,lambdaBs,lambdaSs,B_s_list,S_s_list,Y_list){
  X_res = X_tot 
  for(i in 1:n_mod_B){
    X_res = X_res - B_s_list[[i]]%*%Y_list[[i]] 
  }
  for(i in 1:n_mod_S){
    X_res = X_res - S_s_list[[i]]
  }
  
  loss = 0.5*sum(X_res^2) 
  for(i in 1:n_mod_B){
    loss = loss + lambdaBs[i]*sum(abs(svd(B_s_list[[i]])$d)) 
  }
  for(i in 1:n_mod_S){
    loss = loss + lambdaSs[i]*sum(abs(svd(S_s_list[[i]])$d))
  }
  return(loss)
}