#' @title Generate ground truth B,S,X,Y for simulation use
#' @param seed A numeric, random seed 
#' @param n_sample A numeric vector, number of samples for each cohort
#' @param p_x A numeric, number of features in the observed data matrices
#' @param p_y A numeric, number of accompanying covariates 
#' @param R_b_true A numeric, true rank for simulated B matrices
#' @param R_s_true A numeric, true rank for simulated S matrices
#' @param correlation A numeric, Pearson's correlation aiming to achieve in the simulated Y matrices
#' @param orth_gen A logical, whether to orthogonalize the simulated Y matrices
#' @param Binvolved A logical, whether to generate covaraite effects, i.e. B and Y
#' @param Sinvolved A logical, whether to generate auxiliary variation structures, i.e. S
#' @param modules_B A list of numeric vectors, indicating which cohort has a covariate effect in each module,
#'               i th "1" in the j th vector means the i th cohort is included in the j th module
#'               e.g. c(1,1) means that it is a joint covariate effect of both Y1 and Y2;
#'                    c(0,1) means that it is an individual covariate effect of Y2
#' @param modules_index_B A list of numeric vectors, indicating the index of Yi included
#'               e.g. c(1,2) will corresponds to c(1,1) in "modules_B";
#'               e.g. c(1,3,5) will corresponds to c(1,0,1,0,1) in "modules_B"
#' @param sd_B A numeric vector, standard deviation of covariate coefficient matrices (B)
#' @param n_mod_B A numeric, number of covariate-related modules included 
#' @param modules_S A list of numeric vectors, indicating which cohort has an auxiliary effect in each module,
#'               i th "1" in the j th vector means the i th cohort is included in the j th module
#'               e.g. c(1,1) means that it is a joint auxiliary structure [S1,S2];
#'                    c(0,1) means that it is an individual auxiliary structure S2
#' @param modules_index_S A list of numeric vectors, indicating the index of Si included
#'               e.g. c(1,2) will corresponds to c(1,1) in "modules_S";
#'               e.g. c(1,3,5) will corresponds to c(1,0,1,0,1) in "modules_S"
#' @param sd_S A numeric vector, standard deviation of auxiliary structure matrices (S)
#' @param n_mod_S A numeric, number of covariate-unrelated modules included 
#' @returns:
#' A list of a outcome matrix, covariate matices, covariate effects and auxiliary effects
#' \item{X_tot}{A matrix (#features x #samples), the concatenated version of outcome matrices
#'      e.g. [X1,X2,X3] if there are three cohorts in total}
#' \item{Y_org_list}{A list of matrices, covariate matrices for each cohort}
#' \item{B_list}{A list of matrices, covariate effects of all modules}
#' \item{S_list}{A list of matrices, auxiliary structures of all modules}


data_gen = function(seed = 1,
                    n_sample = c(100,100),
                    p_x = 100,
                    p_y = 10,
                    R_b_true = 1,
                    R_s_true = 5,
                    correlation = 0,
                    orth_gen = FALSE,
                    Binvolved = TRUE,
                    Sinvolved = TRUE,
                    modules_B = list(c(1,1),c(1,0),c(0,1)),
                    modules_index_B = list(c(1,2),c(1),c(2)),
                    sd_B = c(1,1,1),
                    n_mod_B = length(modules_B),
                    modules_S = modules_B,
                    modules_index_S = modules_index_B,
                    sd_S = sd_B,
                    n_mod_S = length(modules_S)){
  
  
  library(MASS)
  set.seed(seed)
  
  len_per_mod = length(modules_B[[1]])
  
  # generate Y only if we include B, but the data generating process should be different from solving
  
  
  # individual Y_i, regressors in each group
  Y_org_list = list()
  sd_Y_org = NULL
  
  if(correlation==0){
    for (i in 1:len_per_mod) {
      Y_org_list[[i]] = matrix(rnorm(n_sample[i]*p_y),ncol = n_sample[i])
      sd_Y_org = c(sd_Y_org,sd(Y_org_list[[i]]))
    }
  }else if(correlation>0){
    sigma = matrix(rep(correlation,p_y^2),ncol = p_y)+diag(1-correlation,p_y)
    for (i in 1:len_per_mod) {
      Y_org_list[[i]] = t(mvrnorm(n_sample[i],rep(0,p_y),Sigma = sigma))
      sd_Y_org = c(sd_Y_org,sd(Y_org_list[[i]]))
    }
  }
  
  
  
  # individually orthogonalize each Yi
  if(orth_gen){
    for(i in 1:len_per_mod){
      Y_org_list[[i]] = t(svd(Y_org_list[[i]])$v)
      # important step
      # rescale orthogonalized Y to have the same sd(Y)
      Y_org_list[[i]] = Y_org_list[[i]]/sd(Y_org_list[[i]])*sd_Y_org[i]
    }
  }
  
  #combined modulized Y to generate X_tot
  # another way to handle this is by generate X sperately and combine together
  Y_gen_list = list()
  
  for (i in 1:n_mod_B) {
    Y_gen_list[[i]] = Y_combiner(Y_org_list,modules_B[[i]])
  }
  
  
  #initialize UB and VB, therefore B
  
  UB_list = list()
  VB_list = list()
  B_list = list()
  
  if(Sinvolved){
    for(i in 1:n_mod_B){
      
      UB_list[[i]] = matrix(rnorm(p_x*R_b_true),ncol = R_b_true)
      VB_list[[i]] = matrix(rnorm(p_y*R_b_true),ncol = R_b_true)
      B_list[[i]] = sd_B[i]*UB_list[[i]]%*%t(VB_list[[i]])/sd(UB_list[[i]]%*%t(VB_list[[i]])%*%Y_gen_list[[i]])
    }
  }else if(!Sinvolved){
    for(i in 1:n_mod_B){
      
      UB_list[[i]] = matrix(rnorm(p_x*R_b_true),ncol = R_b_true)
      VB_list[[i]] = matrix(rnorm(p_y*R_b_true),ncol = R_b_true)
      B_list[[i]] = sd_B[i]*UB_list[[i]]%*%t(VB_list[[i]])/sd(UB_list[[i]]%*%t(VB_list[[i]]))
      #B_list[[i]] = sd_B[i]*UB_list[[i]]%*%t(VB_list[[i]])
    }
  }
  
  # check the correlation after standardization
  # if orth_gen == 0, then it corresponds with the original since multiply a number to the Yi does not affect corr
  # orthonormalize Y will automatically make the correlation to be 0
  #cor(t(Y_org_list[[1]]))
  #cor(t(Y_list[[1]][,1:n_sample[1]]))
  
  # generate S if involved
  
  if(Sinvolved){
    #initialize U and V, therefore S
    
    U_list = list()
    V_list = list()
    S_list = list()
    
    for(i in 1:n_mod_S){
      #U_s is universal for one module
      U_list[[i]] = matrix(rnorm(p_x*R_s_true),ncol = R_s_true)
      #V_s may contain parts of zero matrices
      temp = NULL
      for (j in 1:length(modules_S[[i]])) {
        print(modules_S[[i]][j])
        temp = rbind(temp,V_generator(nc=R_s_true,nr=n_sample[j],nonzero=modules_S[[i]][j]))
      }
      V_list[[i]] = temp
      S_list[[i]] = sd_S[i]*U_list[[i]]%*%t(V_list[[i]])/sd(U_list[[i]]%*%t(V_list[[i]]))
    }
  }else if(!Sinvolved){
    
    U_list = list()
    V_list = list()
    S_list = list()
    
    for(i in 1:n_mod_S){
      U_list[[i]] = matrix(rep(0,p_x*R_s_true),ncol = R_s_true)
      V_list[[i]] = matrix(rep(0,sum(n_sample)*R_s_true),ncol = R_s_true)
      S_list[[i]] = sd_S[i]*U_list[[i]]%*%t(V_list[[i]])
    }
  }
  
  
  X_tot = 0
  for (i in 1:n_mod_B) {
    X_tot = X_tot + B_list[[i]]%*%Y_gen_list[[i]] 
  }
  for(i in 1:n_mod_S){
    X_tot = X_tot + S_list[[i]]
  }
  X_tot = X_tot + matrix(rnorm(p_x*sum(n_sample)),nrow = p_x )
  
  return(list(X_tot=X_tot,Y_org_list=Y_org_list,B_list=B_list,S_list=S_list))
}



#' @title Initialize starting values for optimization
#' @param seed A numeric, random seed 
#' @param n_sample A numeric vector, number of samples for each cohort
#' @param p_x A numeric, number of features in the observed data matrices
#' @param p_y A numeric, number of accompanying covariates 
#' @param R_b A numeric, upper bound for rank of simulated B matrices
#' @param R_s A numeric, upper bound for rank of simulated S matrices
#' @param Binvolved A logical, whether to generate covaraite effects, i.e. B and Y
#' @param Sinvolved A logical, whether to generate auxiliary variation structures, i.e. S
#' @param n_mod_B A numeric, number of covariate-related modules included 
#' @param n_mod_S A numeric, number of covariate-unrelated modules included 
#' @returns:
#' A list of initialized estimates.
#' \item{UB_s_list}{A list of matrices, loadings for covariate effects of all modules}
#' \item{VB_s_list}{A list of matrices, scores for covariate effects of all modules}
#' \item{B_s_list}{A list of matrices, covariate effects of all modules}
#' \item{U_s_list}{A list of matrices, loadings for auxiliary structures of all modules}
#' \item{V_s_list}{A list of matrices, scores for auxiliary structures of all modules}
#' \item{S_s_list}{A list of matrices, auxiliary structures of all modules}



ini_gen = function(seed = 1,
                   n_sample = c(100,100),
                   p_x = 100,
                   p_y = 10,
                   R_b = 10,
                   R_s = 10,
                   n_mod_B = 3,
                   n_mod_S = n_mod_B,
                   Binvolved = TRUE,
                   Sinvolved = TRUE){

  set.seed(seed)
  
  UB_s_list = list()
  VB_s_list = list()
  B_s_list = list()
  U_s_list = list()
  V_s_list = list()
  S_s_list = list()
  
  if(Binvolved){
    for(i in 1:n_mod_B){
      UB_s_list[[i]] = matrix(rnorm(p_x*R_b),ncol = R_b)
      VB_s_list[[i]] = matrix(rnorm(p_y*R_b),ncol = R_b)
      B_s_list[[i]] = UB_s_list[[i]] %*% t(VB_s_list[[i]])
    }
  }else if(!Binvolved){
    for(i in 1:n_mod_B){
      #U_s is universal for one module
      UB_s_list[[i]] = matrix(rep(0,p_x*R_b),ncol = R_b)
      VB_s_list[[i]] = matrix(rep(0,p_y*R_b),ncol = R_b)
      B_s_list[[i]] = UB_s_list[[i]] %*% t(VB_s_list[[i]])
    }
  }
  
  
  
  if(Sinvolved){
    for(i in 1:n_mod_S){
      #U_s is universal for one module
      U_s_list[[i]] = matrix(rnorm(p_x*R_s),ncol = R_s)
      V_s_list[[i]] = matrix(rnorm(sum(n_sample)*R_s),ncol = R_s)
      S_s_list[[i]] = U_s_list[[i]]%*%t(V_s_list[[i]])
    }
  }else if(!Sinvolved){
    for(i in 1:n_mod_S){
      #U_s is universal for one module
      U_s_list[[i]] = matrix(rep(0,p_x*R_s),ncol = R_s)
      V_s_list[[i]] = matrix(rep(0,sum(n_sample)*R_s),ncol = R_s)
      S_s_list[[i]] = U_s_list[[i]]%*%t(V_s_list[[i]])
    }
  }
  
  return(list(UB_s_list,VB_s_list,B_s_list,U_s_list,V_s_list,S_s_list))
  
}