
#' @title Missing Imputation Based on Multiple Augmented Reduced Rank Regression via Soft-thresholding
#' @description  Impute the missing values of multi-cohort data with the aid of another multi-cohort data based on
#' alternating least square and soft-thresholding. It requires semi-orthogonality of predictors.
#' Optimization based on algorithm 2 from maRRR paper.
#' @author Jiuzhou Wang 
#' @param X_tot_incom A matrix (#features x #samples), the concatenated version of 
#'             outcome matrices with NAs, e.g. [X1,X2,X3] if there are three cohorts in total.
#' @param Y_org_list A list of matrices, covariate matrices for each cohort.
#' @param n_sample A numeric vector, number of samples for each cohort.
#' @param col_index A list of numeric vectors, i th vector represents 
#'               the index for i th cohort in terms of Xtot.
#' @param orth_sol A logical, whether to orthogonalize covariate matrices before optimization
#'              must set as "TRUE".
#' @param loss_comp A logical, whether to compute loss for each epoch.
#' @param lambdaBs A numeric vector, penalties for each B module.
#' @param lambdaSS A numeric vector, penalties for each S module.
#' @param bound A numeric, convergence criteria for squared norm difference of B/S.
#' @param max_iter A numeric, max number of epochs for optimization.
#' @param Binvolved A logical, whether to generate covaraite effects, i.e. B and Y.
#' @param Sinvolved A logical, whether to generate auxiliary variation structures, i.e. S.
#' @param modules_B A list of numeric vectors, indicating which cohort has a covariate effect in each module,
#'               i th "1" in the j th vector means the i th cohort is included in the j th module.
#'               e.g. c(1,1) means that it is a joint covariate effect of both Y1 and Y2;
#'                    c(0,1) means that it is an individual covariate effect of Y2.
#' @param n_mod_B A numeric, number of covariate-related modules included.
#' @param modules_S A list of numeric vectors, indicating which cohort has an auxiliary effect in each module,
#'               i th "1" in the j th vector means the i th cohort is included in the j th module.
#'               e.g. c(1,1) means that it is a joint auxiliary structure [S1,S2];
#'                    c(0,1) means that it is an individual auxiliary structure S2.
#' @param n_mod_S A numeric, number of covariate-unrelated modules included. Initial values for optimization.
#' @param UB_s_list A list of matrices, loadings for covariate effects of all modules. Initial values for optimization.
#' @param VB_s_list A list of matrices, scores for covariate effects of all modules. Initial values for optimization.
#' @param B_s_list A list of matrices, covariate effects of all modules. Initial values for optimization.
#' @param U_s_list A list of matrices, loadings for auxiliary structures of all modules. Initial values for optimization.
#' @param V_s_list A list of matrices, scores for auxiliary structures of all modules. Initial values for optimization.
#' @param S_s_list A list of matrices, auxiliary structures of all modules. Initial values for optimization.
#' @returns:
#' A list contains the following:
#' @item {X_tot}{A matrix, completed version of X_tot_incom.}
#' @item {diff_list}{A numeric vector, sum of squared difference between the current 
#'              and previous imputation.}
#' @item {time}{A numeric, total optimization time.} 
#' @item {B_s_list}{A list of matrices, estimated covariate effects of all modules.}
#' @item {S_s_list}{A list of matrices, estimated auxiliary structures of all modules.}
#' @item {loss_list}{A numeric vector, loss for each epoch.}
#' @examples
#' ALS_BS_missing_impute(X_tot,Y_org_list,n_mod_B,B_s_list,UB_s_list,
#' VB_s_list,modules_B,n_mod_S=n_mod_B,S_s_list,U_s_list,V_s_list,
#' modules_S=modules_B,n_sample,lambdaBs,lambdaSs,bound,max_iter,
#' Binvolved,Sinvolved,orth_sol,loss_comp,col_index)


ALS_BS_missing_impute = function(X_tot_incom,Y_org_list,n_mod_B,B_s_list,UB_s_list,
                                 VB_s_list,modules_B,n_mod_S=n_mod_B,S_s_list,U_s_list,V_s_list,
                                 modules_S=modules_B,modules_index_S,
                                 n_sample,lambdaBs,lambdaSs,bound,max_iter,
                                 Binvolved,Sinvolved,orth_sol,loss_comp,col_index){
  
  # get the matrix index of missing positions
  missing_pos = is.na(X_tot_incom)
  
  n_tot = dim(X_tot_incom)[2]
  R_b = dim(UB_s_list[[1]])[2]
  p_y = dim(VB_s_list[[1]])[1]
  p_x = dim(UB_s_list[[1]])[1]
  R_s = dim(U_s_list[[1]])[2]
  
  Y_list = list()
  
  for (i in 1:n_mod_B) {
    Y_list[[i]] = Y_combiner(Y_org_list,modules_B[[i]])
  }
  
  Y_SVD_list = list()
  
  for(i in 1:length(Y_list)){
    Y_SVD_list[[i]] = svd(Y_list[[i]])
  }
  
  # standardize Y and B if do not orthogonalize data
  # aim to make each row of Y has norm one
  if(!orth_sol){
    trans4B = list()
    for(i in 1:n_mod_B){
      reg_const = sqrt(rowSums(Y_list[[i]]^2))
      #B_list[[i]] = B_list[[i]]%*%diag(reg_const,ncol = length(reg_const))
      Y_list[[i]] = diag(1/reg_const,ncol = length(reg_const))%*%Y_list[[i]]
      # The estimated B should multiply the same scale
      trans4B[[i]] = diag(1/reg_const,ncol = length(reg_const))
    }
  }else if(orth_sol){
    trans4B = list()
    for(i in 1:length(Y_list)){
      Y_list[[i]] = t(Y_SVD_list[[i]]$v)
      trans4B[[i]] = diag(Y_SVD_list[[i]]$d^(-1))%*%t(Y_SVD_list[[i]]$u)
    }
  }
  
  # 
  X_tot = X_tot_incom
  
  # initialize some values for missing
  # fill the NA with row means
  for(i in 1:nrow(X_tot)){
    X_tot[i,is.na(X_tot[i,])] <- mean(X_tot[i,], na.rm = TRUE)
  }
  # current values that fill NA
  current_filling = X_tot[missing_pos]
  
  time = 0
  loss_list = NULL
  diff_list = NULL
  S_s_list_prev = S_s_list
  B_s_list_prev = B_s_list
  
  
  for (iter in 1:max_iter) {
    time = time + 1
    criter = 0
    
    #construct residual
    X_res = X_tot
    
    for(i in 1:n_mod_B){
      X_res = X_res - B_s_list[[i]]%*%Y_list[[i]] 
    }
    for (i in 1:n_mod_S) {
      X_res = X_res - S_s_list[[i]]
    }
    
    # update B
    if(Binvolved){
      for(i in 1:n_mod_B){
        X_res = X_res + B_s_list[[i]]%*%Y_list[[i]]
        #X_res does not need mod_refine since B does not have zero blocks
        XresY_svd = svd(X_res%*%t(Y_list[[i]]))
        B_s_list[[i]] = XresY_svd$u%*%diag(thres(XresY_svd$d-lambdaBs[i]))%*%
          t(XresY_svd$v)
        
        criter = criter + sum((B_s_list[[i]]-B_s_list_prev[[i]])^2)
        
        X_res = X_res - B_s_list[[i]]%*%Y_list[[i]]
      }
    }
    
    
    # update S
    if(Sinvolved){
      for (i in 1:n_mod_S) {
        
        X_res = X_res + S_s_list[[i]]
        
        # only SVD for parameters of interests
        index_list = c()
        for(cohort in modules_index_S[[i]]){
          index_list = c(index_list, col_index[[cohort]])
        }
        
        X_res2compose = X_res[,index_list]
        S_new_svd_2 = svd(X_res2compose)
        S_new2 = S_new_svd_2$u%*%diag(thres(S_new_svd_2$d - lambdaSs[i]))%*%t(S_new_svd_2$v)
        
        S_new2back = matrix(0, nrow = p_x, ncol = n_tot)
        
        S_new2back[,index_list] = S_new2
        S_s_list[[i]] = S_new2back
        
        X_res = X_res - S_s_list[[i]]
        
        criter = criter + sum((S_s_list[[i]]-S_s_list_prev[[i]])^2)
      }
    }
    
    if(loss_comp){
      loss_list = c(loss_list,loss_calculator(X_tot,n_mod_B,
                                              n_mod_S,lambdaBs,lambdaSs,B_s_list,S_s_list,Y_list))
    }
    S_s_list_prev = S_s_list
    B_s_list_prev = B_s_list
    
    # generate new X_tot
    X_tot_new = 0
    if(Binvolved){
      for (i in 1:n_mod_B) {
        X_tot_new = X_tot_new + B_s_list[[i]]%*%Y_list[[i]] 
      }
    }
    if(Sinvolved){
      for(i in 1:n_mod_S){
        X_tot_new = X_tot_new + S_s_list[[i]]
      }
    }
    
    # replace the NA with updates
    X_tot[missing_pos] = X_tot_new[missing_pos]
    # current difference with previous NA filling values
    diff_filling = sum((current_filling-X_tot_new[missing_pos])^2)
    diff_list = c(diff_list, diff_filling)
    # update the current filling values 
    current_filling = X_tot_new[missing_pos]
    
    if(max(criter,diff_filling)<bound){
      break
    }
    
  }
  
  # transform back to access the true scale of B
  for(i in 1:length(Y_list)){
    B_s_list[[i]] = B_s_list[[i]]%*%trans4B[[i]]
  }
  
  return(list(X_tot,diff_list,time,B_s_list,S_s_list,loss_list))
}


#' @title Missing Imputation Based on Multiple Augmented Reduced Rank Regression via Bilinear Decomposition
#' @description  Impute the missing values of multi-cohort data with the aid of another multi-cohort data based on
#' alternating least square and UV decomposition. 
#' Optimization based on algorithm 1 from maRRR paper.
#' @author Jiuzhou Wang 
#' @param X_tot_incom A matrix (#features x #samples), the concatenated version of 
#'             outcome matrices with NAs, e.g. [X1,X2,X3] if there are three cohorts in total.
#' @param Y_org_list A list of matrices, covariate matrices for each cohort.
#' @param n_sample A numeric vector, number of samples for each cohort.
#' @param orth_sol A logical, whether to orthogonalize covariate matrices before optimization
#'              must set as "TRUE".
#' @param loss_comp A logical, whether to compute loss for each epoch.
#' @param lambdaBs A numeric vector, penalties for each B module.
#' @param lambdaSS A numeric vector, penalties for each S module.
#' @param bound A numeric, convergence criteria for squared norm difference of B/S.
#' @param max_iter A numeric, max number of epochs for optimization.
#' @param Binvolved A logical, whether to generate covaraite effects, i.e. B and Y.
#' @param Sinvolved A logical, whether to generate auxiliary variation structures, i.e. S.
#' @param modules_B A list of numeric vectors, indicating which cohort has a covariate effect in each module,
#'               i th "1" in the j th vector means the i th cohort is included in the j th module.
#'               e.g. c(1,1) means that it is a joint covariate effect of both Y1 and Y2;
#'                    c(0,1) means that it is an individual covariate effect of Y2.
#' @param n_mod_B A numeric, number of covariate-related modules included.
#' @param modules_S A list of numeric vectors, indicating which cohort has an auxiliary effect in each module,
#'               i th "1" in the j th vector means the i th cohort is included in the j th module.
#'               e.g. c(1,1) means that it is a joint auxiliary structure [S1,S2];
#'                    c(0,1) means that it is an individual auxiliary structure S2.
#' @param n_mod_S A numeric, number of covariate-unrelated modules included. Initial values for optimization.
#' @param UB_s_list A list of matrices, loadings for covariate effects of all modules. Initial values for optimization.
#' @param VB_s_list A list of matrices, scores for covariate effects of all modules. Initial values for optimization.
#' @param B_s_list A list of matrices, covariate effects of all modules. Initial values for optimization.
#' @param U_s_list A list of matrices, loadings for auxiliary structures of all modules. Initial values for optimization.
#' @param V_s_list A list of matrices, scores for auxiliary structures of all modules. Initial values for optimization.
#' @param S_s_list A list of matrices, auxiliary structures of all modules. Initial values for optimization.
#' @returns:
#' A list contains the following:
#' @item {X_tot}{A matrix, completed version of X_tot_incom.}
#' @item {diff_list}{A numeric vector, sum of squared difference between the current 
#'              and previous imputation.}
#' @item {time}{A numeric, total optimization time.} 
#' @item {B_s_list}{A list of matrices, estimated covariate effects of all modules.}
#' @item {S_s_list}{A list of matrices, estimated auxiliary structures of all modules.}
#' @item {loss_list}{A numeric vector, loss for each epoch.}
#' @examples
#' ALS_UV_missing_impute(X_tot,Y_org_list,n_mod_B,B_s_list,UB_s_list,
#' VB_s_list,modules_B,n_mod_S=n_mod_B,S_s_list,U_s_list,V_s_list,
#' modules_S=modules_B,n_sample,lambdaBs,lambdaSs,bound,max_iter,
#' Binvolved,Sinvolved,orth_sol,loss_comp)

ALS_UV_missing_impute = function(X_tot_incom,Y_org_list,n_mod_B,B_s_list,UB_s_list,
                                 VB_s_list,modules_B,n_mod_S=n_mod_B,S_s_list,U_s_list,V_s_list,
                                 modules_S=modules_B,n_sample,lambdaBs,lambdaSs,bound,max_iter,
                                 Binvolved,Sinvolved,orth_sol,loss_comp){

  # get the matrix index of missing positions
  missing_pos = is.na(X_tot_incom)
  
  R_b = dim(UB_s_list[[1]])[2]
  p_y = dim(VB_s_list[[1]])[1]
  p_x = dim(UB_s_list[[1]])[1]
  R_s = dim(U_s_list[[1]])[2]
  
  Y_list = list()
  
  for (i in 1:n_mod_B) {
    Y_list[[i]] = Y_combiner(Y_org_list,modules_B[[i]])
  }
  
  Y_SVD_list = list()
  
  for(i in 1:length(Y_list)){
    Y_SVD_list[[i]] = svd(Y_list[[i]])
  }
  
  # standardize Y and B if do not orthogonalize data
  # aim to make each row of Y has norm one
  if(!orth_sol){
    trans4B = list()
    for(i in 1:n_mod_B){
      reg_const = sqrt(rowSums(Y_list[[i]]^2))
      #B_list[[i]] = B_list[[i]]%*%diag(reg_const,ncol = length(reg_const))
      Y_list[[i]] = diag(1/reg_const,ncol = length(reg_const))%*%Y_list[[i]]
      # The estimated B should multiply the same scale
      trans4B[[i]] = diag(1/reg_const,ncol = length(reg_const))
    }
  }else if(orth_sol){
    trans4B = list()
    for(i in 1:length(Y_list)){
      Y_list[[i]] = t(Y_SVD_list[[i]]$v)
      trans4B[[i]] = diag(Y_SVD_list[[i]]$d^(-1))%*%t(Y_SVD_list[[i]]$u)
    }
  }
  
  
  # 
  X_tot = X_tot_incom
  
  # initialize some values for missing
  # fill the NA with row means
  for(i in 1:nrow(X_tot)){
    X_tot[i,is.na(X_tot[i,])] <- mean(X_tot[i,], na.rm = TRUE)
  }
  # current values that fill NA
  current_filling = X_tot[missing_pos]
  
  time = 0
  loss_list = NULL
  diff_list = NULL
  S_s_list_prev = S_s_list
  B_s_list_prev = B_s_list
  
  
  for (iter in 1:max_iter) {
    time = time + 1
    criter = 0
    
    #construct residual
    X_res = X_tot
    
    for(i in 1:n_mod_B){
      X_res = X_res - B_s_list[[i]]%*%Y_list[[i]] 
    }
    for (i in 1:n_mod_S) {
      X_res = X_res - S_s_list[[i]]
    }
    
    # update B
    if(Binvolved){
      for(i in 1:n_mod_B){
        X_res = X_res + B_s_list[[i]]%*%Y_list[[i]]
        UB_s_list[[i]] = X_res%*%t(Y_list[[i]])%*%VB_s_list[[i]]%*%solve(t(VB_s_list[[i]])%*%
                                                                           crossprod(t(Y_list[[i]]))%*%VB_s_list[[i]]+lambdaBs[i]*diag(R_b))
        
        inv = solve(kronecker(crossprod(UB_s_list[[i]]),crossprod(t(Y_list[[i]]))) + lambdaBs[i]*diag(R_b*p_y))
        V_s_vec = inv%*%as.vector(Y_list[[i]]%*%t(X_res)%*%UB_s_list[[i]])
        VB_s_list[[i]] = matrix(V_s_vec, ncol = R_b, byrow = FALSE)
        
        
        B_s_list[[i]] = UB_s_list[[i]]%*%t(VB_s_list[[i]])
        
        criter = criter + sum((B_s_list[[i]]-B_s_list_prev[[i]])^2)
        
        X_res = X_res - B_s_list[[i]]%*%Y_list[[i]]
      }
    }
    
    
    # update S
    if(Sinvolved){
      for (i in 1:n_mod_S) {
        
        X_res = X_res + S_s_list[[i]]
        X_res_re = mod_refine(X_res,mod = modules_S[[i]],n_sample = n_sample)
        
        # refine X_s to satisfy the module structure, ie. some parts to be 0
        U_s_list[[i]] = X_res_re%*%V_s_list[[i]]%*%solve(t(V_s_list[[i]])%*%V_s_list[[i]] + lambdaSs[i]*diag(R_s))
        
        V_s_list[[i]] = t(X_res_re)%*%U_s_list[[i]]%*%solve(t(U_s_list[[i]])%*%U_s_list[[i]] + lambdaSs[i]*diag(R_s))
        S_s_list[[i]] = U_s_list[[i]]%*%t(V_s_list[[i]])
        
        X_res = X_res - S_s_list[[i]]
        
        criter = criter + sum((S_s_list[[i]]-S_s_list_prev[[i]])^2)
      }
    }
    
    if(loss_comp){
      loss_list = c(loss_list,loss_calculator(X_tot,n_mod_B,
                                              n_mod_S,lambdaBs,lambdaSs,B_s_list,S_s_list,Y_list))
    }
    S_s_list_prev = S_s_list
    B_s_list_prev = B_s_list
    
    # generate new X_tot
    X_tot_new = 0
    if(Binvolved){
      for (i in 1:n_mod_B) {
        X_tot_new = X_tot_new + B_s_list[[i]]%*%Y_list[[i]] 
      }
    }
    if(Sinvolved){
      for(i in 1:n_mod_S){
        X_tot_new = X_tot_new + S_s_list[[i]]
      }
    }
    
    # replace the NA with updates
    X_tot[missing_pos] = X_tot_new[missing_pos]
    # current difference with previous NA filling values
    diff_filling = sum((current_filling-X_tot_new[missing_pos])^2)
    diff_list = c(diff_list, diff_filling)
    # update the current filling values 
    current_filling = X_tot_new[missing_pos]
    
    if(max(criter,diff_filling)<bound){
      break
    }
    
  }
  
  # transform back to access the true scale of B
  for(i in 1:length(Y_list)){
    B_s_list[[i]] = B_s_list[[i]]%*%trans4B[[i]]
  }
  
  return(list(X_tot,diff_list,time,B_s_list,S_s_list,loss_list))
}

