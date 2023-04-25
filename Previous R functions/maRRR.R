# all functions for maRRR


####################### helper functions #######################

# generate V or V-shaped zero matrix
V_generator = function(nc=10,nr=100,nonzero=1){
  if(nonzero==1){
    V = matrix(rnorm(nc*nr),ncol = nc)
  }else if(nonzero==0){
    V = matrix(rep(0,nc*nr),ncol = nc)
  }
  return(V)
}

# generate Y matrix for each module 
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


# generate Y matrix without zero-like matrix
Y_combiner_nonzero = function(Yi_list,mod){
  Y = NULL
  for (i in 1:length(mod)) {
    if(mod[i]==1){
      Y = cbind(Y,Yi_list[[i]])
    }
  }
  return(Y)
}

# concatenate all matrices in a module, including zeros

mod_refine = function(X,mod,n_sample){
  ###########################################################################################################################
  ## @parameters: 
  ##    n_sample: numeric vector, number of samples for each cohort
  ##    X: matrix, full concatenated outcome matrix without zeros 
  ##    mod: numeric vector, current module information, e.g. c(1,0,1)
  ## ESTIMATION:
  ##    No estimation in the data generation function.
  ## @returns:
  ##    temp: matrix, concatenated outcome matrix for current module with zeros
  ###########################################################################################################################  
  
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

# replace all negative numbers with zeros in a vector

thres = function(vector){
  new_vec = NULL
  for (i in vector) {
    new_vec = c(new_vec,max(0,i))
  }
  return(new_vec)
}

# calculate current loss
# computationally slow

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

####################### main functions #######################

# generate ground truth B,Y,S for simulation use
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
  
  
  ###########################################################################################################################
  ## @parameters: 
  ##    seed: numeric, random seed 
  ##    n_sample: numeric vector, number of samples for each cohort
  ##    p_x: numeric, number of features in the observed data matrices
  ##    p_y: numeric, number of accompanying covariates 
  ##    R_b_true: numeric, true rank for simulated B matrices
  ##    R_s_true: numeric, true rank for simulated S matrices
  ##    correlation: numeric, Pearson's correlation aiming to achieve in the simulated Y matrices
  ##    orth_gen: logical, whether to orthogonalize the simulated Y matrices
  ##    Binvolved: logical, whether to generate covaraite effects, i.e. B and Y
  ##    Sinvolved: logical, whether to generate auxiliary variation structures, i.e. S
  ##    modules_B: list of numeric vectors, indicating which cohort has a covariate effect in each module,
  ##               i th "1" in the j th vector means the i th cohort is included in the j th module
  ##               e.g. c(1,1) means that it is a joint covariate effect of both Y1 and Y2;
  ##                    c(0,1) means that it is an individual covariate effect of Y2
  ##    modules_index_B: list of numeric vectors, indicating the index of Yi included
  ##               e.g. c(1,2) will corresponds to c(1,1) in "modules_B";
  ##               e.g. c(1,3,5) will corresponds to c(1,0,1,0,1) in "modules_B"
  ##    sd_B: numeric vector, standard deviation of covariate coefficient matrices (B)
  ##    n_mod_B: numeric, number of covariate-related modules included 
  ##    modules_S: list of numeric vectors, indicating which cohort has an auxiliary effect in each module,
  ##               i th "1" in the j th vector means the i th cohort is included in the j th module
  ##               e.g. c(1,1) means that it is a joint auxiliary structure [S1,S2];
  ##                    c(0,1) means that it is an individual auxiliary structure S2
  ##    modules_index_S: list of numeric vectors, indicating the index of Si included
  ##               e.g. c(1,2) will corresponds to c(1,1) in "modules_S";
  ##               e.g. c(1,3,5) will corresponds to c(1,0,1,0,1) in "modules_S"
  ##    sd_S: numeric vector, standard deviation of auxiliary structure matrices (S)
  ##    n_mod_S: numeric, number of covariate-unrelated modules included 
  ## ESTIMATION:
  ##    No estimation in the data generation function.
  ## @returns:
  ##    X_tot: matrix (#features x #samples), the concatenated version of outcome matrices
  ##           e.g. [X1,X2,X3] if there are three cohorts in total
  ##    Y_org_list: list of matrices, covariate matrices for each cohort
  ##    B_list: list of matrices, covariate effects of all modules
  ##    S_list: list of matrices, auxiliary structures of all modules
  ###########################################################################################################################  
  
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
  
  return(list(X_tot,Y_org_list,B_list,S_list))
}


# initialize starting values for estimation
# # of B_s = # of modules, can be initialized to be the same

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
  ###########################################################################################################################
  ## @parameters: 
  ##    seed: numeric, random seed 
  ##    n_sample: numeric vector, number of samples for each cohort
  ##    p_x: numeric, number of features in the observed data matrices
  ##    p_y: numeric, number of accompanying covariates 
  ##    R_b: numeric, upper bound for rank of simulated B matrices
  ##    R_s: numeric, upper bound for rank of simulated S matrices
  ##    Binvolved: logical, whether to generate covaraite effects, i.e. B and Y
  ##    Sinvolved: logical, whether to generate auxiliary variation structures, i.e. S
  ##    n_mod_B: numeric, number of covariate-related modules included 
  ##    n_mod_S: numeric, number of covariate-unrelated modules included 
  ## ESTIMATION:
  ##    No estimation in the initialization function.
  ## @returns:
  ##    UB_s_list: list of matrices, loadings for covariate effects of all modules
  ##    VB_s_list: list of matrices, scores for covariate effects of all modules
  ##    B_s_list: list of matrices, covariate effects of all modules
  ##    U_s_list: list of matrices, loadings for auxiliary structures of all modules
  ##    V_s_list: list of matrices, scores for auxiliary structures of all modules
  ##    S_s_list: list of matrices, auxiliary structures of all modules
  ###########################################################################################################################  
  
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



# optimization based on algorithm 1

ALS_UV = function(X_tot,Y_org_list,n_mod_B,B_s_list,UB_s_list,
                  VB_s_list,modules_B,n_mod_S=n_mod_B,S_s_list,U_s_list,V_s_list,
                  modules_S=modules_B,n_sample,lambdaBs,lambdaSs,bound,max_iter,
                  Binvolved,Sinvolved,orth_sol,loss_comp){
  ###########################################################################################################################
  ## @parameters: 
  ##    X_tot: matrix (#features x #samples), the concatenated version of outcome matrices
  ##           e.g. [X1,X2,X3] if there are three cohorts in total
  ##    Y_org_list: list of matrices, covariate matrices for each cohort
  ##    n_sample: numeric vector, number of samples for each cohort
  ##    orth_sol: logical, whether to orthogonalize covariate matrices before optimization
  ##    loss_comp: logical, whether to compute loss for each epoch
  ##    lambdaBs: numeric vector, penalties for each B module
  ##    lambdaSS: numeric vector, penalties for each S module
  ##    bound: numeric, convergence criteria for squared norm difference of B/S
  ##    max_iter: numeric, max number of epochs for optimization
  ##    Binvolved: logical, whether to generate covaraite effects, i.e. B and Y
  ##    Sinvolved: logical, whether to generate auxiliary variation structures, i.e. S
  ##    modules_B: list of numeric vectors, indicating which cohort has a covariate effect in each module,
  ##               i th "1" in the j th vector means the i th cohort is included in the j th module
  ##               e.g. c(1,1) means that it is a joint covariate effect of both Y1 and Y2;
  ##                    c(0,1) means that it is an individual covariate effect of Y2
  ##    n_mod_B: numeric, number of covariate-related modules included 
  ##    modules_S: list of numeric vectors, indicating which cohort has an auxiliary effect in each module,
  ##               i th "1" in the j th vector means the i th cohort is included in the j th module
  ##               e.g. c(1,1) means that it is a joint auxiliary structure [S1,S2];
  ##                    c(0,1) means that it is an individual auxiliary structure S2
  ##    n_mod_S: numeric, number of covariate-unrelated modules included 
  ## ESTIMATION:
  ##    UB_s_list: list of matrices, loadings for covariate effects of all modules
  ##    VB_s_list: list of matrices, scores for covariate effects of all modules
  ##    B_s_list: list of matrices, covariate effects of all modules
  ##    U_s_list: list of matrices, loadings for auxiliary structures of all modules
  ##    V_s_list: list of matrices, scores for auxiliary structures of all modules
  ##    S_s_list: list of matrices, auxiliary structures of all modules
  ## @returns:
  ##    time: numeric, total optimization time
  ##    B_s_list: list of matrices, estimated covariate effects of all modules
  ##    S_s_list: list of matrices, estimated auxiliary structures of all modules
  ##    loss_list: numeric vector, loss for each epoch
  ###########################################################################################################################  
  
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
  
  
  
  time = 0
  loss_list = NULL
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
    
    if(criter<bound){
      break
    }
    
  }
  
  # transform back to access the true scale of B
  for(i in 1:length(Y_list)){
    B_s_list[[i]] = B_s_list[[i]]%*%trans4B[[i]]
  }
  
  return(list(time,B_s_list,S_s_list,loss_list))
}

# optimization based on algorithm 2

ALS_BS = function(X_tot,Y_org_list,n_mod_B,B_s_list,UB_s_list,
                  VB_s_list,modules_B,n_mod_S=n_mod_B,S_s_list,U_s_list,V_s_list,
                  modules_S=modules_B,modules_index_S,
                  n_sample,lambdaBs,lambdaSs,bound,max_iter,
                  Binvolved,Sinvolved,orth_sol,loss_comp,col_index){
  ###########################################################################################################################
  ## @parameters: 
  ##    col_index: list of numeric vectors, i th vector represents 
  ##               the index for i th cohort in terms of Xtot
  ##    X_tot: matrix (#features x #samples), the concatenated version of outcome matrices
  ##           e.g. [X1,X2,X3] if there are three cohorts in total
  ##    Y_org_list: list of matrices, covariate matrices for each cohort
  ##    n_sample: numeric vector, number of samples for each cohort
  ##    orth_sol: logical, whether to orthogonalize covariate matrices before optimization
  ##              must set as "TRUE"
  ##    loss_comp: logical, whether to compute loss for each epoch
  ##    lambdaBs: numeric vector, penalties for each B module
  ##    lambdaSS: numeric vector, penalties for each S module
  ##    bound: numeric, convergence criteria for squared norm difference of B/S
  ##    max_iter: numeric, max number of epochs for optimization
  ##    Binvolved: logical, whether to generate covaraite effects, i.e. B and Y
  ##    Sinvolved: logical, whether to generate auxiliary variation structures, i.e. S
  ##    modules_B: list of numeric vectors, indicating which cohort has a covariate effect in each module,
  ##               i th "1" in the j th vector means the i th cohort is included in the j th module
  ##               e.g. c(1,1) means that it is a joint covariate effect of both Y1 and Y2;
  ##                    c(0,1) means that it is an individual covariate effect of Y2
  ##    n_mod_B: numeric, number of covariate-related modules included 
  ##    modules_S: list of numeric vectors, indicating which cohort has an auxiliary effect in each module,
  ##               i th "1" in the j th vector means the i th cohort is included in the j th module
  ##               e.g. c(1,1) means that it is a joint auxiliary structure [S1,S2];
  ##                    c(0,1) means that it is an individual auxiliary structure S2
  ##    modules_index_S: list of numeric vectors, indicating the index of Si included
  ##               e.g. c(1,2) will corresponds to c(1,1) in "modules_S";
  ##               e.g. c(1,3,5) will corresponds to c(1,0,1,0,1) in "modules_S"
  ##    n_mod_S: numeric, number of covariate-unrelated modules included 
  ## ESTIMATION:
  ##    UB_s_list: list of matrices, loadings for covariate effects of all modules
  ##    VB_s_list: list of matrices, scores for covariate effects of all modules
  ##    B_s_list: list of matrices, covariate effects of all modules
  ##    U_s_list: list of matrices, loadings for auxiliary structures of all modules
  ##    V_s_list: list of matrices, scores for auxiliary structures of all modules
  ##    S_s_list: list of matrices, auxiliary structures of all modules
  ## @returns:
  ##    time: numeric, total optimization time
  ##    B_s_list: list of matrices, estimated covariate effects of all modules
  ##    S_s_list: list of matrices, estimated auxiliary structures of all modules
  ##    loss_list: numeric vector, loss for each epoch
  ###########################################################################################################################  
  
  n_tot = dim(X_tot)[2]
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
  
  
  
  time = 0
  loss_list = NULL
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
    
    if(criter<bound){
      break
    }
    
  }
  
  # transform back to access the true scale of B
  for(i in 1:length(Y_list)){
    B_s_list[[i]] = B_s_list[[i]]%*%trans4B[[i]]
  }
  
  return(list(time,B_s_list,S_s_list,loss_list))
}



# missing imputation based on algorithm 2


ALS_BS_missing_impute = function(X_tot_incom,Y_org_list,n_mod_B,B_s_list,UB_s_list,
                                 VB_s_list,modules_B,n_mod_S=n_mod_B,S_s_list,U_s_list,V_s_list,
                                 modules_S=modules_B,modules_index_S,
                                 n_sample,lambdaBs,lambdaSs,bound,max_iter,
                                 Binvolved,Sinvolved,orth_sol,loss_comp,col_index){
    ###########################################################################################################################
  ## @parameters: 
  ##    col_index: list of numeric vectors, i th vector represents 
  ##               the index for i th cohort in terms of Xtot
  ##    X_tot_incom: matrix (#features x #samples), the concatenated version of 
  ##             outcome matrices with NAs, e.g. [X1,X2,X3] if there are three cohorts in total
  ##    Y_org_list: list of matrices, covariate matrices for each cohort
  ##    n_sample: numeric vector, number of samples for each cohort
  ##    orth_sol: logical, whether to orthogonalize covariate matrices before optimization
  ##              must set as "TRUE"
  ##    loss_comp: logical, whether to compute loss for each epoch
  ##    lambdaBs: numeric vector, penalties for each B module
  ##    lambdaSS: numeric vector, penalties for each S module
  ##    bound: numeric, convergence criteria for squared norm difference of B/S
  ##    max_iter: numeric, max number of epochs for optimization
  ##    Binvolved: logical, whether to generate covaraite effects, i.e. B and Y
  ##    Sinvolved: logical, whether to generate auxiliary variation structures, i.e. S
  ##    modules_B: list of numeric vectors, indicating which cohort has a covariate effect in each module,
  ##               i th "1" in the j th vector means the i th cohort is included in the j th module
  ##               e.g. c(1,1) means that it is a joint covariate effect of both Y1 and Y2;
  ##                    c(0,1) means that it is an individual covariate effect of Y2
  ##    n_mod_B: numeric, number of covariate-related modules included 
  ##    modules_S: list of numeric vectors, indicating which cohort has an auxiliary effect in each module,
  ##               i th "1" in the j th vector means the i th cohort is included in the j th module
  ##               e.g. c(1,1) means that it is a joint auxiliary structure [S1,S2];
  ##                    c(0,1) means that it is an individual auxiliary structure S2
  ##    modules_index_S: list of numeric vectors, indicating the index of Si included
  ##               e.g. c(1,2) will corresponds to c(1,1) in "modules_S";
  ##               e.g. c(1,3,5) will corresponds to c(1,0,1,0,1) in "modules_S"
  ##    n_mod_S: numeric, number of covariate-unrelated modules included 
  ## ESTIMATION:
  ##    UB_s_list: list of matrices, loadings for covariate effects of all modules
  ##    VB_s_list: list of matrices, scores for covariate effects of all modules
  ##    B_s_list: list of matrices, covariate effects of all modules
  ##    U_s_list: list of matrices, loadings for auxiliary structures of all modules
  ##    V_s_list: list of matrices, scores for auxiliary structures of all modules
  ##    S_s_list: list of matrices, auxiliary structures of all modules
  ## @returns:
  ##    X_tot: matrix, completed version of X_tot_incom
  ##    diff_list: numeric vector, sum of squared difference between the current 
  ##              and previous imputation
  ##    time: numeric, total optimization time
  ##    B_s_list: list of matrices, estimated covariate effects of all modules
  ##    S_s_list: list of matrices, estimated auxiliary structures of all modules
  ##    loss_list: numeric vector, loss for each epoch
  ###########################################################################################################################  
  
  
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

# missing imputation based on algorithm 1

ALS_UV_missing_impute = function(X_tot_incom,Y_org_list,n_mod_B,B_s_list,UB_s_list,
                                 VB_s_list,modules_B,n_mod_S=n_mod_B,S_s_list,U_s_list,V_s_list,
                                 modules_S=modules_B,n_sample,lambdaBs,lambdaSs,bound,max_iter,
                                 Binvolved,Sinvolved,orth_sol,loss_comp){
    ###########################################################################################################################
  ## @parameters: 
  ##    X_tot_incom: matrix (#features x #samples), the concatenated version of 
  ##             outcome matrices with NAs, e.g. [X1,X2,X3] if there are three cohorts in total
  ##    Y_org_list: list of matrices, covariate matrices for each cohort
  ##    n_sample: numeric vector, number of samples for each cohort
  ##    orth_sol: logical, whether to orthogonalize covariate matrices before optimization
  ##              must set as "TRUE"
  ##    loss_comp: logical, whether to compute loss for each epoch
  ##    lambdaBs: numeric vector, penalties for each B module
  ##    lambdaSS: numeric vector, penalties for each S module
  ##    bound: numeric, convergence criteria for squared norm difference of B/S
  ##    max_iter: numeric, max number of epochs for optimization
  ##    Binvolved: logical, whether to generate covaraite effects, i.e. B and Y
  ##    Sinvolved: logical, whether to generate auxiliary variation structures, i.e. S
  ##    modules_B: list of numeric vectors, indicating which cohort has a covariate effect in each module,
  ##               i th "1" in the j th vector means the i th cohort is included in the j th module
  ##               e.g. c(1,1) means that it is a joint covariate effect of both Y1 and Y2;
  ##                    c(0,1) means that it is an individual covariate effect of Y2
  ##    n_mod_B: numeric, number of covariate-related modules included 
  ##    modules_S: list of numeric vectors, indicating which cohort has an auxiliary effect in each module,
  ##               i th "1" in the j th vector means the i th cohort is included in the j th module
  ##               e.g. c(1,1) means that it is a joint auxiliary structure [S1,S2];
  ##                    c(0,1) means that it is an individual auxiliary structure S2
  ##    modules_index_S: list of numeric vectors, indicating the index of Si included
  ##               e.g. c(1,2) will corresponds to c(1,1) in "modules_S";
  ##               e.g. c(1,3,5) will corresponds to c(1,0,1,0,1) in "modules_S"
  ##    n_mod_S: numeric, number of covariate-unrelated modules included 
  ## ESTIMATION:
  ##    UB_s_list: list of matrices, loadings for covariate effects of all modules
  ##    VB_s_list: list of matrices, scores for covariate effects of all modules
  ##    B_s_list: list of matrices, covariate effects of all modules
  ##    U_s_list: list of matrices, loadings for auxiliary structures of all modules
  ##    V_s_list: list of matrices, scores for auxiliary structures of all modules
  ##    S_s_list: list of matrices, auxiliary structures of all modules
  ## @returns:
  ##    X_tot: matrix, completed version of X_tot_incom
  ##    diff_list: numeric vector, sum of squared difference between the current 
  ##              and previous imputation
  ##    time: numeric, total optimization time
  ##    B_s_list: list of matrices, estimated covariate effects of all modules
  ##    S_s_list: list of matrices, estimated auxiliary structures of all modules
  ##    loss_list: numeric vector, loss for each epoch
  ###########################################################################################################################  
  
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




