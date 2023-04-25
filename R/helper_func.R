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

# generate multi-cohort Y matrix for each module from each cohort
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