# here we demonstrate how to use the maRRR function for missing imputation
# see detailed definitions for each function 

# import the required functions
source("maRRR.R")

# simulate complete data

## true model:
## X1 = B1Y1 + BY1 + S1j + S1 (+ E1)
## X2 = B2Y2 + BY2 + S2j + S2 (+ E2)
## this equal to
## [X1, X2] = B[Y1,Y2] + B1[Y1, 0] + B2[0, Y2]
##            + [S1j, S2j] + [S1, 0] + [0, S2] (+ [E1, E2])
## B is the shared covariate effect
## [S1j, S2j] = U[V1j V2j] is the joint auxiliary structure
## there are 3 covariate-related modules and 3 auxiliary structures
## each of them has standard deviation as 1

dat = data_gen(seed = 20230226,
               n_sample = c(100,100),
               p_x = 100,
               p_y = 10,
               R_b_true = 1,
               R_s_true = 5,
               correlation = 0,
               orth_gen = FALSE,
               Sinvolved = TRUE,
               # correlates to generation of B 
               modules_B = list(c(1,1),c(1,0),c(0,1)),
               modules_index_B = list(c(1,2),c(1),c(2)),
               sd_B = c(1,1,1),
               #n_mod_B = length(modules_B),
               n_mod_B = 3,
               # correlates to generation of S
               # by default, make the structure of S same as that of B
               # if Sinvolved = FALSE, this follow can be set as arbitrary
               modules_S = list(c(1,1),c(1,0),c(0,1)),
               modules_index_S = list(c(1,2),c(1),c(2)),
               sd_S = c(1,1,1),
               n_mod_S = 3)


# assign random missingness
# NOTE: if the data is incomplete, do not need to add extra missingness
Xtot = dat[[1]]
missing_prop = 0.05
n_features = dim(Xtot)[1]
n_tot_samples = dim(Xtot)[2]
n_entries = n_features*n_tot_samples

missing_index_entries = sample(1:n_entries, 
                               size = missing_prop*n_entries,
                               replace = FALSE)

# record the true values 
missing_index = missing_index_entries
true_missing = Xtot[missing_index]
Xtot_incom = Xtot
Xtot_incom[missing_index] = NA


# initialize estimation
# make sure "n_sample, p_x, p_y" match the those in "ini_gen"
# "n_mod_B" depends on the number of modules that you assume/detect
# if any, set "Binvolved = TRUE"
# so does "n_mod_S"

ini_dat = ini_gen(seed = 20230227,
                  n_sample = c(100,100),
                  p_x = 100,
                  p_y = 10,
                  R_b = 10,
                  R_s = 10,
                  n_mod_B = 3,
                  n_mod_S = 3,
                  Binvolved = TRUE,
                  Sinvolved = TRUE)

# generate best penalties based on detected modules

lambdaSs = lambda_S_gen(p_x=100,
                        modules_index_S=list(c(1,2),c(1),c(2)),
                        n_sample=c(100,100))
lambdaBs = lambda_B_gen(B_list = dat[[3]])


# impute with algorithm 1

maRRR_impute_result_UV = 
  ALS_UV_missing_impute(X_tot=Xtot_incom,
                        Y_org_list=dat[[2]],
                        n_mod_B=3,
                        B_s_list=ini_dat[[3]],
                        UB_s_list=ini_dat[[1]],
                        VB_s_list=ini_dat[[2]],
                        modules_B=list(c(1,1),c(1,0),c(0,1)),
                        n_mod_S=3,
                        S_s_list=ini_dat[[6]],
                        U_s_list=ini_dat[[4]],
                        V_s_list=ini_dat[[5]],
                        modules_S=list(c(1,1),c(1,0),c(0,1)),
                        n_sample=c(100,100),
                        lambdaBs,
                        lambdaSs,
                        10^(-1),
                        max_iter=30,
                        Binvolved=TRUE,
                        Sinvolved=TRUE,
                        orth_sol=FALSE,
                        loss_comp = FALSE)


Xtot_imputed = maRRR_impute_result_UV[[1]]
imputed_missing = Xtot_imputed[missing_index]

# relative MSE 
RSE_maRRR_UV = sum((true_missing - imputed_missing)^2)/sum((true_missing)^2)



# create the parameters for ALS_BS
n_samples = c(100,100)
col_index = list()
cur_index = 0
for (i in 1:length(n_samples)) {
  col_index[[i]] = (cur_index+1):(cur_index+n_samples[i])
  cur_index = cur_index + n_samples[i]
}

# impute with algorithm 1

maRRR_impute_result_BS = 
  ALS_BS_missing_impute(X_tot=Xtot_incom,
                        Y_org_list=dat[[2]],
                        n_mod_B=3,
                        B_s_list=ini_dat[[3]],
                        UB_s_list=ini_dat[[1]],
                        VB_s_list=ini_dat[[2]],
                        modules_B=list(c(1,1),c(1,0),c(0,1)),
                        n_mod_S=3,
                        S_s_list=ini_dat[[6]],
                        U_s_list=ini_dat[[4]],
                        V_s_list=ini_dat[[5]],
                        modules_S=list(c(1,1),c(1,0),c(0,1)),
                        modules_index_S = list(c(1,2),c(1),c(2)),
                        n_sample=c(100,100),
                        lambdaBs,
                        lambdaSs,
                        10^(-1),
                        max_iter=30,
                        Binvolved=TRUE,
                        Sinvolved=TRUE,
                        orth_sol=TRUE,
                        loss_comp = FALSE,
                        col_index=col_index)


Xtot_imputed = maRRR_impute_result_BS[[1]]
imputed_missing = Xtot_imputed[missing_index]

# relative MSE 
RSE_maRRR_BS = sum((true_missing - imputed_missing)^2)/sum((true_missing)^2)


## Both algorithms give similar estimations 
