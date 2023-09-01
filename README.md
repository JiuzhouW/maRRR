# Multiple Augmented Reduced Rank Regression (maRRR)

This folder contains R package and a simple demostration of usage for the **multiple augmented reduced rank regression** (maRRR) method presented in the article:

Jiuzhou Wang and Eric F. Lock. [Multiple Augmented Reduced Rank Regression for Pan-Cancer Analysis](https://arxiv.org/pdf/2308.16333.pdf), _Biometrics_, to appear, 2023.

This contains codes for fitting models by the two algorithms in Section 5.1, imputing missing values as in Section 6.2, generating penalties as in Section 4, and generating data as in Section 6.1.  This R package aims to predict one multi-cohort multi-variate data from another multi-cohort data. It is capable of doing missing data imputation as well. For an up-to-date and user-friendly R functions for the methods described see [here](https://github.com/JiuzhouW/maRRR/).  

One can use the method by loading R scripts `maRRR.R` in this [folder](https://github.com/JiuzhouW/maRRR/tree/main/Previous%20R%20functions). The folder contains two demos for complete data estimation and missing data imputation.

However, the R package is recommended :). 

Besides, we provide the data and all the model estimates that support our findings in the [paper](https://arxiv.org/pdf/2308.16333.pdf) as [Rdata file](https://www.dropbox.com/s/ub9zu5inxlbh6x5/30grps_50mods_est_0418_version2.RData?dl=0) with detailed [notation explanations](https://www.dropbox.com/s/af743f8ocucz7p6/readme_modelFit.txt?dl=0) and heatmaps for all module estimates in an [online file](https://www.dropbox.com/s/891gvukhymr4a4b/30grps_50mods_est_0509.csv?dl=0).

Now please let me walk through the method with you in 10 minutes.

## R package installation

This version of the `maRRR` package can be installed, directly from GitHub, using the devtools library:

```
install.packages("devtools")
library(devtools)
install_github("JiuzhouW/maRRR")
library(maRRR)
```

The details of the method can be found in the [manuscript](https://github.com/JiuzhouW/maRRR/blob/main/maRRR_manuscripts.pdf).


## Example 1: estimate covariate effects (B) and auxiliary structures (S)

All code chunks are arranged sequentially and are suggested to run in order.

### Model setup

Let's assume our true model is as follows:
```math
\displaylines{X_1 = B_1Y_1 + BY_1 + S_{1j} + S_1 (+ E_1) \\ X_2 = B_2Y_2 + BY_2 + S_{2j} + S_2 (+ E_2)}
```
The outcomes $X_1, X_2$ and covariates $Y_1,Y_2$ are observed from two cohorts. There are 3 covariate-related modules (one joint $B$ and two individual $B1,B2$). There are 3 auxiliary structures (one joint $[S_{1j}, S_{2j}]= U[V_{1j},V_{2j}]$, where $U$ is the joint loading matrix, and two individual $S_1,S_2$) to estimate. $E_1, E_2$ are random errors. Each of them has standard deviation as 1 by default.

If we concatenate all data from multiple cohorts, this model is equal to

```math
[X_1, X_2] = B[Y_1,Y_2] + B_1[Y_1, 0] + B_2[0, Y_2] + [S_{1j}, S_{2j}] + [S_1, 0] + [0, S_2] (+ [E_1, E_2])

```
 
We assume each cohort has 100 samples. Each sample has 100 outcomes and 10 covariates. Now we can generate the data.

```
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
               # by default, make the structure of S same as that of B, but free to make them different 
               # if Sinvolved = FALSE, this follow can be set as arbitrary
               modules_S = list(c(1,1),c(1,0),c(0,1)),
               modules_index_S = list(c(1,2),c(1),c(2)),
               sd_S = c(1,1,1),
               n_mod_S = 3)
```

### Prepare to estimate

To estimate all the effects of interest, we intialize our estimates with random numbers with the correct matrix structures.

```
# make sure "n_sample, p_x, p_y" match the those in "ini_gen"
# "n_mod_B" depends on the number of modules that you assume/detect
# if any, set "Binvolved = TRUE"
# so does "n_mod_S"

ini_dat = ini_gen(seed = 20230225,
                  n_sample = c(100,100),
                  p_x = 100,
                  p_y = 10,
                  R_b = 10,
                  R_s = 10,
                  n_mod_B = 3,
                  n_mod_S = 3,
                  Binvolved = TRUE,
                  Sinvolved = TRUE)
```

We generate best penalties based on random matrix theory mentioned in the paper but feel free to try different ones.

```
lambdaSs = lambda_S_gen(p_x=100,
                        modules_index_S=list(c(1,2),c(1),c(2)),
                        n_sample=c(100,100))
lambdaBs = lambda_B_gen(B_list = dat[[3]])

```
### maRRR version 1

Finally, let's run the `maRRR` method based on bilinear decomposition!

```
result_ALSUV = ALS_UV(X_tot=dat[[1]],
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
                      10^(-10),
                      max_iter=100,
                      Binvolved=TRUE,
                      Sinvolved=TRUE,
                      orth_sol=FALSE,
                      loss_comp=FALSE)

# examine estimations

R_b_true = 1
# relative mse for each B
mse_Bs_ALSUV = NULL
B_s_svd_list = list()
# estimated rank of B
rank_Bs_ALSUV = NULL
for (i in 1:3) {
  mse_Bs_ALSUV = c(mse_Bs_ALSUV,norm(result_ALSUV[[2]][[i]]-dat[[3]][[i]],"F")^2/norm(dat[[3]][[i]],"F")^2)
  B_s_svd_list[[i]] = svd(dat[[3]][[i]])
  rank_Bs_ALSUV = c(rank_Bs_ALSUV,sum(B_s_svd_list[[i]]$d>=0.1))
}

print(paste("The relative mean square errors for the joint covariate effect and two individual effects are", mse_Bs_ALSUV[1],", ", mse_Bs_ALSUV[2], ", and", mse_Bs_ALSUV[3], "."))
```

### maRRR version 2

There is another version of `maRRR` which is based on soft-thresholding. This requires semi-orthogonality in covariates (simply set `orth_sol=TRUE`) and specifying one more parameter `col_index`, which stores the sample index for each cohort.

```
# create the parameters for ALS_BS
n_samples = c(100,100)
col_index = list()
cur_index = 0
for (i in 1:length(n_samples)) {
  col_index[[i]] = (cur_index+1):(cur_index+n_samples[i])
  cur_index = cur_index + n_samples[i]
}

result_ALSBS = ALS_BS(X_tot=dat[[1]],
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
                      10^(-10),
                      max_iter=100,
                      Binvolved=TRUE,
                      Sinvolved=TRUE,
                      orth_sol=TRUE,
                      loss_comp=FALSE,
                      col_index=col_index)

# examine estimations

R_b_true = 1
# relative mse for each B
mse_Bs_ALSBS = NULL
B_s_svd_list = list()
# estimated rank of B
rank_Bs_ALSBS = NULL
for (i in 1:3) {
  mse_Bs_ALSBS = c(mse_Bs_ALSBS,norm(result_ALSBS[[2]][[i]]-dat[[3]][[i]],"F")^2/norm(dat[[3]][[i]],"F")^2)
  B_s_svd_list[[i]] = svd(dat[[3]][[i]])
  rank_Bs_ALSBS = c(rank_Bs_ALSBS,sum(B_s_svd_list[[i]]$d>=0.1))
}


print(paste("The relative mean square errors for the joint covariate effect and two individual effects are", mse_Bs_ALSBS[1],", ", mse_Bs_ALSBS[2], ", and", mse_Bs_ALSBS[3], "."))
```

Two version of `maRRR` methods should give you similar estiates!



## Example 2: missing imputation via maRRR

Just pour the incomplete outcome matrix and complete covariates information. It will return you the full complete outcome matrix! The model assumption is the same as it in Example 1 for simiplicity.

### Generate data and assign missingness
```
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
```

### parameter initialization

```
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
```

### Missing imputation based on maRRR version 1

```
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
```
The full complete outcome matrix is available at `Xtot_imputed`.

### Missing imputation based on maRRR version 2

```
# create the parameters for ALS_BS
n_samples = c(100,100)
col_index = list()
cur_index = 0
for (i in 1:length(n_samples)) {
  col_index[[i]] = (cur_index+1):(cur_index+n_samples[i])
  cur_index = cur_index + n_samples[i]
}

# impute with algorithm 2

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
```
The full complete outcome matrix is available at `Xtot_imputed`.
