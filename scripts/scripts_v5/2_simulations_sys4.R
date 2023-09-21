## ----load_lib, include = FALSE, warning=FALSE, message=FALSE, echo=FALSE-------------------------------------------------------------
library(here)
library(data.table)
library(dplyr)
library(tidyr)
library(foreach)

library(stringr)
library(glmnet)

library(origami)
library(hal9001)
library(tictoc)

library(R.utils)

library(pROC)

library(ggplot2)
library(cowplot)
library(gridExtra)
library(grid)

library(mvtnorm)

## ----setup, include = FALSE----------------------------------------------------------------------------------------------------------
plotFolder <- here("results","images", "v5")
if(!file.exists(plotFolder)) dir.create(plotFolder,recursive=TRUE)

knitr::opts_chunk$set(
  cache=FALSE, autodep=FALSE, warning=FALSE, message=FALSE, echo=FALSE,
  results = 'markup', dev='png', dpi=150, fig.align = "center", fig.path=paste0(plotFolder, "/"),
  cache.path=".cache/",
  duplicate.label="allow"
)

source(here("scripts", "scripts_v5", "1_hal_functions.R"))
source(here("scripts", "scripts_v5", "1_simu_functions.R"))



## ----check_sys4----------------------------------------------------------------------------------------------------------------------
generate_data_4 <- function(n, a=NA){
  # exogenous variables
  U_W_sigma = c(1, 1, 1.5, 0.5, 0.3, 1, 1, 5, 0.5, 2)
  
  U_W_r = matrix(c(1,    0.8,  0.9, 0.2, 0.1, 0,   0,   0.3, 0.1, 0.96,
                   0.8,  1,    0.8, 0.3, 0.1, 0,   0,   0.8, 0.1, 0.8,   
                   0.9,  0.8,  1,   0.1, 0.5, 0.7, 0.1, 0.7, 0.8, 0.9, 
                   0.2,  0.3,  0.1, 1,   0.5, 0.3, 0,   0,   1,   1,   
                   0.1,  0.1,  0.5, 0.5, 1,   0,   0,   0,   0.1, 0.1,
                   0,    0,    0.7, 0.3, 0,   1,   0.8, 0,   0.4, 0.1,
                   0,    0,    0.1, 0,   0,   0.8, 1,   0.9, 0.8, 0.5,
                   0.3,  0.8,  0.7, 0,   0,   0,   0.9, 1,   0.7, 0.9,
                   0.1,  0.1,  0.8, 1,   0.1, 0.4, 0.8, 0.7, 1,   0.8,
                   0.96, 0.8,  0.9, 1,   0.1, 0.1, 0.5, 0.9, 0.8, 1), nrow = 10)
  
  U_W_var = matrix(nrow = 10, ncol = 10)
  
  for (i in 1:10) {
    for (j in 1:10) {
      U_W_var[i,j] = U_W_r[i,j] * U_W_sigma[i] * U_W_sigma[j]
    }
  }
  
  U_W <- rmvnorm(n, mean = c(0, 0, 5, 1, 1, 0.3, 1, 0, 3, 2), sigma = U_W_var)
  U_A <- rnorm(n, 0, 1)
  U_Y <- runif(n, 0, 1)
  
  # endogenous variables
  W <- data.frame(U_W)
  colnames(W) <- paste0('W', 1:ncol(W))
  
  if(is.na(a)){
    A <-  (0.1*W$W1 + 0.2*W$W2  + 0.5*W$W3 + 0.15*W$W4 - 0.05*W$W5 - 0.01*W$W3*W$W5 + U_A)
    A[A<=0] = 0
    A[A>=5] = 5
  } else {
    A <- rep(a, n)
  }

  Y <- as.numeric(U_Y < plogis(-7 - W$W1 + 2*W$W2 - 1*W$W1 * W$W3 - 0.5*W$W4 - 0.2 * W$W6^2 + 0.01 * W$W8 *W$W10 + 4*A + 0.5*A*W$W2))
  mean(Y)
  
  # data frame
  O <- W
  O$A = A
  O$Y = Y
  return(O)
}


## ----true_psi_sys4-------------------------------------------------------------------------------------------------------------------
# Getting trul value of psi
#------------------------------------------------------------------------------------
# 
a_vec <- seq(0,5,0.1)
psi0_a <- c()
psi0_a <- c()

N = 1e+07

for (i in 1:length(a_vec)) {
  a <- a_vec[i]

  data_a <- generate_data_4(n=N, a=a)
  psi0_a[i] <- mean(data_a$Y) # - data_0$Y
}

psi0_line <- data.frame(a=a_vec, psi0 = psi0_a)
plot(psi0_line$a, psi0_line$psi0)

eval_points = seq(0,5,0.5)
psi0_pnt <- psi0_line[psi0_line$a %in% eval_points,]

save.image(file=here("data", "rdata", "02_simu_V5_sys4_psi0.RData"))



#======================================================================
# Simulations

## -----------------------------------------------------------------------------------------------------------------------
load(here("data", "rdata", "02_simu_V5_sys4_psi0.RData"))

source(here("scripts", "scripts_v5", "1_hal_functions.R"))
source(here("scripts", "scripts_v5", "1_simu_functions.R"))

#======================================================================
n = 200

## -----------------------------------------------------------------------------------------------------------------------
set.seed(123)

results <- run_simu_rep(generate_data_4, eval_points, y_type = "binomial", n=n, rounds=500, return_all_rslts=T)
save.image(file=here("data", "rdata", "02_simu_V5_sys4_200.RData"))

## -----------------------------------------------------------------------------------------------------------------------
rm(results)

set.seed(123)

results_0 <- run_simu_rep(generate_data_4, eval_points, y_type = "binomial", n=n, rounds=500, return_all_rslts=T, defualt_setting = T)
save.image(file=here("data", "rdata", "02_simu_V5_sys4_200_default.RData"))

## -----------------------------------------------------------------------------------------------------------------------
rm(results_0)

set.seed(123)
results_grid <- run_simu_scaled_rep(generate_data_4, eval_points, y_type = "binomial", n=n, rounds=500, return_all_rslts=T)

save.image(file=here("data", "rdata", "02_simu_V5_sys4_200_grid.RData"))

## -----------------------------------------------------------------------------------------------------------------------
rm(results_grid)

set.seed(123)
results_adapt <- run_simu_smoothness_adaptive_HAL_rep(generate_data_4, eval_points, y_type = "binomial", n=n, rounds=500, return_all_rslts=T)
save.image(file=here("data", "rdata", "02_simu_V5_sys4_200_adapt.RData"))



## -----------------------------------------------------------------------------------------------------------------------
rm(results_adapt)

source(here("scripts", "scripts_v5", "1_simu_functions_noHAL.R"))

set.seed(123)
results_gam <- run_simu_gam_poly_rep(generate_data_4, eval_points, y_type = "binomial", n=n, rounds=500, return_all_rslts=T, method = "GAM")
save.image(file=here("data", "rdata", "02_simu_V5_sys4_200_GAM.RData"))

rm(results_gam)
set.seed(123)
results_poly <- run_simu_gam_poly_rep(generate_data_4, eval_points, y_type = "binomial", n=n, rounds=500, return_all_rslts=T, method = "POLY")
save.image(file=here("data", "rdata", "02_simu_V5_sys4_200_poly.RData"))





#======================================================================
rm(results_poly)

n = 500

## -----------------------------------------------------------------------------------------------------------------------
set.seed(123)

results <- run_simu_rep(generate_data_4, eval_points, y_type = "binomial", n=n, rounds=500, return_all_rslts=T)
save.image(file=here("data", "rdata", "02_simu_V5_sys4_500.RData"))

## -----------------------------------------------------------------------------------------------------------------------
rm(results)

set.seed(123)

results_0 <- run_simu_rep(generate_data_4, eval_points, y_type = "binomial", n=n, rounds=500, return_all_rslts=T, defualt_setting = T)
save.image(file=here("data", "rdata", "02_simu_V5_sys4_500_default.RData"))

## -----------------------------------------------------------------------------------------------------------------------
rm(results_0)

set.seed(123)
results_grid <- run_simu_scaled_rep(generate_data_4, eval_points, y_type = "binomial", n=n, rounds=500, return_all_rslts=T)

save.image(file=here("data", "rdata", "02_simu_V5_sys4_500_grid.RData"))

## -----------------------------------------------------------------------------------------------------------------------
rm(results_grid)

set.seed(123)
results_adapt <- run_simu_smoothness_adaptive_HAL_rep(generate_data_4, eval_points, y_type = "binomial", n=n, rounds=500, return_all_rslts=T)
save.image(file=here("data", "rdata", "02_simu_V5_sys4_500_adapt.RData"))



## -----------------------------------------------------------------------------------------------------------------------
rm(results_adapt)

source(here("scripts", "scripts_v5", "1_simu_functions_noHAL.R"))

set.seed(123)
results_gam <- run_simu_gam_poly_rep(generate_data_4, eval_points, y_type = "binomial", n=n, rounds=500, return_all_rslts=T, method = "GAM")
save.image(file=here("data", "rdata", "02_simu_V5_sys4_500_GAM.RData"))

rm(results_gam)
set.seed(123)
results_poly <- run_simu_gam_poly_rep(generate_data_4, eval_points, y_type = "binomial", n=n, rounds=500, return_all_rslts=T, method = "POLY")
save.image(file=here("data", "rdata", "02_simu_V5_sys4_500_poly.RData"))





#======================================================================
rm(results_poly)

n = 1000

## -----------------------------------------------------------------------------------------------------------------------
set.seed(123)

results <- run_simu_rep(generate_data_4, eval_points, y_type = "binomial", n=n, rounds=500, return_all_rslts=T)
save.image(file=here("data", "rdata", "02_simu_V5_sys4_1000.RData"))

## -----------------------------------------------------------------------------------------------------------------------
rm(results)

set.seed(123)

results_0 <- run_simu_rep(generate_data_4, eval_points, y_type = "binomial", n=n, rounds=500, return_all_rslts=T, defualt_setting = T)
save.image(file=here("data", "rdata", "02_simu_V5_sys4_1000_default.RData"))

## -----------------------------------------------------------------------------------------------------------------------
rm(results_0)

set.seed(123)
results_grid <- run_simu_scaled_rep(generate_data_4, eval_points, y_type = "binomial", n=n, rounds=500, return_all_rslts=T)

save.image(file=here("data", "rdata", "02_simu_V5_sys4_1000_grid.RData"))

## -----------------------------------------------------------------------------------------------------------------------
rm(results_grid)

set.seed(123)
results_adapt <- run_simu_smoothness_adaptive_HAL_rep(generate_data_4, eval_points, y_type = "binomial", n=n, rounds=500, return_all_rslts=T)
save.image(file=here("data", "rdata", "02_simu_V5_sys4_1000_adapt.RData"))



## -----------------------------------------------------------------------------------------------------------------------
rm(results_adapt)

source(here("scripts", "scripts_v5", "1_simu_functions_noHAL.R"))

set.seed(123)
results_gam <- run_simu_gam_poly_rep(generate_data_4, eval_points, y_type = "binomial", n=n, rounds=500, return_all_rslts=T, method = "GAM")
save.image(file=here("data", "rdata", "02_simu_V5_sys4_1000_GAM.RData"))

rm(results_gam)
set.seed(123)
results_poly <- run_simu_gam_poly_rep(generate_data_4, eval_points, y_type = "binomial", n=n, rounds=500, return_all_rslts=T, method = "POLY")
save.image(file=here("data", "rdata", "02_simu_V5_sys4_1000_poly.RData"))
