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



## ----check_sys6----------------------------------------------------------------------------------------------------------------------
expit <- function(x){
  return(exp(x)/(1+exp(x)))
}

generate_data_6 <- function(n, a=NA){
  W1 <- runif(n)
  W2 <- rbinom(n, size=1, prob =0.7)
  W3 <- rnorm(W1, 0.25*exp(2*W1))
  
  v_W <- exp(1 + 2*W1*expit(W3))
  mu_W <- expit(0.03 - 0.8*log(1+W2) + 0.9*exp(W1)*W2 - 0.4*atan(W3+2)*W2*W1)
  if(is.na(a)){
    A <- rbeta(n, v_W*mu_W, v_W*(1-mu_W))
  } else {
    A <- a
  }
  
  
  Q <- expit(-2 + 1.5*A + 5*A^3 - 2.5*W1 + 0.5*A*W2 - log(A)*W1*W2 + 0.5*(A^(3/4))*W1*W3)
  Y <- rbinom(n, size = 1, prob = Q)
  
  O <- data.frame(W1, W2,W3, A, Y)
  return(O)
}

## ----true_psi_sys6-------------------------------------------------------------------------------------------------------------------
# Getting trul value of psi
#------------------------------------------------------------------------------------

a_vec <- seq(0.05, 1, 0.01)
psi0_a <- c()
psi0_a <- c()

N = 1e+07

for (i in 1:length(a_vec)) {
  a <- a_vec[i]

  data_a <- generate_data_6(n=N, a=a)
  psi0_a[i] <- mean(data_a$Y) # - data_0$Y
}

psi0_line <- data.frame(a=a_vec, psi0 = psi0_a)
plot(psi0_line$a, psi0_line$psi0)

eval_points = seq(0.05,1,0.05)

psi0_pnt <- psi0_line[psi0_line$a %in% eval_points,]

save.image(file=here("data", "rdata", "02_simu_V5_sys6_psi0.RData"))



#======================================================================
# Simulations

## -----------------------------------------------------------------------------------------------------------------------
load(here("data", "rdata", "02_simu_V5_sys6_psi0.RData"))

source(here("scripts", "scripts_v5", "1_hal_functions.R"))
source(here("scripts", "scripts_v5", "1_simu_functions.R"))

#======================================================================
n = 200

## -----------------------------------------------------------------------------------------------------------------------
set.seed(123)

results <- run_simu_rep(generate_data_6, eval_points, y_type = "binomial", n=n, rounds=500, return_all_rslts=T)
save.image(file=here("data", "rdata", "02_simu_V5_sys6_200.RData"))

## -----------------------------------------------------------------------------------------------------------------------
rm(results)

set.seed(123)

results_0 <- run_simu_rep(generate_data_6, eval_points, y_type = "binomial", n=n, rounds=500, return_all_rslts=T, defualt_setting = T)
save.image(file=here("data", "rdata", "02_simu_V5_sys6_200_default.RData"))

## -----------------------------------------------------------------------------------------------------------------------
rm(results_0)

set.seed(123)
results_grid <- run_simu_scaled_rep(generate_data_6, eval_points, y_type = "binomial", n=n, rounds=500, return_all_rslts=T)

save.image(file=here("data", "rdata", "02_simu_V5_sys6_200_grid.RData"))

## -----------------------------------------------------------------------------------------------------------------------
rm(results_grid)

set.seed(123)
results_adapt <- run_simu_smoothness_adaptive_HAL_rep(generate_data_6, eval_points, y_type = "binomial", n=n, rounds=500, return_all_rslts=T)
save.image(file=here("data", "rdata", "02_simu_V5_sys6_200_adapt.RData"))



## -----------------------------------------------------------------------------------------------------------------------
rm(results_adapt)

source(here("scripts", "scripts_v5", "1_simu_functions_noHAL.R"))

set.seed(123)
results_gam <- run_simu_gam_poly_rep(generate_data_6, eval_points, y_type = "binomial", n=n, rounds=500, return_all_rslts=T, method = "GAM")
save.image(file=here("data", "rdata", "02_simu_V5_sys6_200_GAM.RData"))

rm(results_gam)
set.seed(123)
results_poly <- run_simu_gam_poly_rep(generate_data_6, eval_points, y_type = "binomial", n=n, rounds=500, return_all_rslts=T, method = "POLY")
save.image(file=here("data", "rdata", "02_simu_V5_sys6_200_poly.RData"))





#======================================================================
rm(results_poly)

n = 500

## -----------------------------------------------------------------------------------------------------------------------
set.seed(123)

results <- run_simu_rep(generate_data_6, eval_points, y_type = "binomial", n=n, rounds=500, return_all_rslts=T)
save.image(file=here("data", "rdata", "02_simu_V5_sys6_500.RData"))

## -----------------------------------------------------------------------------------------------------------------------
rm(results)

set.seed(123)

results_0 <- run_simu_rep(generate_data_6, eval_points, y_type = "binomial", n=n, rounds=500, return_all_rslts=T, defualt_setting = T)
save.image(file=here("data", "rdata", "02_simu_V5_sys6_500_default.RData"))

## -----------------------------------------------------------------------------------------------------------------------
rm(results_0)

set.seed(123)
results_grid <- run_simu_scaled_rep(generate_data_6, eval_points, y_type = "binomial", n=n, rounds=500, return_all_rslts=T)

save.image(file=here("data", "rdata", "02_simu_V5_sys6_500_grid.RData"))

## -----------------------------------------------------------------------------------------------------------------------
rm(results_grid)

set.seed(123)
results_adapt <- run_simu_smoothness_adaptive_HAL_rep(generate_data_6, eval_points, y_type = "binomial", n=n, rounds=500, return_all_rslts=T)
save.image(file=here("data", "rdata", "02_simu_V5_sys6_500_adapt.RData"))



## -----------------------------------------------------------------------------------------------------------------------
rm(results_adapt)

source(here("scripts", "scripts_v5", "1_simu_functions_noHAL.R"))

set.seed(123)
results_gam <- run_simu_gam_poly_rep(generate_data_6, eval_points, y_type = "binomial", n=n, rounds=500, return_all_rslts=T, method = "GAM")
save.image(file=here("data", "rdata", "02_simu_V5_sys6_500_GAM.RData"))

rm(results_gam)
set.seed(123)
results_poly <- run_simu_gam_poly_rep(generate_data_6, eval_points, y_type = "binomial", n=n, rounds=500, return_all_rslts=T, method = "POLY")
save.image(file=here("data", "rdata", "02_simu_V5_sys6_500_poly.RData"))





#======================================================================
rm(results_poly)

n = 1000

## -----------------------------------------------------------------------------------------------------------------------
set.seed(123)

results <- run_simu_rep(generate_data_6, eval_points, y_type = "binomial", n=n, rounds=500, return_all_rslts=T)
save.image(file=here("data", "rdata", "02_simu_V5_sys6_1000.RData"))

## -----------------------------------------------------------------------------------------------------------------------
rm(results)

set.seed(123)

results_0 <- run_simu_rep(generate_data_6, eval_points, y_type = "binomial", n=n, rounds=500, return_all_rslts=T, defualt_setting = T)
save.image(file=here("data", "rdata", "02_simu_V5_sys6_1000_default.RData"))

## -----------------------------------------------------------------------------------------------------------------------
rm(results_0)

set.seed(123)
results_grid <- run_simu_scaled_rep(generate_data_6, eval_points, y_type = "binomial", n=n, rounds=500, return_all_rslts=T)

save.image(file=here("data", "rdata", "02_simu_V5_sys6_1000_grid.RData"))

## -----------------------------------------------------------------------------------------------------------------------
rm(results_grid)

set.seed(123)
results_adapt <- run_simu_smoothness_adaptive_HAL_rep(generate_data_6, eval_points, y_type = "binomial", n=n, rounds=500, return_all_rslts=T)
save.image(file=here("data", "rdata", "02_simu_V5_sys6_1000_adapt.RData"))



## -----------------------------------------------------------------------------------------------------------------------
rm(results_adapt)

source(here("scripts", "scripts_v5", "1_simu_functions_noHAL.R"))

set.seed(123)
results_gam <- run_simu_gam_poly_rep(generate_data_6, eval_points, y_type = "binomial", n=n, rounds=500, return_all_rslts=T, method = "GAM")
save.image(file=here("data", "rdata", "02_simu_V5_sys6_1000_GAM.RData"))

rm(results_gam)
set.seed(123)
results_poly <- run_simu_gam_poly_rep(generate_data_6, eval_points, y_type = "binomial", n=n, rounds=500, return_all_rslts=T, method = "POLY")
save.image(file=here("data", "rdata", "02_simu_V5_sys6_1000_poly.RData"))
