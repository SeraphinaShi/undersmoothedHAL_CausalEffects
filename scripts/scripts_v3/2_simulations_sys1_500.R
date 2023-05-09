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
plotFolder <- here("results","images", "v3")
if(!file.exists(plotFolder)) dir.create(plotFolder,recursive=TRUE)

knitr::opts_chunk$set(
  cache=FALSE, autodep=FALSE, warning=FALSE, message=FALSE, echo=FALSE,
  results = 'markup', dev='png', dpi=150, fig.align = "center", fig.path=paste0(plotFolder, "/"),
  cache.path=".cache/",
  duplicate.label="allow"
)

source(here("scripts", "scripts_v3", "1_simu_functions_hal9001.R"))
source(here("scripts", "scripts_v3", "1_simu_functions.R"))


## ----check_sys1----------------------------------------------------------------------------------------------------------------------
generate_data_1 <- function(n, a=NA, z=NA){
  # exogenous variables
  U_W <- rnorm(n, 0, 1)
  U_A <- rnorm(n, 0, 2)
  U_Z <- runif(n, 0, 1)
  U_Y <- runif(n, 0, 1)
  
  # endogenous variables
  W <- U_W
  
  if(is.na(a)){
    A <-  2 - 0.5*W + U_A
    A[A<=0] = 0
    A[A>=5] = 5
  } else {
    A <- rep(a, n)
  }
  
  if(is.na(z)){
    Z <- as.numeric(U_Z < plogis(2-W-A))
  } else {
    Z <- rep(z, n)
  }
  
  Y <- as.numeric(U_Y < plogis(W + 5*A + Z - 0.5 * W * A - 8))
  
  # data frame
  O <- data.frame(W, A, Z, Y)
  return(O)
}

## ----true_psi_sys1-------------------------------------------------------------------------------------------------------------------
# Getting trul value of psi
#------------------------------------------------------------------------------------
# a_vec <- seq(0,5,0.1)
# psi0_a_0 <- c()
# psi0_a_1 <- c()
# 
# N = 1e+07
# data_0_0 <- generate_data_1(n=N, a=0, z=0)
# data_0_1 <- generate_data_1(n=N, a=0, z=1)
# 
# for (i in 1:length(a_vec)) {
#   a <- a_vec[i]
# 
#   data_a_0 <- generate_data_1(n=N, a=a, z=0)
#   psi0_a_0[i] <- mean(data_a_0$Y - data_0_0$Y)
# 
#   data_a_1 <- generate_data_1(n=N, a=a, z=1)
#   psi0_a_1[i] <- mean(data_a_1$Y - data_0_1$Y)
# }
# 
# psi0_pnt <- data.frame(a=rep(a_vec, 2), z=c(rep(1,length(a_vec)), rep(0,length(a_vec))), psi0 = c(psi0_a_1,psi0_a_0))
# psi0_10pnt <- psi0_pnt[psi0_pnt$a %in% seq(0.5,5,0.5),]
# save.image(file=here("data", "rdata", "02_simu_V3_sys1_psi0.RData"))
#------------------------------------------------------------------------------------
load(file=here("data", "rdata", "02_simu_V3_sys1_psi0.RData"))
source(here("scripts", "scripts_v3", "1_simu_functions_hal9001.R"))
source(here("scripts", "scripts_v3", "1_simu_functions.R"))
#------------------------------------------------------------------------------------


# 
# ## ----simu_sys1_n500-------------------------------------------------------------------------------------------------------------------
nn=500
# 
# 
# 
# ## ----simu_sys1_n500_1_cv, fig.width=6, fig.height=4-----------------------------------------------------------------------------------
# set.seed(123)
# results <- run_simu_1round(generate_data_1, n=nn)
# psi_10pnt <- merge(as.data.frame(psi0_10pnt), as.data.frame(results), by=c("a", "z"))
# 
# 
# ## ----simu_sys1_n500_B_cv, fig.width=6, fig.height=7-----------------------------------------------------------------------------------
# set.seed(123)
# simu_results <- run_simu_rep(generate_data_1, n=nn, B=1000, return_all_rslts=T)
# 
# save.image(file=here("data", "rdata", "02_simu_V3_sys1_500_CV.RData"))
# 
# 
# ## ----simu_sys1_n500_1_u, fig.width=6, fig.height=4------------------------------------------------------------------------------------
# set.seed(123)
# n = nn
# results_under <- run_simu_1round(generate_data_1, n=nn, undersmooth=T)
# psi_10pnt <- merge(as.data.frame(psi0_10pnt), as.data.frame(results_under), by=c("a", "z"))
# 
# 
# ## ----simu_sys1_n500_B_u, fig.width=6, fig.height=7------------------------------------------------------------------------------------
# set.seed(123)
# simu_results <- run_simu_rep(generate_data_1, n=nn, B=1000, return_all_rslts=T,  undersmooth=T)
# 
# save.image(file=here("data", "rdata", "02_simu_V3_sys1_500_U.RData"))
# 
# 

## ----simu_sys1_n500_1_u_local, fig.width=6, fig.height=4------------------------------------------------------------------------------------
set.seed(123)
n = nn
results_under <- run_simu_1round(generate_data_1, n=nn, undersmooth='local')
psi_10pnt <- merge(as.data.frame(psi0_10pnt), as.data.frame(results_under), by=c("a", "z"))


## ----simu_sys1_n500_B_u_local, fig.width=6, fig.height=7------------------------------------------------------------------------------------
set.seed(123)
simu_results <- run_simu_rep(generate_data_1, n=nn, B=20, return_all_rslts=T,  undersmooth='local')
# [1] "round 1"
# [1] "  CV lambda: 0.003049"
# [1] "  globally u lambda: 0.000239"
# [1] "  locally u lambdas: 0.003049, 0.003049, 0.003049, 0.003049, 0.003049"
# [1] "round 2"
# [1] "  CV lambda: 0.002055"
# [1] "  globally u lambda: 0.000026"
# [1] "  locally u lambdas: 0.002055, 0.002055, 0.002055, 0.002055, 0.002055"
# [1] "round 3"
# [1] "  CV lambda: 0.002323"
# [1] "  globally u lambda: 0.000377"
# [1] "  locally u lambdas: 0.002323, 0.002323, 0.002323, 0.002323, 0.002323"
# [1] "round 4"
# [1] "  CV lambda: 0.003342"
# [1] "  globally u lambda: 0.001615"
# [1] "  locally u lambdas: 0.003342, 0.003342, 0.003342, 0.003342, 0.003342"
# [1] "round 5"
# [1] "  CV lambda: 0.006453"
# [1] "  globally u lambda: 0.001048"
# [1] "  locally u lambdas: 0.006453, 0.006453, 0.006453, 0.006453, 0.006453"
# [1] "round 6"
# [1] "  CV lambda: 0.003260"
# [1] "  globally u lambda: 0.000761"
# [1] "  locally u lambdas: 0.00326, 0.00326, 0.00326, 0.00326, 0.00326"
# [1] "round 7"
# [1] "  CV lambda: 0.003335"
# [1] "  globally u lambda: 0.000541"
# [1] "  locally u lambdas: 0.003335, 0.003335, 0.003335, 0.003335, 0.003335"
# [1] "round 8"
# [1] "  CV lambda: 0.004956"
# [1] "  globally u lambda: 0.001665"
# [1] "  locally u lambdas: 0.004956, 0.004956, 0.004956, 0.004956, 0.004956"
# [1] "round 9"
# [1] "  CV lambda: 0.003956"
# [1] "  globally u lambda: 0.000447"
# [1] "  locally u lambdas: 0.003956, 0.003956, 0.003956, 0.003956, 0.003956"
# [1] "round 10"
# [1] "  CV lambda: 0.004587"
# [1] "  globally u lambda: 0.000250"
# [1] "  locally u lambdas: 0.004587, 0.004587, 0.004587, 0.004587, 0.004587"
# [1] "round 11"
# [1] "  CV lambda: 0.008346"
# [1] "  globally u lambda: 0.001949"
# [1] "  locally u lambdas: 0.008346, 0.008346, 0.008346, 0.008346, 0.008346"
# [1] "round 12"
# [1] "  CV lambda: 0.002642"
# [1] "  globally u lambda: 0.000888"
# [1] "  locally u lambdas: 0.002642, 0.002642, 0.002642, 0.002642, 0.002642"
# [1] "round 13"
# [1] "  CV lambda: 0.004751"
# [1] "  globally u lambda: 0.001596"
# [1] "  locally u lambdas: 0.004751, 0.004751, 0.004751, 0.004751, 0.004751"
# [1] "round 14"
# [1] "  CV lambda: 0.000404"
# [1] "  globally u lambda: 0.000022"
# [1] "  locally u lambdas: 0.000404, 0.000404, 0.000404, 0.000404, 0.000404"
# [1] "round 15"
# [1] "  CV lambda: 0.001618"
# [1] "  globally u lambda: 0.000378"
# [1] "  locally u lambdas: 0.001618, 0.001618, 0.001618, 0.001618, 0.001618"
# [1] "round 16"
# [1] "  CV lambda: 0.002955"
# [1] "  globally u lambda: 0.000690"
# [1] "  locally u lambdas: 0.002955, 0.002955, 0.002955, 0.002955, 0.002955"
# [1] "round 17"
# [1] "  CV lambda: 0.000885"
# [1] "  globally u lambda: 0.000207"
# [1] "  locally u lambdas: 0.000885, 0.000885, 0.000885, 0.000885, 0.000885"
# [1] "round 18"
# [1] "  CV lambda: 0.007975"
# [1] "  globally u lambda: 0.001295"
# [1] "  locally u lambdas: 0.007975, 0.007975, 0.007975, 0.007975, 0.007975"
# [1] "round 19"
# [1] "  CV lambda: 0.007726"
# [1] "  globally u lambda: 0.001805"
# [1] "  locally u lambdas: 0.007726, 0.007726, 0.007726, 0.007726, 0.007726"
# [1] "round 20"
# [1] "  CV lambda: 0.002753"
# [1] "  globally u lambda: 0.001331"
# [1] "  locally u lambdas: 0.002753, 0.002753, 0.002753, 0.002753, 0.002753"
# # save.image(file=here("data", "rdata", "02_simu_V3_sys1_500_U_l.RData"))
# 
# 
# ## ----simu_sys1_n500_B_grid------------------------------------------------------------------------------------------------------------
# set.seed(123)
# 
# lambda_scalers <- c(1.2, 1.1, 10^seq(from=0, to=-3, length=30))
# 
# simu_results_lists <- list()
# for(i in 1:length(lambda_scalers)){
#   scaler = lambda_scalers[i]
#   simu_results_lists[[i]] <- run_simu_rep(generate_data_1, n=nn, B=1000, lambda_scaler=scaler, return_all_rslts=F,  undersmooth=F)
# }
# simu_results_all <- do.call("rbind", simu_results_lists) %>% as.data.frame()
# 
# save.image(file=here("data", "rdata", "02_simu_V3_sys1_500_grid.RData"))
# 
