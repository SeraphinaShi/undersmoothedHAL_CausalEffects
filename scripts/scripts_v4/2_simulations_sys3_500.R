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



## ----check_sys3----------------------------------------------------------------------------------------------------------------------

generate_data_3 <- function(n, a=NA, z=NA){
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
  
  Y <- as.numeric(U_Y < plogis(-10 - 3*W + 2*A + Z * (5 + 2*sin(A^2) -20*as.numeric(A > 4)) ))
  
  # data frame
  O <- data.frame(W, A, Z, Y)
  return(O)
}

## ----true_psi_sys3-------------------------------------------------------------------------------------------------------------------
# Getting trul value of psi
#------------------------------------------------------------------------------------
# 
# a_vec <- seq(0,5,0.1)
# psi0_a_0 <- c()
# psi0_a_1 <- c()
# 
# N = 1e+07
# data_0_0 <- generate_data_3(n=N, a=0, z=0)
# data_0_1 <- generate_data_3(n=N, a=0, z=1)
# 
# for (i in 1:length(a_vec)) {
#   a <- a_vec[i]
#   
#   data_a_0 <- generate_data_3(n=N, a=a, z=0)
#   psi0_a_0[i] <- mean(data_a_0$Y - data_0_0$Y)
#   
#   data_a_1 <- generate_data_3(n=N, a=a, z=1)
#   psi0_a_1[i] <- mean(data_a_1$Y - data_0_1$Y)
# }
# 
# psi0_pnt <- data.frame(a=rep(a_vec, 2), z=c(rep(1,length(a_vec)), rep(0,length(a_vec))), psi0 = c(psi0_a_1,psi0_a_0))
# psi0_10pnt <- psi0_pnt[psi0_pnt$a %in% seq(0.5,5,0.5),]
# 
# save.image(file=here("data", "rdata", "02_simu_V3_sys3_psi0.RData"))
#------------------------------------------------------------------------------------ 
load(file=here("data", "rdata", "02_simu_V3_sys3_psi0.RData"))
source(here("scripts", "scripts_v3", "1_simu_functions.R"))
source(here("scripts", "scripts_v3", "1_simu_functions_hal9001.R"))
#------------------------------------------------------------------------------------


## ----simu_sys3_n500-------------------------------------------------------------------------------------------------------------------
nn=500
# 
# ## ----simu_sys3_n500_1_cv, fig.width=6, fig.height=4-----------------------------------------------------------------------------------
# set.seed(123)
# results <- run_simu_1round(generate_data_3, n=nn)
# psi_10pnt <- merge(as.data.frame(psi0_10pnt), as.data.frame(results), by=c("a", "z"))
# 
# 
# ## ----simu_sys3_n500_B_cv, fig.width=6, fig.height=7-----------------------------------------------------------------------------------
# set.seed(123)
# simu_results <- run_simu_rep(generate_data_3, n=nn, B=1000, return_all_rslts=T)
# 
# save.image(file=here("data", "rdata", "02_simu_V3_sys3_500_CV.RData"))
# 
# ## ----simu_sys3_n500_1_u, fig.width=6, fig.height=4------------------------------------------------------------------------------------
# set.seed(123)
# n = nn
# results_under <- run_simu_1round(generate_data_3, n=nn, undersmooth=T)
# 
# psi_10pnt <- merge(as.data.frame(psi0_10pnt), as.data.frame(results_under), by=c("a", "z"))
# cat(paste0("Undersmoothed lambda: ", unique(psi_10pnt$lambda), "\n which is ", unique(psi_10pnt$lambda_scaler), " * lambda_CV"))
# 
# 
# ## ----simu_sys3_n500_B_u, fig.width=6, fig.height=7------------------------------------------------------------------------------------
# set.seed(123)
# simu_results <- run_simu_rep(generate_data_3, n=nn, B=1000, return_all_rslts=T,  undersmooth=T)
# 
# save.image(file=here("data", "rdata", "02_simu_V3_sys3_500_U.RData"))
# 
# 

## ----simu_sys3_n500_1_u_local, fig.width=6, fig.height=4------------------------------------------------------------------------------------
set.seed(123)
n = nn
results_under <- run_simu_1round(generate_data_3, n=nn, undersmooth='local')
psi_10pnt <- merge(as.data.frame(psi0_10pnt), as.data.frame(results_under), by=c("a", "z"))


## ----simu_sys3_n500_B_u_local, fig.width=6, fig.height=7------------------------------------------------------------------------------------
set.seed(123)
simu_results <- run_simu_rep(generate_data_3, n=nn, B=20, return_all_rslts=T,  undersmooth='local')
# [1] "round 1"
# [1] "  CV lambda: 0.001797"
# [1] "  globally u lambda: 0.000292"
# [1] "  locally u lambdas: 0.001797, 0.001797, 0.001797, 0.001797, 0.001797"
# [1] "round 2"
# [1] "  CV lambda: 0.000763"
# [1] "  globally u lambda: 0.000178"
# [1] "  locally u lambdas: 0.000369, 0.000369, 0.000531, 0.000763, 0.000369" （0.5, 1.5, 2.5, 3.5, 4.5）
# [1] "round 3"
# [1] "  CV lambda: 0.002268"
# [1] "  globally u lambda: 0.000368"
# [1] "  locally u lambdas: 0.001577, 0.001577, 0.001577, 0.002268, 0.002268"
# [1] "round 4"
# [1] "  CV lambda: 0.000856"
# [1] "  globally u lambda: 0.000097"
# [1] "  locally u lambdas: 0.000139, 0.000139, 0.000139, 0.000856, 0.000139"
# [1] "round 5"
# [1] "  CV lambda: 0.002401"
# [1] "  globally u lambda: 0.000271"
# [1] "  locally u lambdas: 0.001669, 0.001669, 0.001669, 0.002401, 0.002401"
# [1] "round 6"
# [1] "  CV lambda: 0.001694"
# [1] "  globally u lambda: 0.000569"
# [1] "  locally u lambdas: 0.001694, 0.001694, 0.001694, 0.001694, 0.001694"
# [1] "round 7"
# [1] "  CV lambda: 0.000592"
# [1] "  globally u lambda: 0.000096"
# [1] "  locally u lambdas: 0.000199, 0.000286, 0.000592, 0.000592, 0.000286"
# [1] "round 8"
# [1] "  CV lambda: 0.001835"
# [1] "  globally u lambda: 0.000144"
# [1] "  locally u lambdas: 0.001835, 0.001835, 0.001835, 0.001835, 0.001835"
# [1] "round 9"
# [1] "  CV lambda: 0.001392"
# [1] "  globally u lambda: 0.000325"
# [1] "  locally u lambdas: 0.001392, 0.001392, 0.001392, 0.001392, 0.001392"
# [1] "round 10"
# [1] "  CV lambda: 0.001487"
# [1] "  globally u lambda: 0.000347"
# [1] "  locally u lambdas: 0.001487, 0.001487, 0.001487, 0.001487, 0.001487"
# [1] "round 11"
# [1] "  CV lambda: 0.002186"
# [1] "  globally u lambda: 0.000510"
# [1] "  locally u lambdas: 0.002186, 0.002186, 0.002186, 0.002186, 0.002186"
# [1] "round 12"
# [1] "  CV lambda: 0.001252"
# [1] "  globally u lambda: 0.000203"
# [1] "  locally u lambdas: 0.001252, 0.001252, 0.001252, 0.001252, 0.001252"
# [1] "round 13"
# [1] "  CV lambda: 0.001503"
# [1] "  globally u lambda: 0.000351"
# [1] "  locally u lambdas: 0.001503, 0.001503, 0.001503, 0.001503, 0.001503"
# [1] "round 14"
# [1] "  CV lambda: 0.001981"
# [1] "  globally u lambda: 0.000155"
# [1] "  locally u lambdas: 0.001981, 0.001981, 0.001981, 0.001981, 0.001981"
# [1] "round 15"
# [1] "  CV lambda: 0.005359"
# [1] "  globally u lambda: 0.000421"
# [1] "  locally u lambdas: 0.0018, 0.0018, 0.00259, 0.00259, 0.00259"
# [1] "round 16"
# [1] "  CV lambda: 0.002337"
# [1] "  globally u lambda: 0.000264"
# [1] "  locally u lambdas: 0.002337, 0.002337, 0.002337, 0.002337, 0.002337"
# [1] "round 17"
# [1] "  CV lambda: 0.000290"
# [1] "  globally u lambda: 0.000016"
# [1] "  locally u lambdas: 0.00029, 0.00029, 0.00029, 0.00029, 0.00029"
# [1] "round 18"
# [1] "  CV lambda: 0.003251"
# [1] "  globally u lambda: 0.000367"
# [1] "  locally u lambdas: 0.003251, 0.003251, 0.003251, 0.003251, 0.003251"
# [1] "round 19"
# [1] "  CV lambda: 0.000987"
# [1] "  globally u lambda: 0.000332"
# [1] "  locally u lambdas: 0.000987, 0.000987, 0.000987, 0.000987, 0.000987"
# [1] "round 20"
# [1] "  CV lambda: 0.002674"
# [1] "  globally u lambda: 0.000625"
# [1] "  locally u lambdas: 0.002674, 0.002674, 0.002674, 0.002674, 0.002674"
save.image(file=here("data", "rdata", "02_simu_V3_sys3_500_U_l.RData"))

# 
# ## ----simu_sys3_n500_B_grid------------------------------------------------------------------------------------------------------------
# set.seed(123)
# 
# lambda_scalers <- c(1.2, 1.1, 10^seq(from=0, to=-3, length=30))
# 
# simu_results_lists <- list()
# for(i in 1:length(lambda_scalers)){
#   scaler = lambda_scalers[i]
#   simu_results_lists[[i]] <- run_simu_rep(generate_data_3, n=nn, B=1000, lambda_scaler=scaler, return_all_rslts=F,  undersmooth=F)
# }
# simu_results_all <- do.call("rbind", simu_results_lists) %>% as.data.frame()
# 
# save.image(file=here("data", "rdata", "02_simu_V3_sys3_500_grid.RData"))
