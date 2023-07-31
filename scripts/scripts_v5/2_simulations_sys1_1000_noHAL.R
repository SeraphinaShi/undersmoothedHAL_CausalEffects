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


# ----check_sys1----------------------------------------------------------------------------------------------------------------------
generate_data_1 <- function(n, a=NA){
  # exogenous variables
  U_W <- rnorm(n, 0, 1)
  U_A <- rnorm(n, 0, 2)
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
  
  
  Y <- as.numeric(U_Y < plogis(-5 + W + 2.25*A - 0.5 * W * A ))
  
  # data frame
  O <- data.frame(W, A, Y)
  return(O)
}

## ----true_psi_sys1-------------------------------------------------------------------------------------------------------------------
# Getting trul value of psi
#------------------------------------------------------------------------------------

# a_vec <- seq(0,5,0.1)
# psi0_a <- c()
# psi0_a <- c()
# 
# N = 1e+07
# data_0 <- generate_data_1(n=N, a=0)
# 
# for (i in 1:length(a_vec)) {
#   a <- a_vec[i]
# 
#   data_a <- generate_data_1(n=N, a=a)
#   psi0_a[i] <- mean(data_a$Y) # - data_0$Y
# }
# 
# psi0_line <- data.frame(a=a_vec, psi0 = psi0_a)
# 
# # plot(x= psi0_line$a, psi0_line$psi0)
# 
# eval_points = seq(0, 5, 0.5)
# psi0_pnt <- psi0_line[psi0_line$a %in% eval_points,]
# 
# save.image(file=here("data", "rdata", "02_simu_v5_sys1_psi0.RData"))


## -----------------------------------------------------------------------------------------------------------------------
load(here("data", "rdata", "02_simu_v5_sys1_psi0.RData"))

source(here("scripts", "scripts_v5", "1_hal_functions.R"))
source(here("scripts", "scripts_v5", "1_simu_functions.R"))

n = 1000

## -----------------------------------------------------------------------------------------------------------------------
set.seed(123)
results_so <- run_simu_smooth_orders_rep(generate_data_1, eval_points, y_type = "binomial", n=n, rounds=500, return_all_rslts=T)

save.image(file=here("data", "rdata", "02_simu_V5_sys1_1000_SO.RData"))
## -----------------------------------------------------------------------------------------------------------------------
# 
# source(here("scripts", "scripts_v5", "1_simu_functions_noHAL.R"))
# 
# set.seed(123)
# results_gam <- run_simu_gam_poly_rep(generate_data_1, eval_points, y_type = "binomial", n=n, rounds=500, return_all_rslts=T, method = "GAM")
# save.image(file=here("data", "rdata", "02_simu_V5_sys1_1000_GAM.RData"))
# 
# rm(results_gam)
# set.seed(123)
# results_poly <- run_simu_gam_poly_rep(generate_data_1, eval_points, y_type = "binomial", n=n, rounds=500, return_all_rslts=T, method = "POLY")
# save.image(file=here("data", "rdata", "02_simu_V5_sys1_1000_poly.RData"))




