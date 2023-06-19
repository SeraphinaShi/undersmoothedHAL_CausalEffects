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



## ----check_sys5----------------------------------------------------------------------------------------------------------------------

generate_data_5 <- function(n, a=NA){
  # exogenous variables
  U_W_var <- matrix(c(1,        0.8,      0.9*1.5,     0.2*0.5,      0.1*0.3,
                      0.8,      1,        0.8*1.5,     0.3*0.5,      0.1*0.3, 
                      0.9*1.5,  0.8*1.5,  1.5^2,       0.1*0.5*1.5,  0.5*0.3*1.5,
                      0.2*0.5,  0.3*0.5,  0.1*0.5*1.5, 0.5^2,        0.75*0.3*0.5,
                      0.1*0.3,  0.1*0.3,  0.5*0.3*1.5, 0.75*0.3*0.5, 0.3^2), nrow = 5)
  U_W <- rmvnorm(n, mean = c(0, 0, 5, 1, 1), sigma = U_W_var)
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

  Y <- as.numeric(U_Y < plogis(-10 - W$W1 + 2*W$W2 - 0.5*W$W1 * W$W3 - 0.5*W$W4 + 4*A + 5*sin((0.8*A)^2-2.6)*as.numeric(A > 2)) )
  mean(Y)
  
  # data frame
  O <- W
  O$A = A
  O$Y = Y
  return(O)
}

## ----true_psi_sys5-------------------------------------------------------------------------------------------------------------------
# Getting trul value of psi
#------------------------------------------------------------------------------------

# a_vec <- seq(0,5,0.1)
# psi0_a <- c()
# psi0_a <- c()
# 
# N = 1e+07
# 
# for (i in 1:length(a_vec)) {
#   a <- a_vec[i]
# 
#   data_a <- generate_data_5(n=N, a=a)
#   psi0_a[i] <- mean(data_a$Y) # - data_0$Y
# }
# 
# psi0_line <- data.frame(a=a_vec, psi0 = psi0_a)
# # plot(psi0_line$a, psi0_line$psi0)
# 
# eval_points = seq(0,5,0.5)
# psi0_pnt <- psi0_line[psi0_line$a %in% eval_points,]
# 
# save.image(file=here("data", "rdata", "02_simu_V5_sys5_psi0.RData"))



## -----------------------------------------------------------------------------------------------------------------------
load(here("data", "rdata", "02_simu_V5_sys5_psi0.RData"))
source(here("scripts", "scripts_v5", "1_hal_functions.R"))
source(here("scripts", "scripts_v5", "1_simu_functions.R"))


n = 500

set.seed(123)

results <- run_simu_rep(generate_data_5, n=n, rounds=500, return_all_rslts=T)
save.image(file=here("data", "rdata", "02_simu_V5_sys5_500.RData"))


## ----simu_sys2_n200-------------------------------------------------------------------------------------------------------------------
rm(results)

set.seed(123)
results_grid <- run_simu_scaled_rep(generate_data_5, n=n, rounds=500, return_all_rslts=T)

save.image(file=here("data", "rdata", "02_simu_V5_sys5_500_grid.RData"))



