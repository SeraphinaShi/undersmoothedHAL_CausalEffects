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
  
  Y <- as.numeric(U_Y < plogis(-10 - 3*W + 2*A + Z * (5 + 2*sin(A^2) -20*as.numeric(a > 4)) ))
  
  # data frame
  O <- data.frame(W, A, Z, Y)
  return(O)
}

## ----true_psi_sys3-------------------------------------------------------------------------------------------------------------------
# Getting trul value of psi
#------------------------------------------------------------------------------------

a_vec <- seq(0,5,0.1)
psi0_a_0 <- c()
psi0_a_1 <- c()

N = 1e+07
data_0_0 <- generate_data_3(n=N, a=0, z=0)
data_0_1 <- generate_data_3(n=N, a=0, z=1)

for (i in 1:length(a_vec)) {
  a <- a_vec[i]
  
  data_a_0 <- generate_data_3(n=N, a=a, z=0)
  psi0_a_0[i] <- mean(data_a_0$Y - data_0_0$Y)
  
  data_a_1 <- generate_data_3(n=N, a=a, z=1)
  psi0_a_1[i] <- mean(data_a_1$Y - data_0_1$Y)
}

psi0_pnt <- data.frame(a=rep(a_vec, 2), z=c(rep(1,length(a_vec)), rep(0,length(a_vec))), psi0 = c(psi0_a_1,psi0_a_0))
psi0_10pnt <- psi0_pnt[psi0_pnt$a %in% seq(0.5,5,0.5),]

save.image(file=here("data", "rdata", "02_simu_V3_sys3_psi0.RData"))
#------------------------------------------------------------------------------------ 
# load(file=here("data", "rdata", "02_simu_V3_sys3_psi0.RData"))
#------------------------------------------------------------------------------------


## ----simu_sys3_n500-------------------------------------------------------------------------------------------------------------------
nn=500

## ----simu_sys3_n500_1_cv, fig.width=6, fig.height=4-----------------------------------------------------------------------------------
set.seed(123)
results <- run_simu_1round(generate_data_3, n=nn)
psi_10pnt <- merge(as.data.frame(psi0_10pnt), as.data.frame(results), by=c("a", "z"))


## ----simu_sys3_n500_B_cv, fig.width=6, fig.height=7-----------------------------------------------------------------------------------
set.seed(123)
simu_results <- run_simu_rep(generate_data_3, n=nn, B=1000, return_all_rslts=T)

save.image(file=here("data", "rdata", "02_simu_V3_sys3_500_CV.RData"))

## ----simu_sys3_n500_1_u, fig.width=6, fig.height=4------------------------------------------------------------------------------------
set.seed(123)
n = nn
results_under <- run_simu_1round(generate_data_3, n=nn, undersmooth=T)

psi_10pnt <- merge(as.data.frame(psi0_10pnt), as.data.frame(results_under), by=c("a", "z"))
cat(paste0("Undersmoothed lambda: ", unique(psi_10pnt$lambda), "\n which is ", unique(psi_10pnt$lambda_scaler), " * lambda_CV"))


## ----simu_sys3_n500_B_u, fig.width=6, fig.height=7------------------------------------------------------------------------------------
set.seed(123)
simu_results <- run_simu_rep(generate_data_3, n=nn, B=1000, return_all_rslts=T,  undersmooth=T)

save.image(file=here("data", "rdata", "02_simu_V3_sys3_500_U.RData"))



## ----simu_sys3_n500_B_grid------------------------------------------------------------------------------------------------------------
set.seed(123)

lambda_scalers <- c(1.2, 1.1, 10^seq(from=0, to=-3, length=30))

simu_results_lists <- list()
for(i in 1:length(lambda_scalers)){
  scaler = lambda_scalers[i]
  simu_results_lists[[i]] <- run_simu_rep(generate_data_3, n=nn, B=1000, lambda_scaler=scaler, return_all_rslts=F,  undersmooth=F)
}
simu_results_all <- do.call("rbind", simu_results_lists) %>% as.data.frame()

save.image(file=here("data", "rdata", "02_simu_V3_sys3_500_grid.RData"))
