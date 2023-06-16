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



## ----check_sys6----------------------------------------------------------------------------------------------------------------------

generate_data_6 <- function(n, a=NA){
  # exogenous variables

  # data frame
  O <- W
  O$A = A
  O$Y = Y
  return(O)
}

## ----true_psi_sys6-------------------------------------------------------------------------------------------------------------------
# Getting trul value of psi
#------------------------------------------------------------------------------------

a_vec <- seq(0,5,0.1)
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

eval_points = seq(0,5,0.5)
psi0_pnt <- psi0_line[psi0_line$a %in% eval_points,]

save.image(file=here("data", "rdata", "02_simu_V5_sys6_psi0.RData"))



## ----simu_sys1_n200-------------------------------------------------------------------------------------------------------------------
n = 200

set.seed(123)

results <- run_simu_rep(generate_data_6, n=n, rounds=500, return_all_rslts=T)
save.image(file=here("data", "rdata", "02_simu_V5_sys6_200.RData"))


## ----simu_sys2_n200-------------------------------------------------------------------------------------------------------------------
rm(results)

set.seed(123)
results_grid <- run_simu_scaled_rep(generate_data_6, n=n, rounds=500, return_all_rslts=T)

save.image(file=here("data", "rdata", "02_simu_V5_sys6_200_grid.RData"))



