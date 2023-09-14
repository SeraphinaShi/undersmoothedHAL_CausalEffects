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


# ----check_sys7----------------------------------------------------------------------------------------------------------------------
generate_data_7 <- function(n, a=NA){
  # exogenous variables
  U_W <- rnorm(n, 0, 1)
  U_A <- rnorm(n, 0, 1)
  U_Y <- runif(n, 0, 1)
  
  # endogenous variables
  W <- U_W
  
  if(is.na(a)){
    A <-  2.5 - 0.5*W + U_A
    A[A<=0] = 0
    A[A>=5] = 5
  } else {
    A <- rep(a, n)
  }
  
  Y <- as.numeric(U_Y < plogis(-6 + W + 3.5*A*as.numeric(A >= 2) - 4*A*as.numeric(A >= 4) - 0.5 * W * A ))
  
  # data frame
  O <- data.frame(W, A, Y)
  return(O)
}

# df = generate_data_7(10000)
# hist(df$A)
# hist(df$W)
# mean(df$Y)

## ----true_psi_sys7-------------------------------------------------------------------------------------------------------------------
# Getting trul value of psi
#------------------------------------------------------------------------------------

a_vec <- seq(0,5,0.1)
psi0_a <- c()

for (i in 1:length(a_vec)) {
  a <- a_vec[i]
  if(a < 2 | a >= 4){
    psi0_a[i] = 0
  } else {
    EW = 0
    psi0_a[i] = plogis(-6 + EW + 3.5*a*as.numeric(a >= 2) - 4*a*as.numeric(a >= 4) - 0.5 * EW * a )
  }
  
}

psi0_line <- data.frame(a=a_vec, psi0 = psi0_a)

eval_points = seq(0, 5, 0.5)
psi0_pnt <- psi0_line[psi0_line$a %in% eval_points,]

# ## -----------------------------------------------------------------------------------------------------------------------
# load(here("data", "rdata", "02_simu_V5_sys7_psi0.RData"))

source(here("scripts", "scripts_v5_final", "1_hal_functions.R"))
source(here("scripts", "scripts_v5_final", "1_simu_functions.R"))

n = 1000

# -----------------------------------------------------------------------------------------------------------------------
# 
# set.seed(123)
# #
# results_0 <- run_simu_rep(generate_data_7, eval_points, y_type = "binomial", n=n, rounds=500, return_all_rslts=T, defualt_setting = T)
# save.image(file=here("data", "rdata", "02_simu_v5_sys7_1000_default.RData"))
# p <- plot_performences_cv_ug_alla_noBT(df = results$result_summary)

# -----------------------------------------------------------------------------------------------------------------------
# rm(results_0)

set.seed(123)
results_adapt <- run_simu_smoothness_adaptive_HAL_rep(generate_data_7, eval_points, y_type = "binomial", n=n, rounds=500, return_all_rslts=T)
save.image(file=here("data", "rdata", "02_simu_v5_sys7_1000_adapt.RData"))
#p <- plot_performences_adapt(df = results_adapt$result_summary)

# -----------------------------------------------------------------------------------------------------------------------
rm(results_adapt)

set.seed(123)

results <- run_simu_rep(generate_data_7, eval_points, y_type = "binomial", n=n, rounds=500, return_all_rslts=T)
save.image(file=here("data", "rdata", "02_simu_V5_sys7_1000.RData"))


# ## -----------------------------------------------------------------------------------------------------------------------
rm(results)

set.seed(123)
results_grid <- run_simu_scaled_rep(generate_data_7, eval_points, y_type = "binomial", n=n, rounds=500, return_all_rslts=T)

save.image(file=here("data", "rdata", "02_simu_V5_sys7_1-00_grid.RData"))














