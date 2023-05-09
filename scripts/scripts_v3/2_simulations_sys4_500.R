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

generate_data_4 <- function(n, a=NA, z=NA){
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
  
  Y <- as.numeric(U_Y < plogis(-10 - 3*W + 2*A + Z * (5 + 2*sin(0.5*A^2)*as.numeric(A > 2)) ))
  
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
data_0_0 <- generate_data_4(n=N, a=0, z=0)
data_0_1 <- generate_data_4(n=N, a=0, z=1)

for (i in 1:length(a_vec)) {
  a <- a_vec[i]

  data_a_0 <- generate_data_4(n=N, a=a, z=0)
  psi0_a_0[i] <- mean(data_a_0$Y) # - data_0_0$Y

  data_a_1 <- generate_data_4(n=N, a=a, z=1)
  psi0_a_1[i] <- mean(data_a_1$Y) #  - data_0_1$Y
}

psi0_pnt <- data.frame(a=rep(a_vec, 2), z=c(rep(1,length(a_vec)), rep(0,length(a_vec))), psi0 = c(psi0_a_1,psi0_a_0))
psi0_10pnt <- psi0_pnt[psi0_pnt$a %in% seq(0.5,5,0.5),]

ggplot() +
  geom_line(data=psi0_pnt %>% mutate(z = as.factor(z)), aes(x=a, y=psi0, color=z, group=z)) + 
  geom_point(data=psi0_10pnt %>% mutate(z = as.factor(z)), aes(x=a, y=psi0, color=z, group=z)) + 
  labs(x="a", y="ATE",
       title = "True Average Treatment Effect \n  P_0(E[Y|a,z])") +
  theme(plot.title = element_text(hjust = 0.5), #size=8.6, 
        #axis.title.y = element_blank(),
        panel.grid.major.x = element_blank(),
        #axis.text.y = element_text(size=8),
        #axis.title.x = element_text(size=8),
        axis.text = element_text(size=7)) 
# 
