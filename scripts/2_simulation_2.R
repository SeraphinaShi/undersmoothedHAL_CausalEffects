## ----load_lib, include = FALSE, warning=FALSE, message=FALSE, echo=FALSE-----------------------------------------
library(here)
library(data.table)
library(dplyr)
library(foreach)

library(stringr)
library(glmnet)

library(origami)
library(hal9001)
library(tictoc)

library(R.utils)


## ----setup, include = FALSE--------------------------------------------------------------------------------------
plotFolder <- here("results","images")
if(!file.exists(plotFolder)) dir.create(plotFolder,recursive=TRUE)

knitr::opts_chunk$set(
  cache=FALSE, autodep=FALSE, warning=FALSE, message=FALSE, echo=FALSE,
  results = 'markup', dev='png', dpi=150, fig.align = "center", fig.path=plotFolder,
  cache.path=".cache/",
  duplicate.label="allow"
)

source(here("scripts", "1_simu_functions_hal9001.R"))
source(here("scripts", "1_simu_functions.R"))


## ----------------------------------------------------------------------------------------------------------------
generate_data2 <- function(n, a=NA, z=NA){
  # exogenous variables
  U_W <- runif(n, 0, 1)
  U_A <- rnorm(n,0,4)
  U_Z <- runif(n, 0, 1)
  U_Y <- runif(n, 0, 1)
  
  
  # endogenous variables
  W <- U_W
  
  if(is.na(a)){
    A <-  10 - 10*W  + U_A
    A[A<=0] = 0
    A[A>=10] = 10
  } else {
    A <- rep(a, n)
  }
  
  if(is.na(z)){
    Z <- as.numeric(U_Z < plogis(0.2*(W+A)))
  } else {
    Z <- rep(z, n)
  }
  
  Y <- as.numeric(U_Z < plogis(0.05 * (-W + 5*A - 3 * W * A - Z)))
  

  # data frame
  O <- data.frame(W, A, Z, Y)
  return(O)
}


obs <- generate_data2(n=10000)
summary(obs)

# check positivity violations
cat("Summary of A given W < 0.25:")
summary(obs$A[obs$W < 0.25])
cat("Summary of A given 0.25 < W <= 0.5:")
summary(obs$A[0.25 <= obs$W & obs$W < 0.5])
cat("Summary of A given 0.5 < W <= 0.75:")
summary(obs$A[0.5 <= obs$W & obs$W < 0.75])
cat("Summary of A given 0.75 < W:")
summary(obs$A[0.75 <= obs$W])

par(mfrow=c(2,3))
plot(obs$W,obs$A)
plot(obs$A,obs$Z)
plot(obs$Z,obs$Y)
plot(obs$W,obs$Z)
plot(obs$W,obs$Y)
plot(obs$A,obs$Y)



## ----get_true_psis-----------------------------------------------------------------------------------------------
# Getting trul value of psi
a_vec <- seq(0.5,10,0.5)
psi0_a_0 <- c()
psi0_a_1 <- c()

N = 1e+07
data_0_0 <- generate_data2(n=N, a=0, z=0)
data_0_1 <- generate_data2(n=N, a=0, z=1)

for (i in 1:length(a_vec)) {
  a <- a_vec[i]

  data_a_0 <- generate_data2(n=N, a=a, z=0)
  psi0_a_0[i] <- mean(data_a_0$Y - data_0_0$Y)

  data_a_1 <- generate_data2(n=N, a=a, z=1)
  psi0_a_1[i] <- mean(data_a_1$Y - data_0_1$Y)
}


## ----simu_n500---------------------------------------------------------------------------------------------------
set.seed(123)
n=500
B=1000
results_500 = run_simu(generate_data = generate_data2, n=n, B=B)
save.image(file=here("data", "rdata", "02_simulation_2_500.RData"))


## ----------------------------------------------------------------------------------------------------------------
load(here("data", "rdata", "02_simulation_2_500.RData"))

results_500[[1]] <- results_500[[1]] %>% mutate_if(is.numeric, round, digits=4)
results_500[[2]] <- results_500[[2]] %>% mutate_if(is.numeric, round, digits=4)

print(results_500)


## ----simu_n1000--------------------------------------------------------------------------------------------------
set.seed(234)
n=1000
B=1000
results_1000 = run_simu(generate_data = generate_data2, n=n, B=B)
save.image(file=here("data", "rdata", "02_simulation_2_1000.RData"))


## ----------------------------------------------------------------------------------------------------------------
load(here("data", "rdata", "02_simulation_2_1000.RData"))

results_1000[[1]] <- results_1000[[1]] %>% mutate_if(is.numeric, round, digits=4)
results_1000[[2]] <- results_1000[[2]] %>% mutate_if(is.numeric, round, digits=4)

print(results_1000)


## ----simu_n5000--------------------------------------------------------------------------------------------------
set.seed(345)
n=5000
B=1000
results_5000 = run_simu(generate_data = generate_data2, n=n, B=B)
save.image(file=here("data", "rdata", "02_simulation_2_5000.RData"))


## ----------------------------------------------------------------------------------------------------------------
load(here("data", "rdata", "02_simulation_2_5000.RData"))

results_5000[[1]] <- results_5000[[1]] %>% mutate_if(is.numeric, round, digits=4)
results_5000[[2]] <- results_5000[[2]] %>% mutate_if(is.numeric, round, digits=4)

print(results_5000)


## ----------------------------------------------------------------------------------------------------------------
knitr::purl(here("scripts", "2_simulation_2.Rmd"))

