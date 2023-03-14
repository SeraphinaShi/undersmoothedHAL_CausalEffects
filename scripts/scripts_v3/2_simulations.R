## ----load_lib, include = FALSE, warning=FALSE, message=FALSE, echo=FALSE----------------------------------------------------
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
library(gridExtra)
library(grid)


## ----setup, include = FALSE-------------------------------------------------------------------------------------------------
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


## ---------------------------------------------------------------------------------------------------------------------------
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

obs <- generate_data_1(n=10000)
print(summary(obs))
  
# check positivity violations
cat("Summary of A given W < -1:")
summary(obs$A[obs$W < -1])
cat("Summary of A given -1 < W <= 0:")
summary(obs$A[-1 <= obs$W & obs$W < 0])
cat("Summary of A given 0 < W <= 1:")
summary(obs$A[0 <= obs$W & obs$W < 1])
cat("Summary of A given 1 < W:")
summary(obs$A[1 <= obs$W])

par(mfrow=c(2,3))
plot(obs$W,obs$A)
plot(obs$A,obs$Z)
plot(obs$Z,obs$Y)
plot(obs$W,obs$Z)
plot(obs$W,obs$Y)
plot(obs$A,obs$Y)

glm_fit <- glm(formula = Y ~ W + A + W*A + Z,
               family = binomial,
               data = obs)
summary(glm_fit)
y_preds <- predict(glm_fit, type = "response")
mse <- sum((y_preds - obs$Y)^2)
auc <- auc(obs$Y, y_preds)
print(paste0("    MSE: ", round(mse, 4), ", AUC: ", round(auc, 4)))


## ----get_true_psis_1--------------------------------------------------------------------------------------------------------
# Getting trul value of psi
a_vec <- seq(0,5,0.1)
psi0_a_0 <- c()
psi0_a_1 <- c()

N = 1e+07
data_0_0 <- generate_data_1(n=N, a=0, z=0)
data_0_1 <- generate_data_1(n=N, a=0, z=1)

for (i in 1:length(a_vec)) {
  a <- a_vec[i]

  data_a_0 <- generate_data_1(n=N, a=a, z=0)
  psi0_a_0[i] <- mean(data_a_0$Y - data_0_0$Y)

  data_a_1 <- generate_data_1(n=N, a=a, z=1)
  psi0_a_1[i] <- mean(data_a_1$Y - data_0_1$Y)
}

psi0_pnt <- data.frame(a=rep(a_vec, 2), z=c(rep(1,length(a_vec)), rep(0,length(a_vec))), psi0 = c(psi0_a_1,psi0_a_0))
psi0_10pnt <- psi0_pnt[psi0_pnt$a %in% seq(0.5,5,0.5),]

p <- ggplot() +
    geom_line(data=psi0_pnt %>% mutate(z = as.factor(z)), aes(x=a, y=psi0, color=z, group=z)) + 
    geom_point(data=psi0_10pnt %>% mutate(z = as.factor(z)), aes(x=a, y=psi0, color=z, group=z)) + 
    labs(x="a", y="ATE",
         title = "True Average Treatment Effect \n  P_0(E[Y|a,z] - E[Y|0,z])") +
        theme(plot.title = element_text(hjust = 0.5), #size=8.6, 
              #axis.title.y = element_blank(),
              panel.grid.major.x = element_blank(),
              #axis.text.y = element_text(size=8),
              #axis.title.x = element_text(size=8),
              axis.text = element_text(size=7)) 

p

save.image(file=here("data", "rdata", "02_simulation_V3_1_200_psi0.RData"))


## ----simu_1_n200------------------------------------------------------------------------------------------------------------
set.seed(123)
nn=200

results_200 <- run_simu_1round(generate_data_1, n=nn)
psi_10pnt <- merge(as.data.frame(psi0_10pnt), as.data.frame(results_200), by=c("a", "z"))


## ----simu_1_n200_rslt, fig.width=6, fig.height=4----------------------------------------------------------------------------
p1 <- plot_perforences_1lambda_alla(psi_10pnt, z_para=1, est_plot_only = T) +
  labs(title = "z=1") + theme(plot.title = element_text(hjust = 0.5))
p0 <- plot_perforences_1lambda_alla(psi_10pnt, z_para=0, est_plot_only = T) +
  labs(title = "z=0") + theme(plot.title = element_text(hjust = 0.5))

p <- grid.arrange(p1, p0, nrow=1,
                  top = textGrob(paste0("Average Treatment Effect,  E[Y|a,z] - E[Y|0,z]",
                                         gp=gpar(fontface = 'bold')))) # gp=gpar(fontsize=11, fontface = 'bold')


## ----simu_B_n200_rslt, fig.width=6, fig.height=7----------------------------------------------------------------------------
set.seed(123)
simu_results <- run_simu_rep(generate_data_1, n=nn, B=1000, return_all_rslts=T)
cat("z=1:")
p1_list <- plot_perforences_1lambda_alla(simu_results$result_summary, z_para=1, plot_list=F)
cat("z=0:")
p0_list <- plot_perforences_1lambda_alla(simu_results$result_summary, z_para=0, plot_list=F)

save.image(file=here("data", "rdata", "02_simulation_V3_1_200_CV.RData"))


## ----simu_1_n200_u----------------------------------------------------------------------------------------------------------
set.seed(123)
n=200

results_200_under <- run_simu_1round(generate_data_1, n=nn, undersmooth=T)
psi_10pnt <- merge(as.data.frame(psi0_10pnt), as.data.frame(results_200_under), by=c("a", "z"))


## ----simu_1_n200_rslt_u, fig.width=6, fig.height=4--------------------------------------------------------------------------
p1 <- plot_perforences_1lambda_alla(psi_10pnt, z_para=1, est_plot_only = T) +
  labs(title = "z=1") + theme(plot.title = element_text(hjust = 0.5))
p0 <- plot_perforences_1lambda_alla(psi_10pnt, z_para=0, est_plot_only = T) +
  labs(title = "z=0") + theme(plot.title = element_text(hjust = 0.5))

p <- grid.arrange(p1, p0, nrow=1,
                  top = textGrob(paste0("Average Treatment Effect,  E[Y|a,z] - E[Y|0,z]",
                                         gp=gpar(fontface = 'bold')))) # gp=gpar(fontsize=11, fontface = 'bold')


## ----simu_B_n200_rslt_u, fig.width=6, fig.height=7--------------------------------------------------------------------------
set.seed(123)
simu_results <- run_simu_rep(generate_data_1, n=nn, B=1000, return_all_rslts=T,  undersmooth=T)
cat("z=1:")
p1_list <- plot_perforences_1lambda_alla(simu_results$result_summary, z_para=1, plot_list=F)
cat("z=0:")
p0_list <- plot_perforences_1lambda_alla(simu_results$result_summary, z_para=0, plot_list=F)
save.image(file=here("data", "rdata", "02_simulation_V3_1_200_U.RData"))


## ----simu_B_n200_rslt_grid--------------------------------------------------------------------------------------------------
set.seed(123)
lambda_scalers <- c(1.2, 1.1, seq(from=1, to=0.02, length=28))
simu_results_lists <- list()
for(i in 1:length(lambda_scalers)){
  scaler = lambda_scalers[i]
  simu_results_lists[[i]] <- run_simu_rep(generate_data_1, n=nn, B=1000, lambda_scaler=scaler, return_all_rslts=F,  undersmooth=F)
}
simu_results_all <- do.call("rbind", simu_results_lists) %>% as.data.frame()

save.image(file=here("data", "rdata", "02_simulation_V3_1_200_grid.RData"))


## ----simu_B_n200_rslt_grid_05_1---------------------------------------------------------------------------------------------
plot_perforences_alllambda_1a(simu_results_all, a = 0.5, z = 1)


## ---------------------------------------------------------------------------------------------------------------------------
knitr::purl("2_simulations.Rmd")

