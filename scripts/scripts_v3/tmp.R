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



load(here("data", "rdata", "02_simulation_1_3_200.RData"))

results_200[[1]] <- results_200[[1]]   %>% 
  select(targ_par, psi0, mean_under, mean_init, bias_u, bias_i, sd_u, sd_i, cov_rate_u, cov_rate_i, mse_u, mse_i) 

results_200[[2]] <- results_200[[2]]   %>% 
  select(targ_par, psi0, mean_under, mean_init, bias_u, bias_i, sd_u, sd_i, cov_rate_u, cov_rate_i, mse_u, mse_i) 

save.image(file=here("data", "rdata", "02_simulation_1_3_200.RData"))
rm(results_200)



load(here("data", "rdata", "02_simulation_1_3_500.RData"))

results_500[[1]] <- results_500[[1]]   %>% 
  select(targ_par, psi0, mean_under, mean_init, bias_u, bias_i, sd_u, sd_i, cov_rate_u, cov_rate_i, mse_u, mse_i) 

results_500[[2]] <- results_500[[2]]   %>%
  select(targ_par, psi0, mean_under, mean_init, bias_u, bias_i, sd_u, sd_i, cov_rate_u, cov_rate_i, mse_u, mse_i) 

save.image(file=here("data", "rdata", "02_simulation_1_3_500.RData"))
rm(results_500)



load(here("data", "rdata", "02_simulation_1_3_1000.RData"))

results_1000[[1]] <- results_1000[[1]]   %>%
  select(targ_par, psi0, mean_under, mean_init, bias_u, bias_i, sd_u, sd_i, cov_rate_u, cov_rate_i, mse_u, mse_i) 

results_1000[[2]] <- results_1000[[2]]   %>% 
  select(targ_par, psi0, mean_under, mean_init, bias_u, bias_i, sd_u, sd_i, cov_rate_u, cov_rate_i, mse_u, mse_i) 

save.image(file=here("data", "rdata", "02_simulation_1_3_1000.RData"))
rm(results_1000)