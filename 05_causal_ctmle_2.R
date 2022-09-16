## ----load_lib, include = TRUE, warning=FALSE, message=FALSE, echo=FALSE------------------------------------------
# setwd("/home/seraphinashi/MissingDataInTheICU")
# install.packages("polspline")
library(dplyr)
library(data.table)
library(origami)
library(sl3)
library(ggplot2)
library(forcats)
library(here)
library(future)
library(ROCR) 
library(gridExtra)
library(tmle3)
library("ctmle")


## ----setup, include = FALSE--------------------------------------------------------------------------------------
plotFolder <- here("results","images","VIA")
if(!file.exists(plotFolder)) dir.create(plotFolder,recursive=TRUE)

set.seed(123)



## ----------------------------------------------------------------------------------------------------------------
load(here("data", "rdata", "02_2_varFusion.RData"))
df_imp <- read.csv(here("data", "varFusion", "data_imputed.csv"))  %>% dplyr::select(-X)


## ----------------------------------------------------------------------------------------------------------------
subsample <- F

set.seed(123)

if(subsample){
  icus_1 <- unique(df_imp$ICUSTAY_ID[df_imp$Y_t_plus_1 == 1])
  icus_0 <- unique(df_imp$ICUSTAY_ID[df_imp$Y_t_plus_1 == 0])
  
  icus_0 <- sample(icus_0, size = 0.5*length(icus_1))
  
  rows <- df_imp$ICUSTAY_ID %in% c(icus_1, icus_0)
  df_imp <- df_imp[rows,]
}


## ----------------------------------------------------------------------------------------------------------------
dat <- df_imp %>% 
  group_by(ICUSTAY_ID) %>% 
  arrange(ICUSTAY_ID, half_day) %>%
  mutate(A_t_plus_1 = n_BP[row_number()+1],
         Y_t_plus_2 = Y_t_plus_1[row_number()+1]) %>%
  ungroup()
# View(dat %>% select(ICUSTAY_ID, half_day, n_BP, A_t_plus_1, Y_t_plus_1, Y_t_plus_2))

dat <- dat %>% 
  filter(!is.na(Y_t_plus_2))


## ----------------------------------------------------------------------------------------------------------------
# BP measurement rate in each time block
# A = cut(dat$A_t_plus_1,
#         breaks = c(-1,0,12,17,101),
#         labels = c("0","1-12","13-17",">17"))
A = ifelse(dat$A_t_plus_1==0, 0, 1)

# mortality in the next time block
Y = dat$Y_t_plus_2

confounders_names <- list(base = var_name_list$baselines,
                          all_base = c(var_name_list$vitals,
                                       var_name_list$labs,
                                       var_name_list$num_baselines,
                                       var_name_list$num_vitals,
                                       var_name_list$num_labs),
                          all_vital = c(var_name_list$baselines,
                                        var_name_list$labs,
                                        var_name_list$num_baselines,
                                        var_name_list$num_vitals,
                                        var_name_list$num_labs),
                          all_labs = c(var_name_list$baselines,
                                       var_name_list$vitals,
                                       var_name_list$num_baselines,
                                       var_name_list$num_vitals,
                                       var_name_list$num_labs),
                          all_n_base = c(var_name_list$baselines,
                                         var_name_list$vitals,
                                         var_name_list$labs,
                                         var_name_list$num_vitals,
                                         var_name_list$num_labs),
                          all_n_vital = c(var_name_list$baselines,
                                          var_name_list$vitals,
                                          var_name_list$labs,
                                          var_name_list$num_baselines,
                                          var_name_list$num_labs),
                          all_n_lab = c(var_name_list$baselines,
                                        var_name_list$vitals,
                                        var_name_list$labs,
                                        var_name_list$num_baselines,
                                        var_name_list$num_vitals),
                          all = c(var_name_list$baselines,
                                  var_name_list$vitals,
                                  var_name_list$labs,
                                  var_name_list$num_baselines,
                                  var_name_list$num_vitals,
                                  var_name_list$num_labs))

ctmle_discrete_fits <- list()
for (i in 1:length(confounders_names)) {
  print(i)
  confounders_names_i <- confounders_names[[i]]
  confounders <- dat %>% select(all_of(confounders_names_i))
  
  ICUSTAY_ID <- dat$ICUSTAY_ID
  half_day <- as.factor(dat$half_day+1)
  
  df <- data.frame(ICUSTAY_ID, half_day, Y, A, confounders)
  n <- nrow(df)
  print(n)
  print(names(df))
  
  # With initial estimate of Q
  N <- nrow(df)
  Q <- cbind(rep(mean(Y[A == 0]), N), rep(mean(Y[A == 1]), N))
  print(head(Q))
  
  folds <- make_folds(nrow(df), cluster_ids = df$ICUSTAY_ID)
  
  ctmle_discrete_i <- ctmleDiscrete(Y = Y, A = A, 
                                    W = data.frame(confounders), Q = Q,
                                    family = "binomial", 
                                    folds = folds,
                                    preOrder = FALSE, detailed = TRUE)
  
  print(summary(ctmle_discrete_i))
  
  ctmle_discrete_fits[[i]] = ctmle_discrete_i
}


# knitr::purl("scripts/05_causal_ctmle.Rmd")
save.image(file=here("data", "rdata", "05_causal_ctmle.RData"))

