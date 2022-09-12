library(here)
library(data.table)
library(dplyr)
library(foreach)

library(stringr)
library(glmnet)

library(origami)
library(hal9001)
library(tictoc)



generate_data1 <- function(n){
  # exogenous variables
  U_W <- runif(n, 0, 1)
  U_A <- rnorm(n,0,1)
  U_Z <- runif(n, 0, 1)
  U_Y <- runif(n, 0, 1)
  
  # endogenous variables
  W <- U_W
  A <- (10 - 10*W + U_A)*as.numeric(10 - 10*W + U_A >= 0 & 10 - 10*W + U_A <= 10) + 10 * as.numeric(10 - 10*W + U_A > 10)
  Z <- as.numeric(U_Z < plogis(0.2*(W+A)))
  Y <- as.numeric(U_Z < plogis(0.1*W + 0.3*A + 0.05 * W * A + Z))
  
  mean(Y) 
  # data frame
  O <- data.frame(W, A, Z, Y)
  return(O)
}

set.seed(123)

df <- generate_data1(1000)
y <- as.numeric(as.matrix(df %>% select(Y)))
x <- df %>% 
  select(A,Z,Y) %>% 
  mutate_if(sapply(., is.factor), as.numeric)
n <- nrow(df)


tic()
fit <- fit_hal(X = x,
               Y = y, 
               return_x_basis = TRUE,
               family = "binomial",
               num_knots = hal9001:::num_knots_generator(
                 max_degree = 3,
                 smoothness_orders = 1,
                 base_num_knots_0 = 500,
                 base_num_knots_1 = max(100, ceiling(sqrt(n)))
               )
)
toc()
