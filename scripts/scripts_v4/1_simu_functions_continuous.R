
########################
# calculating efficient influence curves for continuous outcome
########################
cal_IC_for_beta_cont <- function(X, Y, Y_hat, beta_n){

  
  if (!is.matrix(X)) X <- as.matrix(X)
  
  # 1. calculate score: X'(Y - Y_hat)
  score <- sweep(X, 1, (Y - Y_hat), `*`)

  # 2. calculate E_{P_n}(X'X)^(-1)
  d_scaler = solve(t(X) %*% X)
  
  # 3. calculate influence curves
  IC <- d_scaler %*% t(score)
  
  return(IC)
}


cal_IC_for_ATE <- function(Xa_X0_diff, IC_beta){
  
  if (!is.matrix(Xa_X0_diff)) Xa_X0_diff <- as.matrix(Xa_X0_diff)
  
  IC <- Xa_X0_diff %*% IC_beta
  
  return(IC)
}


n=250
X1 <- rnorm(n=n)
X2 <- rbinom(n=n, size=1, prob=0.2)
Y <- X1 + 5*X2 + X1*X2
df = data.frame(Y=Y, X1=X1, X2=X2)

glm_fit <- glm(data = df, formula = Y ~ X1 + X2)




beta_n = glm_fit$coefficients
X <- cbind(1, df %>% select(-Y))
Y_hat <- glm_fit$linear.predictors

IC_betas = cal_IC_for_beta_cont(X = X,
                                Y = Y,
                                Y_hat = Y_hat,
                                beta_n = beta_n)


se_IC <- sqrt(apply(IC_betas, 1, var))


library(sandwich)
# compare estimated se from calculated ICs and from glm_fit. 
cbind(se_IC, sqrt(diag(sandwich(glm_fit))))
