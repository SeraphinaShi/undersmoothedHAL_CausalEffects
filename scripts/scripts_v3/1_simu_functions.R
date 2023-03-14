########################
# calculating efficient influence curves
########################
cal_IC_for_beta <- function(X, Y, Y_hat, beta_n){
  n <- dim(X)[1] 
  p <- length(beta_n)
  if (!is.matrix(X)) X <- as.matrix(X)
  
  # 1. calculate score: X(Y - phi(X))
  res <- Y-Y_hat
  score <- sweep(X, 1, res, `*`)
  
  # 2. calculate the derivative of phi:
  d_phi_scaler <- as.vector(exp(- beta_n %*% t(X)) / ((1 + exp(- beta_n %*% t(X)))^2))
  d_phi <- sweep(X, 1, d_phi_scaler, `*`)
  
  # 3. -E_{P_n}(X d_phi)^(-1)
  tmat <- -solve(t(X) %*% d_phi / n)
  
  # 4. calculate influence curves
  IC <- t(tmat %*% t(score))
  
  return(IC)
}


cal_IC_for_EY <- function(X_new, beta_n, IC_beta){
  
  if (!is.matrix(X_new)) X_new <- as.matrix(X_new)
  
  d_phi_scaler_new <- as.vector(exp(- beta_n %*% t(X_new)) / ((1 + exp(- beta_n %*% t(X_new)))^2))
  d_phi_new <- sweep(X_new, 1, d_phi_scaler_new, `*`)
  
  IC = diag(d_phi_new %*% t(IC_beta))
  
  return(IC)
}


cal_IC_for_ATE <- function(X_new_a, X_new_0, beta_n, IC_beta){
  
  if (!is.matrix(X_new_a)) X_new_a <- as.matrix(X_new_a)
  if (!is.matrix(X_new_0)) X_new_0 <- as.matrix(X_new_0)
  
  d_phi_scaler_new_a <- as.vector(exp(- beta_n %*% t(X_new_a)) / ((1 + exp(- beta_n %*% t(X_new_a)))^2))
  d_phi_new_a <- sweep(X_new_a, 1, d_phi_scaler_new_a, `*`)
  
  d_phi_scaler_new_0 <- as.vector(exp(- beta_n %*% t(X_new_0)) / ((1 + exp(- beta_n %*% t(X_new_0)))^2))
  d_phi_new_0 <- sweep(X_new_0, 1, d_phi_scaler_new_0, `*`)
  
  d_phi_new <- d_phi_new_a - d_phi_new_0
  
  IC = diag(d_phi_new %*% t(IC_beta))
  
  return(IC)
}



##############################################################
# This function runs the simulation with given:
#       - data generating function
#       - sample size n

# returns the estimated ATE and empirical 95% CI
##############################################################
run_simu_1round <- function(gen_data_functions, n){
  obs <- gen_data_functions(n)
  
  y_name = "Y"
  x_names = c("W", "A", "Z")
  y_type = "binomial"
  
  Y <- as.numeric(as.matrix(obs %>% select(all_of(y_name))))
  X <- obs %>% 
    select(all_of(x_names)) %>% 
    mutate_if(sapply(., is.factor), as.numeric)
  
  # fitting HAL
  CV_hal <- undersmooth_init(X, Y, family = y_type) # get CV_lambda from HAL_fit with 0 smoothness order
  CV_lambda <- CV_hal$fit_init$lambda_star
  
  CV_hal_fit <- fit_hal(X = X,
                        Y = Y,
                        family = y_type,
                        return_x_basis = TRUE,
                        num_knots = hal9001:::num_knots_generator(
                          max_degree = ifelse(ncol(X) >= 20, 2, 3),
                          smoothness_orders = 1,
                          base_num_knots_0 = 20,
                          base_num_knots_1 = 20
                          #ceiling(sqrt(n))
                        ),
                        fit_control = list(
                          cv_select = FALSE,
                          n_folds = 10,
                          foldid = NULL,
                          use_min = TRUE,
                          lambda.min.ratio = 1e-4,
                          prediction_bounds = "default"
                        ),
                        lambda = CV_lambda
  )
  
  CV_coef <- CV_hal_fit$coefs
  basis_mat <- cbind(1, as.matrix(CV_hal_fit$x_basis))
  
  nonzero_idx <- which(CV_coef != 0)
  CV_coef_nonzero <- CV_coef[nonzero_idx]
  basis_mat_nonzero <- as.matrix(basis_mat[, nonzero_idx])
  
  IC_beta <- cal_IC_for_beta(X = basis_mat_nonzero, 
                             Y = Y, 
                             Y_hat = predict(CV_hal_fit, new_data = X, type = "response"),
                             beta_n = CV_coef_nonzero
  )
  
  # calculate ATE
  psi_hat_10pnt = matrix(ncol = 6)
  for (z in c(1,0)) {
    for (a in seq(0.5,5,0.5)) {
      X_new <- X
      X_new$A = a
      X_new$Z = z
      
      X_new_0 <- X_new
      X_new_0$A = 0
      
      Ya_hat <- predict(CV_hal_fit, new_data = X_new)
      Y0_hat <- predict(CV_hal_fit, new_data = X_new_0)
      psi_hat <- mean(Ya_hat - Y0_hat)
      
      # efficient influence curve
      x_basis_a <- make_design_matrix(as.matrix(X_new), CV_hal_fit$basis_list, p_reserve = 0.75)
      x_basis_a_nonzero <- as.matrix(cbind(1, x_basis_a)[, nonzero_idx])
      
      x_basis_0 <- make_design_matrix(as.matrix(X_new_0), CV_hal_fit$basis_list, p_reserve = 0.75)
      x_basis_0_nonzero <- as.matrix(cbind(1, x_basis_0)[, nonzero_idx])
      
      IC_ATE <- cal_IC_for_ATE(X_new_a = x_basis_a_nonzero, 
                               X_new_0 = x_basis_0_nonzero, 
                               beta_n = CV_coef_nonzero, IC_beta = IC_beta)
      
      # empirical SE and 95% confidence interval
      
      SE <- sqrt(var(IC_ATE)/n)
      CI_lwr <- psi_hat - 1.96 * SE
      CI_upr <- psi_hat + 1.96 * SE
      
      psi_hat_10pnt <- rbind(psi_hat_10pnt, matrix(c(a, z, psi_hat, SE, CI_lwr, CI_upr), nrow = 1))
    }   
  }
  psi_hat_10pnt <- psi_hat_10pnt[-1, ]
  colnames(psi_hat_10pnt) <- c("a", "z", "psi_hat", "SE", "CI_lwr", "CI_upr")
  
  return(psi_hat_10pnt)
}


##############################################################
# This function runs the simulation for B rounds with given:
#       - data generating function
#       - sample size n

# returns the estimated ATE and empirical 95% CI
##############################################################
run_simu_Bround
