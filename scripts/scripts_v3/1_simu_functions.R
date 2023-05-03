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
  tmat <- t(X) %*% d_phi / n
  if(! is.matrix(try(solve(tmat), silent = TRUE))){
    return(NA)
  }
  tmat <- -solve(tmat)
  
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
run_simu_1round <- function(gen_data_functions, n, undersmooth='none', lambda_scaler = 1){
  obs <- gen_data_functions(n)
  
  start <- Sys.time()
  
  y_name = "Y"
  x_names = c("W", "A", "Z")
  y_type = "binomial"


  Y <- as.numeric(as.matrix(obs %>% select(all_of(y_name))))
  X <- obs %>% 
    select(all_of(x_names)) %>% 
    mutate_if(sapply(., is.factor), as.numeric)
  
  # fitting HAL
  CV_hal <- fit_hal(X = X, Y = Y, family = y_type,
                    return_x_basis = TRUE,
                    num_knots = hal9001:::num_knots_generator(
                      max_degree = ifelse(ncol(X) >= 20, 2, 3),
                      smoothness_orders = 1,
                      base_num_knots_0 = 20,
                      base_num_knots_1 = 20 # max(100, ceiling(sqrt(n)))
                    )
  )
  CV_lambda <- CV_hal$lambda_star
  
  
  if(undersmooth == 'none' & lambda_scaler == 1){
    hal_fit <- CV_hal
    lambda = CV_lambda
    
  } else {
    if(undersmooth == 'none'){
      lambda = lambda_scaler * CV_lambda
    } else if(undersmooth == 'global'){
      CV_nonzero_col <- which(CV_hal$coefs[-1] != 0)
      # if all coefs are zero, skip undersmooth and use the initial fit
      if (length(CV_nonzero_col) == 0){
        lambda = CV_lambda
      }else{
        CV_basis_mat <- as.matrix(CV_hal$x_basis)
        CV_basis_mat <- as.matrix(CV_basis_mat[, CV_nonzero_col])
        
        hal_undersmooth <- undersmooth_hal(X, Y,
                                           fit_init = CV_hal,
                                           basis_mat = CV_basis_mat,
                                           family = y_type)
        lambda = hal_undersmooth$lambda_under
      }
    } else if(undersmooth == 'local'){
      lambdas = rep(NA, 5)
      for(i in 1:5){
        a_seg_min = seq(0,5,1)[i]
        a_seg_max = seq(0,5,1)[i+1]
        obs_local = obs[a_seg_min <= obs$A & obs$A < a_seg_max, ]
        
        Y_seg <- as.numeric(as.matrix(obs_local %>% select(all_of(y_name))))
        X_seg <- obs_local %>% 
          select(all_of(x_names)) %>% 
          mutate_if(sapply(., is.factor), as.numeric)
        
        if((0.1 <= sum(Y_seg==1)/length(Y_seg) & sum(Y_seg==1)/length(Y_seg) <= 0.9) & 
           (sum(Y_seg==1) > 10) & (sum(Y_seg!=1) > 10) &
           (dim(X)[1] > 0)){
          
          # fitting HAL
          CV_hal <- fit_hal(X = X_seg, Y = Y_seg, family = y_type,
                            return_x_basis = TRUE,
                            num_knots = hal9001:::num_knots_generator(
                              max_degree = ifelse(ncol(X) >= 20, 2, 3),
                              smoothness_orders = 1,
                              base_num_knots_0 = 20,
                              base_num_knots_1 = 20 # max(100, ceiling(sqrt(n)))
                            )
          )
          CV_lambda <- CV_hal$lambda_star
          
          CV_nonzero_col <- which(CV_hal$coefs[-1] != 0)
          if (length(CV_nonzero_col) == 0){
            lambda = CV_lambda
          }else{
            CV_basis_mat <- as.matrix(CV_hal$x_basis)
            CV_basis_mat <- as.matrix(CV_basis_mat[, CV_nonzero_col])
            
            hal_undersmooth <- undersmooth_hal(X_seg, Y_seg,
                                               fit_init = CV_hal,
                                               basis_mat = CV_basis_mat,
                                               family = y_type)
            lambda = hal_undersmooth$lambda_under
          }
          lambdas[i] = lambda
        } else {
          lambdas[i] = CV_lambda
        }
      }
      print(lambdas)
      lambda = min(lambdas, na.rm = T)

    }
    
    
    hal_fit <- fit_hal(X = X, Y = Y, family = y_type,
                       return_x_basis = TRUE,
                       num_knots = hal9001:::num_knots_generator(
                         max_degree = ifelse(ncol(X) >= 20, 2, 3),
                         smoothness_orders = 1,
                         base_num_knots_0 = 20, #200
                         base_num_knots_1 = 20 # max(100, ceiling(sqrt(n)))
                       ),
                       fit_control = list(
                         cv_select = FALSE,
                         n_folds = 10,
                         foldid = NULL,
                         use_min = TRUE,
                         lambda.min.ratio = 1e-4,
                         prediction_bounds = "default"
                       ),
                       lambda = lambda)
  }

  
  end <- Sys.time()
  hal_fit_time = end - start
  hal_fit_time_unit = units(hal_fit_time)
  hal_fit_time = as.numeric(hal_fit_time)
  
  coef <- hal_fit$coefs
  basis_mat <- cbind(1, as.matrix(hal_fit$x_basis))
  
  nonzero_idx <- which(coef != 0)
  coef_nonzero <- coef[nonzero_idx]
  basis_mat_nonzero <- as.matrix(basis_mat[, nonzero_idx])
  
  IC_beta <- cal_IC_for_beta(X = basis_mat_nonzero, 
                             Y = Y, 
                             Y_hat = predict(hal_fit, new_data = X, type = "response"),
                             beta_n = coef_nonzero
  )
  
  if(any(is.na(IC_beta))){
    # print("failed to calculate IC_beta because of the given basis matrix is not invertible")
    return(NA)
  }
  
  # calculate ATE
  psi_hat_10pnt = matrix(ncol = 6)
  for (z in c(1,0)) {
    for (a in seq(0.5,5,0.5)) {
      X_new <- X
      X_new$A = a
      X_new$Z = z
      
      X_new_0 <- X_new
      X_new_0$A = 0
      
      Ya_hat <- predict(hal_fit, new_data = X_new)
      Y0_hat <- predict(hal_fit, new_data = X_new_0)
      psi_hat <- mean(Ya_hat - Y0_hat)
      
      # efficient influence curve
      x_basis_a <- make_design_matrix(as.matrix(X_new), hal_fit$basis_list, p_reserve = 0.75)
      x_basis_a_nonzero <- as.matrix(cbind(1, x_basis_a)[, nonzero_idx])
      
      x_basis_0 <- make_design_matrix(as.matrix(X_new_0), hal_fit$basis_list, p_reserve = 0.75)
      x_basis_0_nonzero <- as.matrix(cbind(1, x_basis_0)[, nonzero_idx])
      
      IC_ATE <- cal_IC_for_ATE(X_new_a = x_basis_a_nonzero, 
                               X_new_0 = x_basis_0_nonzero, 
                               beta_n = coef_nonzero, IC_beta = IC_beta)
      
      # empirical SE and 95% confidence interval
      
      SE <- sqrt(var(IC_ATE)/n)
      ci_lwr <- psi_hat - 1.96 * SE
      ci_upr <- psi_hat + 1.96 * SE
      
      psi_hat_10pnt <- rbind(psi_hat_10pnt, matrix(c(a, z, psi_hat, SE, ci_lwr, ci_upr), nrow = 1))
    }   
  }
  psi_hat_10pnt <- psi_hat_10pnt[-1, ]
  colnames(psi_hat_10pnt) <- c("a", "z", "psi_hat", "SE", "ci_lwr", "ci_upr")
  
  psi_hat_10pnt <- cbind(psi_hat_10pnt, lambda, CV_lambda, lambda_scaler, hal_fit_time, hal_fit_time_unit)

  return(psi_hat_10pnt)
}


##############################################################
# This function runs the simulation for B rounds with given:
#       - data generating function
#       - sample size n
#       - number of simulations: B

# returns the estimated ATE, empirical 95% CI, and coverage rates
##############################################################
run_simu_rep <- function(gen_data_functions, n, B, undersmooth='none', lambda_scaler=1, return_all_rslts = F){
  result_list <- list()
  for(b in 1:B){
    result <- run_simu_1round(gen_data_functions, n=n, undersmooth, lambda_scaler)
    while(any(is.na(result))){
      result <- run_simu_1round(gen_data_functions, n=n, undersmooth, lambda_scaler)
    }
    result_list[[b]] <- result
  }
  result_all <-  do.call("rbind", result_list) %>% as.data.frame()
  result_all <- merge(as.data.frame(psi0_10pnt), result_all, by=c("a", "z"))
  
  cols.num <- names(result_all)[names(result_all) != 'hal_fit_time_unit']
  result_all[cols.num] <- sapply(result_all[cols.num],as.numeric)
  
  hal_fit_time_unit = paste(unique(result_all$hal_fit_time_unit), collapse = ' & ')
  
  result_summary <- result_all %>% 
    filter(SE != 0) %>% 
    mutate(bias = abs(psi_hat - psi0),
           bias_se_ratio = bias / SE,
           cover_rate = as.numeric(ci_lwr <= psi0 & psi0 <= ci_upr)) %>% 
    group_by(a, z) %>% 
    mutate(oracal_SE = sqrt(var(psi_hat)),
           oracal_bias_se_ratio = bias / oracal_SE,
           oracal_ci_lwr = psi_hat - 1.96 * oracal_SE,
           oracal_ci_upr = psi_hat + 1.96 * oracal_SE,
           oracal_cover_rate = as.numeric(oracal_ci_lwr <= psi0 & psi0 <= oracal_ci_upr)) %>%
    summarise(across(where(is.numeric), mean)) %>% 
    ungroup() %>%
    mutate(hal_fit_time_unit = hal_fit_time_unit)
  
  if(return_all_rslts){
    return(list(result_summary = result_summary,
                all_results = result_list))
  } else {
    return(result_summary)
  }
}


