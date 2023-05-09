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
run_simu_1round <- function(gen_data_functions, n){
  
  obs <- gen_data_functions(n)
  eval_points = seq(0,5,0.5)
  
  y_name = "Y"
  x_names = c("W", "A")
  y_type = "binomial"


  Y <- as.numeric(as.matrix(obs %>% select(all_of(y_name))))
  X <- obs %>% 
    select(all_of(x_names)) %>% 
    mutate_if(sapply(., is.factor), as.numeric)
  
  
  #================================CV-HAL================================
  start <- Sys.time()
  
  hal_CV <- fit_hal(X = X, Y = Y, family = y_type,
                    return_x_basis = TRUE,
                    num_knots = hal9001:::num_knots_generator(
                      max_degree = ifelse(ncol(X) >= 20, 2, 3),
                      smoothness_orders = 1,
                      base_num_knots_0 = 20,
                      base_num_knots_1 = 20 # max(100, ceiling(sqrt(n)))
                    )
  )
  lambda_CV <- hal_CV$lambda_star
  print(sprintf('  CV lambda: %f', lambda_CV))
  
  end <- Sys.time()
  hal_cv_fit_time = as.numeric(end - start)
  
  #--------------------------estimations---------------------------------
  hal_fit = hal_CV
  
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
  
  # calculate Y
  psi_hat_pnt = matrix(ncol = 5)
  
  for (a in eval_points) {
    X_new <- X
    X_new$A = a
    
    Y_hat <- mean(predict(hal_fit, new_data = X_new))
    
    # efficient influence curve
    x_basis_a <- make_design_matrix(as.matrix(X_new), hal_fit$basis_list, p_reserve = 0.75)
    x_basis_a_nonzero <- as.matrix(cbind(1, x_basis_a)[, nonzero_idx])
    
    IC_EY <- cal_IC_for_EY(X_new = x_basis_a_nonzero, 
                             beta_n = coef_nonzero, IC_beta = IC_beta)
    
    # empirical SE and 95% confidence interval
    
    SE <- sqrt(var(IC_EY)/n)
    ci_lwr <- Y_hat - 1.96 * SE
    ci_upr <- Y_hat + 1.96 * SE
    
    psi_hat_pnt <- rbind(psi_hat_pnt, matrix(c(a, Y_hat, SE, ci_lwr, ci_upr), nrow = 1))
  }
  psi_hat_pnt <- psi_hat_pnt[-1, ]
  colnames(psi_hat_pnt) <- c("a", "y_hat", "SE", "ci_lwr", "ci_upr")
  
  lambda = lambda_CV
  hal_fit_time = hal_cv_fit_time
  lambda_scaler = 1
  psi_hat_pnt <- cbind(psi_hat_pnt, lambda, lambda_scaler, hal_fit_time)
  
  psi_hat_pnt_cv = psi_hat_pnt
  
  
  
  #================================global undersmoothing================================
  start <- Sys.time()
  
  CV_nonzero_col <- which(hal_CV$coefs[-1] != 0)
      # if all coefs are zero, skip undersmooth and use the initial fit
  if (length(CV_nonzero_col) == 0){
    lambda_u_g = lambda_CV
  }else{
    CV_basis_mat <- as.matrix(hal_CV$x_basis)
    CV_basis_mat <- as.matrix(CV_basis_mat[, CV_nonzero_col])
    
    hal_undersmooth <- undersmooth_hal(X, Y,
                                       fit_init = hal_CV,
                                       basis_mat = CV_basis_mat,
                                       family = y_type)
    lambda_u_g = hal_undersmooth$lambda_under
  }
  print(sprintf('  globally u lambdas: %f', lambda_u_g))
  
  hal_u_g <- fit_hal(X = X, Y = Y, family = y_type,
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
                     lambda = lambda_u_g)
  
  end <- Sys.time()
  hal_u_g_fit_time = as.numeric(end - start) + hal_cv_fit_time
  
  
  #--------------------------estimations---------------------------------
  hal_fit = hal_u_g
  
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
  
  # calculate Y
  psi_hat_pnt = matrix(ncol = 5)
  
  for (a in eval_points) {
    X_new <- X
    X_new$A = a
    
    Y_hat <- mean(predict(hal_fit, new_data = X_new))
    
    # efficient influence curve
    x_basis_a <- make_design_matrix(as.matrix(X_new), hal_fit$basis_list, p_reserve = 0.75)
    x_basis_a_nonzero <- as.matrix(cbind(1, x_basis_a)[, nonzero_idx])
    
    IC_EY <- cal_IC_for_EY(X_new = x_basis_a_nonzero, 
                           beta_n = coef_nonzero, IC_beta = IC_beta)
    
    # empirical SE and 95% confidence interval
    
    SE <- sqrt(var(IC_EY)/n)
    ci_lwr <- Y_hat - 1.96 * SE
    ci_upr <- Y_hat + 1.96 * SE
    
    psi_hat_pnt <- rbind(psi_hat_pnt, matrix(c(a, Y_hat, SE, ci_lwr, ci_upr), nrow = 1))
  }
  psi_hat_pnt <- psi_hat_pnt[-1, ]
  colnames(psi_hat_pnt) <- c("a", "y_hat", "SE", "ci_lwr", "ci_upr")
  
  lambda = lambda_u_g
  hal_fit_time = hal_u_g_fit_time
  lambda_scaler = lambda_u_g / lambda_CV
  psi_hat_pnt <- cbind(psi_hat_pnt, lambda, lambda_scaler, hal_fit_time)
  
  psi_hat_pnt_u_g = psi_hat_pnt
  
  
  
  #================================local undersmoothing================================
  start <- Sys.time()
  
  CV_nonzero_col <- which(hal_CV$coefs[-1] != 0)
  # if all coefs are zero, skip undersmooth and use the initial fit
  if (length(CV_nonzero_col) == 0){
    lambda = lambda_CV
  }else{
    
    lambdas_u_l = rep(NA, length(eval_points))
    
    for(i in 1:length(eval_points)){
      X_local <- X
      X_local$A = eval_points[i]
      
      CV_basis_mat_local <- make_design_matrix(as.matrix(X_local), hal_CV$basis_list, p_reserve = 0.75)
      CV_basis_mat_local <- as.matrix(cbind(1, CV_basis_mat_local)[, CV_nonzero_col])
      
      hal_undersmooth_local <- undersmooth_hal(X_local, Y,
                                               fit_init = hal_CV,
                                               basis_mat = CV_basis_mat_local,
                                               family = y_type)
      lambda_local = hal_undersmooth_local$lambda_under
      lambdas_u_l[i] = lambda_local
    }
    print(sprintf('  locally u lambdas: %s', paste(round(lambdas_u_l, 6), collapse = ", ")))
    
    lambda_u_l = unique(lambdas_u_l)
    lambda_u_l_idx = match(lambdas_u_l, lambda_u_l)
    
    hal_u_l = list()
    for(i in 1:length(lambda_u_l)){
      hal_u_l[[i]] <- fit_hal(X = X, Y = Y, family = y_type,
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
                              lambda = lambda_u_l[i])
    }
  }
  
  end <- Sys.time()
  hal_u_l_fit_time = as.numeric(end - start) + hal_cv_fit_time
    
  #--------------------------estimations---------------------------------
  
  psi_hat_pnt = matrix(ncol = 8)
  
  for (i in 1:length(eval_points)) {
    a = eval_points[i]
    u_l_idx = lambda_u_l_idx[i]
    
    hal_fit = hal_u_l[[u_l_idx]]
    lambda = lambda_u_l[u_l_idx]
    lambda_scaler = lambda / lambda_CV
    
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
    
    X_new <- X
    X_new$A = a
    
    Y_hat <- mean(predict(hal_fit, new_data = X_new))
    
    # efficient influence curve
    x_basis_a <- make_design_matrix(as.matrix(X_new), hal_fit$basis_list, p_reserve = 0.75)
    x_basis_a_nonzero <- as.matrix(cbind(1, x_basis_a)[, nonzero_idx])
    
    IC_EY <- cal_IC_for_EY(X_new = x_basis_a_nonzero, 
                           beta_n = coef_nonzero, IC_beta = IC_beta)
    
    # empirical SE and 95% confidence interval
    
    SE <- sqrt(var(IC_EY)/n)
    ci_lwr <- Y_hat - 1.96 * SE
    ci_upr <- Y_hat + 1.96 * SE
    
    hal_fit_time = hal_u_l_fit_time
    psi_hat_pnt <- rbind(psi_hat_pnt, matrix(c(a, Y_hat, SE, ci_lwr, ci_upr, lambda, lambda_scaler, hal_fit_time), nrow = 1))
  }
  
  psi_hat_pnt <- psi_hat_pnt[-1, ]
  colnames(psi_hat_pnt) <- c("a", "y_hat", "SE", "ci_lwr", "ci_upr", 'lambda', 'lambda_scaler', 'hal_fit_time')
  
  psi_hat_pnt_u_l = psi_hat_pnt
  
  
  #====================================================================================
  results <- list(psi_hat_pnt_cv, psi_hat_pnt_u_g, psi_hat_pnt_u_l)
  names(results) = c("CV", "U_G", "U_L")
  return(results)
}


##############################################################
# This function runs the simulation for B rounds with given:
#       - data generating function
#       - sample size n
#       - number of simulations: B

# returns the estimated ATE, empirical 95% CI, and coverage rates
##############################################################
run_simu_rep <- function(gen_data_functions, n, B, return_all_rslts = F){
  result_list <- list()
  for(b in 1:B){
    print(paste0("round ", b))
    result <- tryCatch({
      run_simu_1round(gen_data_functions, n=n)
    }, error = function(e) {
      print(paste0("Error: ", e$message))
      NULL
    })
    
    while(is.null(result)) {
      print('retry with a new generated data')
      result <- tryCatch({
        run_simu_1round(gen_data_functions, n=n)
      }, error = function(e) {
        print(paste0("Error: ", e$message))
        NULL
      })
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


