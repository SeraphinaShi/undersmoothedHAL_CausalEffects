

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
  x_names = names(obs)[names(obs) != 'Y']
  y_type = "binomial"

  Y <- as.numeric(as.matrix(obs %>% select(all_of(y_name))))
  X <- obs %>% 
    select(all_of(x_names)) %>% 
    mutate_if(sapply(., is.factor), as.numeric)
  
  
  fit_hal_all_criteria_rslts <- fit_hal_all_criteria(X, Y, y_type, eval_points)
  
  #================================CV-HAL================================
  lambda_CV <- fit_hal_all_criteria_rslts$lambda_list$lambda_CV
  psi_hat <- sapply(eval_points, function(a){ X_new <- X
                                              X_new$A = a
                                              mean(predict(fit_hal_all_criteria_rslts$hal_fit_list$hal_CV, new_data = X_new)) } )

  psi_hat_pnt_cv <- cbind(eval_points, matrix(psi_hat, ncol=1), lambda_CV, 1, fit_hal_all_criteria_rslts$hal_fit_time_list$hal_cv_fit_time)
  colnames(psi_hat_pnt_cv) <- c("a", "Y_hat", "lambda", "lambda_scaler", "hal_fit_time")

    # IC-based inference
  psi_hat_pnt_cv_se <- IC_based_se(X, Y, fit_hal_all_criteria_rslts$hal_fit_list$hal_CV, eval_points)
  psi_hat_pnt_cv <- as.data.frame(psi_hat_pnt_cv) %>% 
    mutate(SE = psi_hat_pnt_cv_se,
           ci_lwr = Y_hat - 1.96 * SE,
           ci_upr = Y_hat + 1.96 * SE)

  # bootstrap-based inference
  psi_hat_pnt_cv_bt_bds <- bootstrap_inference(X, Y, eval_points, fit_hal_all_criteria_rslts$hal_fit_list$hal_CV, y_type)
  psi_hat_pnt_cv$ci_lwr_bt <- psi_hat_pnt_cv_bt_bds$lower_bd
  psi_hat_pnt_cv$ci_upr_bt <- psi_hat_pnt_cv_bt_bds$upper_bd
  psi_hat_pnt_cv$SE_bt <- psi_hat_pnt_cv_bt_bds$SE
  
  
  #================================global undersmoothing================================
  if(any(fit_hal_all_criteria_rslts$hal_fit_list$hal_u_g$coefs[-1] != 0)){
    
    psi_hat <- sapply(eval_points, function(a){ X_new <- X
    X_new$A = a
    mean(predict(fit_hal_all_criteria_rslts$hal_fit_list$hal_u_g, new_data = X_new)) } )
    
    lambda_scaler = fit_hal_all_criteria_rslts$lambda_list$lambda_u_g / lambda_CV
    psi_hat_pnt_u_g <- cbind(eval_points, matrix(psi_hat, ncol=1), fit_hal_all_criteria_rslts$lambda_list$lambda_u_g, lambda_scaler, fit_hal_all_criteria_rslts$hal_fit_time_list$hal_u_g_fit_time)
    colnames(psi_hat_pnt_u_g) <- c("a", "Y_hat", "lambda", "lambda_scaler", "hal_fit_time")
    
    # IC-based inference
    psi_hat_pnt_u_g_se <- IC_based_se(X, Y, fit_hal_all_criteria_rslts$hal_fit_list$hal_CV, eval_points)
    psi_hat_pnt_u_g <- as.data.frame(psi_hat_pnt_u_g) %>% 
      mutate(SE = psi_hat_pnt_u_g_se,
             ci_lwr = Y_hat - 1.96 * SE,
             ci_upr = Y_hat + 1.96 * SE)
    
    # bootstrap-based inference
    psi_hat_pnt_u_g_bt_bds <- bootstrap_inference(X, Y, eval_points, fit_hal_all_criteria_rslts$hal_fit_list$hal_u_g, y_type)
    psi_hat_pnt_u_g$ci_lwr_bt <- psi_hat_pnt_u_g_bt_bds$lower_bd
    psi_hat_pnt_u_g$ci_upr_bt <- psi_hat_pnt_u_g_bt_bds$upper_bd
    psi_hat_pnt_u_g$SE_bt <- psi_hat_pnt_u_g_bt_bds$SE
    
  } else {
    psi_hat_pnt_u_g <- NA
  }


  #================================local undersmoothing================================

  #--------------------------estimations---------------------------------
  hal_u_l = fit_hal_all_criteria_rslts$hal_fit_list$hal_u_l
  lambda_u_l = fit_hal_all_criteria_rslts$lambda_list$lambda_u_l
  lambda_u_l_idx = fit_hal_all_criteria_rslts$lambda_u_l_idx
  
  psi_hat_pnt_u_l = matrix(ncol = 8, nrow = length(eval_points))
  
  for (i in 1:length(eval_points)) {
    a = eval_points[i]
    u_l_idx = lambda_u_l_idx[i]
    
    hal_fit = hal_u_l[[u_l_idx]]
    lambda = lambda_u_l[u_l_idx]
    lambda_scaler = lambda / lambda_CV
    
    coef <- hal_fit$coefs
    basis_mat <- cbind(1, as.matrix(hal_fit$x_basis))
    
    nonzero_idx <- which(coef != 0)
    
    if(length(nonzero_idx) > 0){
      
      X_new <- X
      X_new$A = a
      Y_hat <- mean(predict(hal_fit, new_data = X_new))

      # IC-based inference
      SE <- IC_based_se_u_l(X, Y, hal_fit, a)
      ci_lwr <- Y_hat - 1.96 * SE
      ci_upr <- Y_hat + 1.96 * SE
      
      psi_hat_pnt_u_l[i,] = c(a, Y_hat, lambda, lambda_scaler, fit_hal_all_criteria_rslts$hal_fit_time_list$hal_u_l_fit_time, SE, ci_lwr, ci_upr)
    } else {
      psi_hat_pnt_u_l[i,] <- c(a, NA, lambda, lambda_scaler, fit_hal_all_criteria_rslts$hal_fit_time_list$hal_u_l_fit_time, NA, NA, NA)
    }
  }
  colnames(psi_hat_pnt_u_l) <- c("a", "y_hat", 'lambda', 'lambda_scaler', 'hal_fit_time', "SE", "ci_lwr", "ci_upr")
  psi_hat_pnt_u_l <- as.data.frame(psi_hat_pnt_u_l)
  
  # bootstrap-based inference
  psi_hat_pnt_u_l_bt_bds <- bootstrap_inference_u_l(X, Y, eval_points, fit_hal_all_criteria_rslts$hal_fit_list$hal_u_l, y_type, lambda_u_l_idx)
  psi_hat_pnt_u_l$ci_lwr_bt <- psi_hat_pnt_u_l_bt_bds$lower_bd
  psi_hat_pnt_u_l$ci_upr_bt <- psi_hat_pnt_u_l_bt_bds$upper_bd
  psi_hat_pnt_u_l$SE_bt <- psi_hat_pnt_u_l_bt_bds$SE
  
  
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
run_simu_rep <- function(gen_data_functions, n, rounds, return_all_rslts = F){
  
  result_list <- list()
  
  for(r in 1:rounds){
    print(paste0("round ", r))
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
    
    result_list[[r]] <- result
  }
  
  results <- list()
  for (method in c("CV", "U_G", "U_L")){
    result_list_method <- lapply(result_list, function(lst) lst[[method]])
    result_all <-  do.call("rbind", result_list_method) %>% as.data.frame()
    result_all <- merge(as.data.frame(psi0_pnt), result_all, by=c("a"))
    
    result_summary <- result_all %>% 
      filter(SE != 0) %>% 
      mutate(bias = abs(y_hat - psi0),
             bias_se_ratio = bias / SE,
             bias_se_ratio_bt = bias / SE_bt,
             cover_rate = as.numeric(ci_lwr <= psi0 & psi0 <= ci_upr),
             cover_rate_bt = as.numeric(ci_lwr_bt <= psi0 & psi0 <= ci_upr_bt)) %>% 
      group_by(a) %>% 
      mutate(oracal_SE = sqrt(var(y_hat)),
             oracal_bias_se_ratio = bias / oracal_SE,
             oracal_ci_lwr = y_hat - 1.96 * oracal_SE,
             oracal_ci_upr = y_hat + 1.96 * oracal_SE,
             oracal_cover_rate = as.numeric(oracal_ci_lwr <= psi0 & psi0 <= oracal_ci_upr)) %>%
      summarise(across(where(is.numeric), mean)) %>% 
      ungroup() %>%
      mutate(hal_fit_time_unit = 'secs',
             method = method)
    
    if(return_all_rslts){
      results[[method]] <- list(result_summary = result_summary,
                                all_results = result_list_method)
    } else {
      results[[method]] <- list(result_summary = result_summary)
    }
  }
  
  results$result_summary <- rbind(results$CV$result_summary, results$U_G$result_summary, results$U_L$result_summary)
  
  return(results)
}







##############################################################


run_simu_1round_scalers <- function(gen_data_functions, n, lambda_scalers){
  
  obs <- gen_data_functions(n)
  eval_points = seq(0,5,0.5)
  
  y_name = "Y"
  x_names = names(obs)[names(obs) != 'Y']
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
  # print(sprintf('  CV lambda: %f', lambda_CV))
  
  end <- Sys.time()
  
  hal_fit_time = end - start
  if(units(hal_fit_time) == 'secs'){
    hal_fit_time = as.numeric(hal_fit_time)
  } else if (units(hal_fit_time) == 'mins'){
    hal_fit_time = as.numeric(hal_fit_time)  * 60
  } else if (units(hal_fit_time) == 'hours'){
    hal_fit_time = as.numeric(hal_fit_time)  * 60 * 24
  }
  
  hal_cv_fit_time = hal_fit_time
  
  
  #================================undersmoothing with scalers================================
  
  results = list()
  for (i in 1:length(lambda_scalers)) {
    start <- Sys.time()
    
    lambda_scaler = lambda_scalers[i]
    # print(sprintf('     %i, lambda_scaler: %f', i, lambda_scaler))
    lambda = lambda_CV * lambda_scaler
    
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
    
    end <- Sys.time()
    
    hal_fit_time = end - start
    if(units(hal_fit_time) == 'secs'){
      hal_fit_time = as.numeric(hal_fit_time)
    } else if (units(hal_fit_time) == 'mins'){
      hal_fit_time = as.numeric(hal_fit_time)  * 60
    } else if (units(hal_fit_time) == 'hours'){
      hal_fit_time = as.numeric(hal_fit_time)  * 60 * 24
    }
    
    hal_scaled_fit_time = hal_cv_fit_time + hal_fit_time
    
    #--------------------------estimations---------------------------------
    coef <- hal_fit$coefs
    basis_mat <- cbind(1, as.matrix(hal_fit$x_basis))
    
    nonzero_idx <- which(coef != 0)
    
    if(length(nonzero_idx) > 0) {
      
      psi_hat <- sapply(eval_points, function(a){ X_new <- X
                                                  X_new$A = a
                                                  mean(predict(hal_fit, new_data = X_new)) } )
      
      psi_hat_pnt_scaled <- cbind(eval_points, matrix(psi_hat, ncol=1), lambda_CV, lambda_scaler, hal_scaled_fit_time)
      colnames(psi_hat_pnt_scaled) <- c("a", "Y_hat", "lambda", "lambda_scaler", "hal_fit_time")
      
      # IC-based inference
      psi_hat_pnt_scaled_se <- IC_based_se(X, Y, hal_fit, eval_points)
      psi_hat_pnt_scaled <- as.data.frame(psi_hat_pnt_scaled) %>% 
        mutate(SE = psi_hat_pnt_scaled_se,
               ci_lwr = Y_hat - 1.96 * SE,
               ci_upr = Y_hat + 1.96 * SE)
      
      
      # bootstrap-based inference
      psi_hat_pnt_scaled_bt_bds <- bootstrap_inference(X, Y, eval_points, hal_fit, y_type)
      psi_hat_pnt_scaled$ci_lwr_bt <- psi_hat_pnt_scaled_bt_bds$lower_bd
      psi_hat_pnt_scaled$ci_upr_bt <- psi_hat_pnt_scaled_bt_bds$upper_bd
      psi_hat_pnt_scaled$SE_bt <- psi_hat_pnt_scaled_bt_bds$SE
      
    } else {
      psi_hat_pnt_scaled = as.data.frame(cbind(eval_points, NA, lambda_CV, lambda_scaler, hal_scaled_fit_time, NA, NA, NA, NA, NA, NA))
      names(psi_hat_pnt_scaled) = c("a", "Y_hat", "lambda", "lambda_scaler", "hal_fit_time", "SE", "ci_lwr", "ci_upr", "ci_lwr_bt", "ci_upr_bt", "SE_bt")
    }
    results[[i]] = psi_hat_pnt_scaled
  }
  
  names(results) = paste0("scale=", round(lambda_scalers, 4))
  
  return(results)
}


##############################################################

run_simu_scaled_rep <- function(gen_data_functions, n, rounds, return_all_rslts = F){
  lambda_scalers = c(1.2, 1.1, 10^seq(from=0, to=-3, length=20))
  result_list <- list()
  for(r in 1:rounds){
    print(paste0("round ", b))
    result <- tryCatch({
      run_simu_1round_scalers(gen_data_functions, n=n, lambda_scalers=lambda_scalers)
    }, error = function(e) {
      print(paste0("Error: ", e$message))
      NULL
    })
    
    while(is.null(result)) {
      print('retry with a new generated data')
      result <- tryCatch({
        run_simu_1round_scalers(gen_data_functions, n=n, lambda_scalers=lambda_scalers)
      }, error = function(e) {
        print(paste0("Error: ", e$message))
        NULL
      })
    }
    
    result_list[[r]] <- result
  }
  
  results <- list()
  no_empirical_CI_proportion <- c()
  
  for (i in 1:length(lambda_scalers)){
    
    lambda_scaler = lambda_scalers[i]
    
    result_list_scale <- lapply(result_list, function(lst) lst[[i]])
    no_empirical_CI_proportion[i] <- mean(sapply(result_list_scale, function(rlt) any(is.na(rlt[,colnames(rlt) == 'SE']))))
    result_all <-  do.call("rbind", result_list_scale) %>% as.data.frame()
    result_all <- merge(as.data.frame(psi0_pnt), result_all, by=c("a"))
    
    result_summary <- result_all %>% 
      filter(SE != 0) %>% 
      mutate(bias = abs(y_hat - psi0),
             bias_se_ratio = bias / SE,
             bias_se_ratio_bt = bias / SE_bt,
             cover_rate = as.numeric(ci_lwr <= psi0 & psi0 <= ci_upr),
             cover_rate_bt = as.numeric(ci_lwr_bt <= psi0 & psi0 <= ci_upr_bt)) %>% 
      group_by(a) %>% 
      mutate(oracal_SE = sqrt(var(y_hat)),
             oracal_bias_se_ratio = bias / oracal_SE,
             oracal_ci_lwr = y_hat - 1.96 * oracal_SE,
             oracal_ci_upr = y_hat + 1.96 * oracal_SE,
             oracal_cover_rate = as.numeric(oracal_ci_lwr <= psi0 & psi0 <= oracal_ci_upr)) %>%
      summarise(across(where(is.numeric), mean)) %>% 
      ungroup() %>%
      mutate(hal_fit_time_unit = 'secs',
             method = 'scale')
    
    if(return_all_rslts){
      results[[paste0("scale=", round(lambda_scaler, 4))]] <- list(result_summary = result_summary,
                                                                   all_results = result_list_scale)
    } else {
      results[[paste0("scale=", round(lambda_scaler, 4))]] <- list(result_summary = result_summary)
    }
  }
  
  result_summary <- results[[1]]$result_summary
  for (i in 1:length(lambda_scalers)) {
    result_summary <- rbind(result_summary, results[[i]]$result_summary)
  }
  results$result_summary <- result_summary
  
  results$no_empirical_CI_proportion <- no_empirical_CI_proportion
  
  return(results)
}

