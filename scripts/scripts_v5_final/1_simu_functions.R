library(sl3)

###############################################################################
#'  fit and return CV HAL, globally undersmoothed HAL, and locally undersmoothed HAL


fit_hal_all_criteria <- function(X, Y, y_type, eval_points){
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
  
  hal_fit_time = end - start
  if(units(hal_fit_time) == 'secs'){
    hal_fit_time = as.numeric(hal_fit_time)
  } else if (units(hal_fit_time) == 'mins'){
    hal_fit_time = as.numeric(hal_fit_time)  * 60
  } else if (units(hal_fit_time) == 'hours'){
    hal_fit_time = as.numeric(hal_fit_time)  * 60 * 24
  }
  
  hal_cv_fit_time = hal_fit_time
  num_basis_CV <- sum(hal_CV$coefs[-1] != 0)
  
  #================================global undersmoothing================================
  start <- Sys.time()
  
  CV_nonzero_col <- which(hal_CV$coefs[-1] != 0)
  # if all coefs are zero, skip undersmooth and use the initial fit
  if (length(CV_nonzero_col) == 0){
    lambda_u_g = lambda_CV
    
    print(sprintf('  globally u lambdas: %f (same as CV-HAL)', lambda_u_g))
    hal_u_g <- hal_CV
  }else{
    CV_basis_mat <- as.matrix(hal_CV$x_basis)
    CV_basis_mat <- as.matrix(CV_basis_mat[, CV_nonzero_col])
    
    hal_undersmooth <- undersmooth_hal(X, Y,
                                       fit_init = hal_CV,
                                       family = y_type)
    lambda_u_g = hal_undersmooth$lambda_under
    
    print(sprintf('  globally u lambdas: %f', lambda_u_g))
    if(is.na(lambda_u_g)){
      lambda_u_g <- lambda_CV
      print(sprintf('  globally u lambdas: %f', lambda_u_g))
      
      hal_u_g <- hal_CV
    } else {
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
    }
  }
  
  end <- Sys.time()
  
  hal_fit_time = end - start
  if(units(hal_fit_time) == 'secs'){
    hal_fit_time = as.numeric(hal_fit_time)
  } else if (units(hal_fit_time) == 'mins'){
    hal_fit_time = as.numeric(hal_fit_time)  * 60
  } else if (units(hal_fit_time) == 'hours'){
    hal_fit_time = as.numeric(hal_fit_time)  * 60 * 24
  }
  
  hal_u_g_fit_time = hal_fit_time + hal_cv_fit_time
  
  num_basis_u_g <- sum(hal_u_g$coefs[-1] != 0)
  
  #================================local undersmoothing================================
  start <- Sys.time()
  
  CV_nonzero_col <- which(hal_CV$coefs[-1] != 0)
  # if all coefs are zero, skip undersmooth and use the initial fit
  if (length(CV_nonzero_col) == 0){
    lambdas_u_l = list()
    hal_u_l = list()
    
    lambdas_u_l[[1]] = lambda_CV
    
    print(sprintf('  locally u lambdas: %s (same as CV-HAL)', lambdas_u_l))
    hal_u_l[[1]] <- hal_CV
    lambda_u_l_idx = rep(1, length(eval_points))
    
  }else{
    
    lambdas_u_l = rep(NA, length(eval_points))
    
    for(i in 1:length(eval_points)){
      X_local <- X
      X_local$A = eval_points[i]
      
      CV_basis_mat_local <- make_design_matrix(as.matrix(X_local), hal_CV$basis_list, p_reserve = 0.75)
      CV_basis_mat_local <- as.matrix(cbind(1, CV_basis_mat_local)[, CV_nonzero_col])
      
      hal_undersmooth_local <- undersmooth_hal(X_local, Y,
                                               fit_init = hal_CV,
                                               family = y_type)
      lambda_local = hal_undersmooth_local$lambda_under
      lambdas_u_l[i] = lambda_local
    }
    print(sprintf('  locally u lambdas: %s', paste(round(lambdas_u_l, 6), collapse = ", ")))
    if(any(is.na(lambdas_u_l))){
      lambdas_u_l[is.na(lambdas_u_l)] = lambda_CV
      print(sprintf('  locally u lambdas: %s', paste(round(lambdas_u_l, 6), collapse = ", ")))
    }
    
    lambda_u_l = unique(lambdas_u_l)
    lambda_u_l_idx = match(lambdas_u_l, lambda_u_l)
    
    hal_u_l = list()
    num_basis_u_l <- c()
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
      num_basis_u_l[i] <- sum(hal_u_l[[i]]$coefs[-1] != 0)
    }
  }
  
  end <- Sys.time()
  
  hal_fit_time = end - start
  if(units(hal_fit_time) == 'secs'){
    hal_fit_time = as.numeric(hal_fit_time)
  } else if (units(hal_fit_time) == 'mins'){
    hal_fit_time = as.numeric(hal_fit_time)  * 60
  } else if (units(hal_fit_time) == 'hours'){
    hal_fit_time = as.numeric(hal_fit_time)  * 60 * 24
  }
  
  hal_u_l_fit_time = hal_fit_time + hal_cv_fit_time
  num_basis_u_l = num_basis_u_l[lambda_u_l_idx]
  
  #================================return================================
  hal_fit_list <- list(hal_CV = hal_CV, hal_u_g = hal_u_g, hal_u_l = hal_u_l)
  lambda_list <- list(lambda_CV = lambda_CV, lambda_u_g = lambda_u_g, lambda_u_l = lambda_u_l)
  hal_fit_time_list <- list(hal_cv_fit_time = hal_cv_fit_time, hal_u_g_fit_time = hal_u_g_fit_time, hal_u_l_fit_time = hal_u_l_fit_time)
  num_basis_list <- list(num_basis_CV = num_basis_CV, num_basis_u_g = num_basis_u_g, num_basis_u_l = num_basis_u_l)
  
  return(list(hal_fit_list = hal_fit_list, lambda_list = lambda_list, hal_fit_time_list = hal_fit_time_list, lambda_u_l_idx = lambda_u_l_idx, num_basis_list = num_basis_list))
}


fit_hal_CV_U <- function(X, Y, y_type, eval_points){
  #================================CV-HAL================================
  print("1")
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
  
  hal_fit_time = end - start
  if(units(hal_fit_time) == 'secs'){
    hal_fit_time = as.numeric(hal_fit_time)
  } else if (units(hal_fit_time) == 'mins'){
    hal_fit_time = as.numeric(hal_fit_time)  * 60
  } else if (units(hal_fit_time) == 'hours'){
    hal_fit_time = as.numeric(hal_fit_time)  * 60 * 24
  }
  
  hal_cv_fit_time = hal_fit_time
  num_basis_CV <- sum(hal_CV$coefs[-1] != 0)
  
  #================================global undersmoothing================================
  start <- Sys.time()
  
  CV_nonzero_col <- which(hal_CV$coefs[-1] != 0)
  # if all coefs are zero, skip undersmooth and use the initial fit
  if (length(CV_nonzero_col) == 0){
    lambda_u_g = lambda_CV
    
    print(sprintf('  globally u lambdas: %f (same as CV-HAL)', lambda_u_g))
    hal_u_g <- hal_CV
  }else{
    CV_basis_mat <- as.matrix(hal_CV$x_basis)
    CV_basis_mat <- as.matrix(CV_basis_mat[, CV_nonzero_col])
    
    hal_undersmooth <- undersmooth_hal(X, Y,
                                       fit_init = hal_CV,
                                       family = y_type)
    lambda_u_g = hal_undersmooth$lambda_under
    
    print(sprintf('  globally u lambdas: %f', lambda_u_g))
    if(is.na(lambda_u_g)){
      lambda_u_g <- lambda_CV
      print(sprintf('  globally u lambdas: %f', lambda_u_g))
      
      hal_u_g <- hal_CV
    } else {
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
    }
  }
  
  end <- Sys.time()
  
  hal_fit_time = end - start
  if(units(hal_fit_time) == 'secs'){
    hal_fit_time = as.numeric(hal_fit_time)
  } else if (units(hal_fit_time) == 'mins'){
    hal_fit_time = as.numeric(hal_fit_time)  * 60
  } else if (units(hal_fit_time) == 'hours'){
    hal_fit_time = as.numeric(hal_fit_time)  * 60 * 24
  }
  
  hal_u_g_fit_time = hal_fit_time + hal_cv_fit_time
  
  num_basis_u_g <- sum(hal_u_g$coefs[-1] != 0)
  
  #================================return================================
  hal_fit_list <- list(hal_CV = hal_CV, hal_u_g = hal_u_g)
  lambda_list <- list(lambda_CV = lambda_CV, lambda_u_g = lambda_u_g)
  hal_fit_time_list <- list(hal_cv_fit_time = hal_cv_fit_time, hal_u_g_fit_time = hal_u_g_fit_time)
  num_basis_list <- list(num_basis_CV = num_basis_CV, num_basis_u_g = num_basis_u_g)
  
  return(list(hal_fit_list = hal_fit_list, lambda_list = lambda_list, hal_fit_time_list = hal_fit_time_list, num_basis_list = num_basis_list))
}

fit_hal_CV_U_0 <- function(X, Y, y_type, eval_points){
  print("0")
  
  #================================CV-HAL================================
  
  start <- Sys.time()
  
  hal_CV <- fit_hal(X = X, Y = Y, family = y_type,
                    smoothness_orders = 0,
                    return_x_basis = TRUE
  )
  lambda_CV <- hal_CV$lambda_star
  print(sprintf('  CV lambda: %f', lambda_CV))
  
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
  num_basis_CV <- sum(hal_CV$coefs[-1] != 0)
  
  #================================global undersmoothing================================
  start <- Sys.time()
  
  CV_nonzero_col <- which(hal_CV$coefs[-1] != 0)
  # if all coefs are zero, skip undersmooth and use the initial fit
  if (length(CV_nonzero_col) == 0){
    lambda_u_g = lambda_CV
    
    print(sprintf('  globally u lambdas: %f (same as CV-HAL)', lambda_u_g))
    hal_u_g <- hal_CV
  }else{
    CV_basis_mat <- as.matrix(hal_CV$x_basis)
    CV_basis_mat <- as.matrix(CV_basis_mat[, CV_nonzero_col])
    
    hal_undersmooth <- undersmooth_hal(X, Y,
                                       fit_init = hal_CV,
                                       family = y_type)
    lambda_u_g = hal_undersmooth$lambda_under
    
    print(sprintf('  globally u lambdas: %f', lambda_u_g))
    if(is.na(lambda_u_g)){
      lambda_u_g <- lambda_CV
      print(sprintf('  globally u lambdas: %f', lambda_u_g))
      
      hal_u_g <- hal_CV
    } else {
      
      hal_u_g <- fit_hal(X = X, Y = Y, family = y_type,
                         smoothness_orders = 0,
                         return_x_basis = TRUE,
                         fit_control = list(
                           cv_select = FALSE,
                           n_folds = 10,
                           foldid = NULL,
                           use_min = TRUE,
                           lambda.min.ratio = 1e-4,
                           prediction_bounds = "default"
                         ),
                         lambda = lambda_u_g)
    }
  }
  
  end <- Sys.time()
  
  hal_fit_time = end - start
  if(units(hal_fit_time) == 'secs'){
    hal_fit_time = as.numeric(hal_fit_time)
  } else if (units(hal_fit_time) == 'mins'){
    hal_fit_time = as.numeric(hal_fit_time)  * 60
  } else if (units(hal_fit_time) == 'hours'){
    hal_fit_time = as.numeric(hal_fit_time)  * 60 * 24
  }
  
  hal_u_g_fit_time = hal_fit_time + hal_cv_fit_time
  
  num_basis_u_g <- sum(hal_u_g$coefs[-1] != 0)
  
  #================================return================================
  hal_fit_list <- list(hal_CV = hal_CV, hal_u_g = hal_u_g)
  lambda_list <- list(lambda_CV = lambda_CV, lambda_u_g = lambda_u_g)
  hal_fit_time_list <- list(hal_cv_fit_time = hal_cv_fit_time, hal_u_g_fit_time = hal_u_g_fit_time)
  num_basis_list <- list(num_basis_CV = num_basis_CV, num_basis_u_g = num_basis_u_g)
  
  return(list(hal_fit_list = hal_fit_list, lambda_list = lambda_list, hal_fit_time_list = hal_fit_time_list, num_basis_list = num_basis_list))
}



run_simu_smoothness_adaptive_HAL_1round <- function(gen_data_func, eval_points, y_type, n){
  
  obs <- gen_data_func(n)

  y_name = "Y"
  x_names = names(obs)[names(obs) != 'Y']
  y_type = "binomial"
  
  Y <- as.numeric(as.matrix(obs %>% select(all_of(y_name))))
  X <- obs %>% 
    select(all_of(x_names)) %>% 
    mutate_if(sapply(., is.factor), as.numeric)
  
  # 1. --- SL
  task <- make_sl3_Task(
    data = obs,
    outcome = y_name,
    covariates = x_names
  )
  
  # num_knots = c(200, 100,  50)
  lrn_hal0 <- Lrnr_hal9001$new(smoothness_orders = 0, family = y_type, return_x_basis = TRUE)
  lrn_hal1 <- Lrnr_hal9001$new(smoothness_orders = 1, family = y_type, return_x_basis = TRUE)
  lrn_hal2 <- Lrnr_hal9001$new(smoothness_orders = 2, family = y_type, return_x_basis = TRUE)
  lrn_hal3 <- Lrnr_hal9001$new(smoothness_orders = 3, family = y_type, return_x_basis = TRUE)
  
  # num_knots = c(20,10,5)
  lrn_hal0_smallerKnots <- Lrnr_hal9001$new(smoothness_orders = 0,
                                            family = y_type, 
                                            return_x_basis = TRUE, 
                                            num_knots = hal9001:::num_knots_generator(
                                              max_degree = ifelse(ncol(X) >= 20, 2, 3), 
                                              smoothness_orders = 0,
                                              base_num_knots_0 = 20,
                                              base_num_knots_1 = 20  
                                            ))
  lrn_hal1_smallerKnots <- Lrnr_hal9001$new(smoothness_orders = 1,
                                            family = y_type, 
                                            return_x_basis = TRUE, 
                                            num_knots = hal9001:::num_knots_generator(
                                              max_degree = ifelse(ncol(X) >= 20, 2, 3), 
                                              smoothness_orders = 1,
                                              base_num_knots_0 = 20,
                                              base_num_knots_1 = 20  
                                            ))
  lrn_hal2_smallerKnots <- Lrnr_hal9001$new(smoothness_orders = 2,
                                            family = y_type, 
                                            return_x_basis = TRUE, 
                                            num_knots = hal9001:::num_knots_generator(
                                              max_degree = ifelse(ncol(X) >= 20, 2, 3), 
                                              smoothness_orders = 2,
                                              base_num_knots_0 = 20,
                                              base_num_knots_1 = 20  
                                            ))
  lrn_hal3_smallerKnots <- Lrnr_hal9001$new(smoothness_orders = 3,
                                            family = y_type, 
                                            return_x_basis = TRUE, 
                                            num_knots = hal9001:::num_knots_generator(
                                              max_degree = ifelse(ncol(X) >= 20, 2, 3), 
                                              smoothness_orders = 3,
                                              base_num_knots_0 = 20,
                                              base_num_knots_1 = 20  
                                            ))
  
  learners <- c(lrn_hal0, lrn_hal0_smallerKnots, lrn_hal1, lrn_hal1_smallerKnots, 
                lrn_hal2, lrn_hal2_smallerKnots, lrn_hal3, lrn_hal3_smallerKnots)
  names(learners) <- c("halcv0", "halcv0_sKnots", "halcv1", "halcv1_sKnots", 
                       "halcv2", "halcv2_sKnots", "halcv3", "halcv3_sKnots")
  stack <- make_learner(Stack, learners)
  
  cv_selector <- Lrnr_cv_selector$new(eval_function = loss_loglik_binomial) # https://tlverse.org/sl3/reference/loss_functions.html
  dSL <- Lrnr_sl$new(learners = stack, metalearner = cv_selector)
  
  dSL_fit <- dSL$train(task)
  
  #--- 2. SL fit results
  Dsl_fit_summary <- dSL_fit$print()
  smooth_orders = c(0,0,1,1,2,2,3,3)
  num_knots = rep(c("default", "smaller"),4)
  
  results <- list()
  
  for (i in 1:length(dSL_fit$learner_fits)) {
    
    hal_CV <- dSL_fit$learner_fits[[i]]$fit_object
    
    # risk
    cv_risk = Dsl_fit_summary[i,3]
    
    # smooth order & if default number of knots of HAL fit
    SO = smooth_orders[i]
    n_knots_default = as.numeric(num_knots[i] == "default")
    
    # sl_pick 
    sl_pick = as.numeric(dSL_fit$coefficients[i]==1)
    
    # lambda
    lambda_CV <- hal_CV$lambda_star
    
    # number of basis
    num_basis_CV <- sum(hal_CV$coefs[-1] != 0)
    
    # prediction
    psi_hat <- sapply(eval_points, function(a){ X_new <- X
    X_new$A = a
    mean(predict(hal_CV, new_data = X_new)) } )
    
    # IC-based inference
    psi_hat_pnt_cv_se <- IC_based_se(X, Y, hal_CV, eval_points)
    
    # returns
    psi_hat_pnt_cv <- cbind(eval_points, matrix(psi_hat, ncol=1), lambda_CV,  1, 
                            psi_hat_pnt_cv_se, 
                            num_basis_CV, cv_risk, SO, n_knots_default, sl_pick)
    
    colnames(psi_hat_pnt_cv) <- c("a", "y_hat", "lambda", "lambda_scaler", "SE",
                                  "n_basis", "cv_risk", "smooth_order", 
                                  "if_n_knots_default", "sl_pick")
    
    psi_hat_pnt_cv <- as.data.frame(psi_hat_pnt_cv) %>% 
      mutate(ci_lwr = y_hat - 1.96 * SE,
             ci_upr = y_hat + 1.96 * SE)
    
    if(y_type == "binomial"){
      bounds <- c(0, 1)
    } else {
      bounds <- c(min(Y), max(Y))
    }
    psi_hat_pnt_cv[,"y_hat"] <- pmax(bounds[1], psi_hat_pnt_cv[,"y_hat"])
    psi_hat_pnt_cv[,"y_hat"] <- pmin(psi_hat_pnt_cv[,"y_hat"], bounds[2])
    
    psi_hat_pnt_cv[,"ci_lwr"] <- pmax(bounds[1], psi_hat_pnt_cv[,"ci_lwr"])
    psi_hat_pnt_cv[,"ci_lwr"] <- pmin(psi_hat_pnt_cv[,"ci_lwr"], bounds[2])
    
    psi_hat_pnt_cv[,"ci_upr"] = pmax(bounds[1], psi_hat_pnt_cv[,"ci_upr"])
    psi_hat_pnt_cv[,"ci_upr"] <- pmin(psi_hat_pnt_cv[,"ci_upr"], bounds[2])
    
    results[[i]] <- psi_hat_pnt_cv
  }
  
  names(results) = names(learners)
  
  #--- 3. SL pick undersmooth
  # undersmo0th hal based on the one sl picked
  sl_pick_idx = which(dSL_fit$coefficients==1)
  
  hal_fit_sl_pick = dSL_fit$learner_fits[[sl_pick_idx]]$fit_object
  SO_pick = smooth_orders[sl_pick_idx]
  n_knots_default_pick = as.numeric(num_knots[sl_pick_idx] == "default")
  
  hal_undersmooth <- undersmooth_hal(X, Y,
                                     fit_init = hal_fit_sl_pick,
                                     family = y_type)
  
  lambda_u_g = hal_undersmooth$lambda_under
  
  if(is.na(lambda_u_g) | hal_undersmooth$lambda_under == hal_undersmooth$lambda_init){
    lambda_u_g <- lambda_CV
    print(sprintf('  globally u lambdas: %f', lambda_u_g))
    
    hal_u_g <- hal_fit_sl_pick
  } else {
    if(n_knots_default_pick == 1){
      hal_u_g <- fit_hal(X = X, Y = Y, 
                         family = y_type,
                         smoothness_orders = SO_pick,
                         return_x_basis = TRUE,
                         fit_control = list(cv_select = FALSE),
                         lambda = lambda_u_g)
    } else {
      hal_u_g <- fit_hal(X = X, Y = Y, 
                         family = y_type,
                         smoothness_orders = SO_pick,
                         return_x_basis = TRUE,
                         fit_control = list(cv_select = FALSE),
                         lambda = lambda_u_g,
                         num_knots = hal9001:::num_knots_generator(
                           max_degree = ifelse(ncol(X) >= 20, 2, 3), 
                           smoothness_orders = SO_pick,
                           base_num_knots_0 = 20,
                           base_num_knots_1 = 20  
                         ))
    }
    
  }
  
  # number of basis
  num_basis_u_g <- sum(hal_u_g$coefs[-1] != 0)
  
  # prediction
  psi_hat_u_g <- sapply(eval_points, function(a){ X_new <- X
  X_new$A = a
  mean(predict(hal_u_g, new_data = X_new)) } )
  
  # IC-based inference
  psi_hat_pnt_u_g_se <- IC_based_se(X, Y, hal_u_g, eval_points)
  
  # returns
  psi_hat_pnt_u_g <- cbind(eval_points, matrix(psi_hat_u_g, ncol=1), lambda_u_g,  lambda_u_g / hal_undersmooth$lambda_init, 
                           psi_hat_pnt_u_g_se, 
                           num_basis_u_g, NA, SO_pick, n_knots_default_pick, 1)
  
  colnames(psi_hat_pnt_u_g) <- c("a", "y_hat", "lambda", "lambda_scaler", "SE",
                                 "n_basis", "cv_risk", "smooth_order", 
                                 "if_n_knots_default", "sl_pick")
  
  psi_hat_pnt_u_g <- as.data.frame(psi_hat_pnt_u_g) %>% 
    mutate(ci_lwr = y_hat - 1.96 * SE,
           ci_upr = y_hat + 1.96 * SE)
  
  if(y_type == "binomial"){
    bounds <- c(0, 1)
  } else {
    bounds <- c(min(Y), max(Y))
  }
  psi_hat_pnt_u_g[,"y_hat"] <- pmax(bounds[1], psi_hat_pnt_u_g[,"y_hat"])
  psi_hat_pnt_u_g[,"y_hat"] <- pmin(psi_hat_pnt_u_g[,"y_hat"], bounds[2])
  
  psi_hat_pnt_u_g[,"ci_lwr"] <- pmax(bounds[1], psi_hat_pnt_u_g[,"ci_lwr"])
  psi_hat_pnt_u_g[,"ci_lwr"] <- pmin(psi_hat_pnt_u_g[,"ci_lwr"], bounds[2])
  
  psi_hat_pnt_u_g[,"ci_upr"] = pmax(bounds[1], psi_hat_pnt_u_g[,"ci_upr"])
  psi_hat_pnt_u_g[,"ci_upr"] <- pmin(psi_hat_pnt_u_g[,"ci_upr"], bounds[2])
  
  results[[length(dSL_fit$learner_fits) + 1]] <- psi_hat_pnt_u_g
  names(results)[length(dSL_fit$learner_fits) + 1] = "U_HAL_SL_pick"
  
  return(results)
}

run_simu_smoothness_adaptive_HAL_rep <- function(gen_data_func, eval_points, y_type, n, rounds, return_all_rslts = F){
  result_list <- list()
  
  for(r in 1:rounds){
    print(paste0("round ", r))
    result <- tryCatch({
      run_simu_smoothness_adaptive_HAL_1round(gen_data_func, eval_points, y_type, n=n)
    }, error = function(e) {
      print(paste0("Error: ", e$message))
      NULL
    })
    
    while(is.null(result)) {
      print('retry with a new generated data')
      result <- tryCatch({
        run_simu_smoothness_adaptive_HAL_1round(gen_data_func, eval_points, y_type, n=n)
      }, error = function(e) {
        print(paste0("Error: ", e$message))
        NULL
      })
    }
    
    result_list[[r]] <- result
  }
  
  result_summaries <- list()
  methods = names(result_list[[1]])
  
  for (method in methods){
    result_list_method <- lapply(result_list, function(lst) lst[[method]])
    result_all <-  do.call("rbind", result_list_method) %>% as.data.frame()
    result_all <- merge(as.data.frame(psi0_pnt), result_all, by=c("a"))
    
    check_n_basis <- function(df) {
      !all(sapply(df$n_basis, is.numeric) & sapply(df$n_basis, length) == 1)
      }
    data_frames_with_varying_n_basis <- result_list_method[sapply(result_list_method, check_n_basis)]
    
    result_summary <- result_all %>% 
      filter(SE != 0) %>% 
      mutate(bias = abs(y_hat - psi0),
             bias_se_ratio = bias / SE,
             cover_rate = as.numeric(ci_lwr <= psi0 & psi0 <= ci_upr)) %>% 
      group_by(a) %>% 
      mutate(oracal_SE = sqrt(var(y_hat)),
             oracal_bias_se_ratio = bias / oracal_SE,
             oracal_ci_lwr = y_hat - 1.96 * oracal_SE,
             oracal_ci_upr = y_hat + 1.96 * oracal_SE,
             oracal_cover_rate = as.numeric(oracal_ci_lwr <= psi0 & psi0 <= oracal_ci_upr)) %>%
      summarise(across(where(is.numeric), mean)) %>% 
      ungroup() 
    
    result_summary$method = method
    result_summary$n_basis = mean(result_all$n_basis)
    result_summary$smooth_order = mean(result_all$smooth_order)
    result_summary$if_n_knots_default = mean(result_all$if_n_knots_default)
    
    result_summaries[[method]] = result_summary
  }
  
  result_summary = do.call("rbind", result_summaries) %>% as.data.frame()
  
  results <- list(result_summary = result_summary, result_list = result_list)
  
  return(results)
}

summary_smoothness_adaptive_HAL <- function(result_list){
  result_summaries <- list()
  methods = names(result_list[[1]])
  
  for (method in methods){
    result_list_method <- lapply(result_list, function(lst) lst[[method]])
    result_all <-  do.call("rbind", result_list_method) %>% as.data.frame()
    result_all <- merge(as.data.frame(psi0_pnt), result_all, by=c("a"))
    
    check_n_basis <- function(df) {
      !all(sapply(df$n_basis, is.numeric) & sapply(df$n_basis, length) == 1)
    }
    data_frames_with_varying_n_basis <- result_list_method[sapply(result_list_method, check_n_basis)]
    
    result_summary <- result_all %>% 
      filter(SE != 0) %>% 
      mutate(bias = abs(y_hat - psi0),
             bias_se_ratio = bias / SE,
             cover_rate = as.numeric(ci_lwr <= psi0 & psi0 <= ci_upr)) %>% 
      group_by(a) %>% 
      mutate(oracal_SE = sqrt(var(y_hat)),
             oracal_bias_se_ratio = bias / oracal_SE,
             oracal_ci_lwr = y_hat - 1.96 * oracal_SE,
             oracal_ci_upr = y_hat + 1.96 * oracal_SE,
             oracal_cover_rate = as.numeric(oracal_ci_lwr <= psi0 & psi0 <= oracal_ci_upr)) %>%
      summarise(across(where(is.numeric), mean)) %>% 
      ungroup() 
    
    result_summary$method = method
    result_summary$n_basis = mean(result_all$n_basis)
    result_summary$smooth_order = mean(result_all$smooth_order)
    result_summary$if_n_knots_default = mean(result_all$if_n_knots_default)
    
    result_summaries[[method]] = result_summary
  }
  
  result_summary = do.call("rbind", result_summaries) %>% as.data.frame()
  
  results <- list(result_summary = result_summary, result_list = result_list)
}

##############################################################
# This function runs the simulation with given:
#       - data generating function
#       - sample size n
# returns the estimated ATE and empirical 95% CI
##############################################################
run_simu_1round <- function(gen_data_func, eval_points, y_type, n, defualt_setting = F){
  
  bootstrap = F
  
  obs <- gen_data_func(n)
  
  y_name = "Y"
  x_names = names(obs)[names(obs) != 'Y']
  y_type = "binomial"
  
  Y <- as.numeric(as.matrix(obs %>% select(all_of(y_name))))
  X <- obs %>% 
    select(all_of(x_names)) %>% 
    mutate_if(sapply(., is.factor), as.numeric)
  
  if(defualt_setting){
    # print("d")
    fit_hal_all_criteria_rslts <- fit_hal_CV_U_0(X, Y, y_type, eval_points)
  } else {
    fit_hal_all_criteria_rslts <- fit_hal_CV_U(X, Y, y_type, eval_points)
  }
  
  
  #================================CV-HAL================================
  lambda_CV <- fit_hal_all_criteria_rslts$lambda_list$lambda_CV
  psi_hat <- sapply(eval_points, function(a){ X_new <- X
  X_new$A = a
  mean(predict(fit_hal_all_criteria_rslts$hal_fit_list$hal_CV, new_data = X_new)) } )
  
  psi_hat_pnt_cv <- cbind(eval_points, matrix(psi_hat, ncol=1), lambda_CV, 1, 
                          fit_hal_all_criteria_rslts$hal_fit_time_list$hal_cv_fit_time,
                          fit_hal_all_criteria_rslts$num_basis_list$num_basis_CV)
  colnames(psi_hat_pnt_cv) <- c("a", "y_hat", "lambda", "lambda_scaler", "hal_fit_time", "n_basis")
  
  
  # IC-based inference
  psi_hat_pnt_cv_se <- IC_based_se(X, Y, fit_hal_all_criteria_rslts$hal_fit_list$hal_CV, eval_points)
  psi_hat_pnt_cv <- as.data.frame(psi_hat_pnt_cv) %>% 
    mutate(SE = psi_hat_pnt_cv_se,
           ci_lwr = y_hat - 1.96 * SE,
           ci_upr = y_hat + 1.96 * SE)
  
  if(y_type == "binomial"){
    bounds <- c(0, 1)
  } else {
    bounds <- c(min(Y), max(Y))
  }
  psi_hat_pnt_cv[,"y_hat"] <- pmax(bounds[1], psi_hat_pnt_cv[,"y_hat"])
  psi_hat_pnt_cv[,"y_hat"] <- pmin(psi_hat_pnt_cv[,"y_hat"], bounds[2])
  
  psi_hat_pnt_cv[,"ci_lwr"] <- pmax(bounds[1], psi_hat_pnt_cv[,"ci_lwr"])
  psi_hat_pnt_cv[,"ci_lwr"] <- pmin(psi_hat_pnt_cv[,"ci_lwr"], bounds[2])
  
  psi_hat_pnt_cv[,"ci_upr"] = pmax(bounds[1], psi_hat_pnt_cv[,"ci_upr"])
  psi_hat_pnt_cv[,"ci_upr"] <- pmin(psi_hat_pnt_cv[,"ci_upr"], bounds[2])
  
  # bootstrap-based inference
  if(bootstrap){
    psi_hat_pnt_cv_bt_bds <- bootstrap_inference(X, Y, eval_points, fit_hal_all_criteria_rslts$hal_fit_list$hal_CV, y_type)
    psi_hat_pnt_cv$ci_lwr_bt <- psi_hat_pnt_cv_bt_bds$lower_bd
    psi_hat_pnt_cv$ci_upr_bt <- psi_hat_pnt_cv_bt_bds$upper_bd
    psi_hat_pnt_cv$SE_bt <- psi_hat_pnt_cv_bt_bds$SE
  }
  
  #================================global undersmoothing================================
  if(any(fit_hal_all_criteria_rslts$hal_fit_list$hal_u_g$coefs[-1] != 0)){
    
    psi_hat <- sapply(eval_points, function(a){ X_new <- X
    X_new$A = a
    mean(predict(fit_hal_all_criteria_rslts$hal_fit_list$hal_u_g, new_data = X_new)) } )
    
    lambda_scaler = fit_hal_all_criteria_rslts$lambda_list$lambda_u_g / lambda_CV
    psi_hat_pnt_u_g <- cbind(eval_points, matrix(psi_hat, ncol=1), fit_hal_all_criteria_rslts$lambda_list$lambda_u_g, lambda_scaler, 
                             fit_hal_all_criteria_rslts$hal_fit_time_list$hal_u_g_fit_time,
                             fit_hal_all_criteria_rslts$num_basis_list$num_basis_u_g)
    
    colnames(psi_hat_pnt_u_g) <- c("a", "y_hat", "lambda", "lambda_scaler", "hal_fit_time", "n_basis")
    
    # IC-based inference
    psi_hat_pnt_u_g_se <- IC_based_se(X, Y, fit_hal_all_criteria_rslts$hal_fit_list$hal_u_g, eval_points)
    psi_hat_pnt_u_g <- as.data.frame(psi_hat_pnt_u_g) %>% 
      mutate(SE = psi_hat_pnt_u_g_se,
             ci_lwr = y_hat - 1.96 * SE,
             ci_upr = y_hat + 1.96 * SE)
    
    if(y_type == "binomial"){
      bounds <- c(0, 1)
    } else {
      bounds <- c(min(Y), max(Y))
    }
    psi_hat_pnt_u_g[,"y_hat"] <- pmax(bounds[1], psi_hat_pnt_u_g[,"y_hat"])
    psi_hat_pnt_u_g[,"y_hat"] <- pmin(psi_hat_pnt_u_g[,"y_hat"], bounds[2])
    
    psi_hat_pnt_u_g[,"ci_lwr"] <- pmax(bounds[1], psi_hat_pnt_u_g[,"ci_lwr"])
    psi_hat_pnt_u_g[,"ci_lwr"] <- pmin(psi_hat_pnt_u_g[,"ci_lwr"], bounds[2])
    
    psi_hat_pnt_u_g[,"ci_upr"] = pmax(bounds[1], psi_hat_pnt_u_g[,"ci_upr"])
    psi_hat_pnt_u_g[,"ci_upr"] <- pmin(psi_hat_pnt_u_g[,"ci_upr"], bounds[2])
    
    
    if(bootstrap){
      # bootstrap-based inference
      psi_hat_pnt_u_g_bt_bds <- bootstrap_inference(X, Y, eval_points, fit_hal_all_criteria_rslts$hal_fit_list$hal_u_g, y_type)
      psi_hat_pnt_u_g$ci_lwr_bt <- psi_hat_pnt_u_g_bt_bds$lower_bd
      psi_hat_pnt_u_g$ci_upr_bt <- psi_hat_pnt_u_g_bt_bds$upper_bd
      psi_hat_pnt_u_g$SE_bt <- psi_hat_pnt_u_g_bt_bds$SE
    }
    
    
  } else {
    psi_hat_pnt_u_g <- NA
  }
  
  #====================================================================================
  results <- list(psi_hat_pnt_cv, psi_hat_pnt_u_g)
  names(results) = c("CV", "U_G")
  return(results)
}


##############################################################
# This function runs the simulation for B rounds with given:
#       - data generating function
#       - sample size n
#       - number of simulations: B

# returns the estimated ATE, empirical 95% CI, and coverage rates
##############################################################
run_simu_rep <- function(gen_data_func, eval_points, y_type, n, rounds, return_all_rslts = F, defualt_setting = F){
  
  bootstrap = F
  
  result_list <- list()
  
  for(r in 1:rounds){
    print(paste0("round ", r))
    result <- tryCatch({
      run_simu_1round(gen_data_func, eval_points, y_type, n=n, defualt_setting)
    }, error = function(e) {
      print(paste0("Error: ", e$message))
      NULL
    })
    
    while(is.null(result)) {
      print('retry with a new generated data')
      result <- tryCatch({
        run_simu_1round(gen_data_func, eval_points, y_type, n=n, defualt_setting)
      }, error = function(e) {
        print(paste0("Error: ", e$message))
        NULL
      })
    }
    
    result_list[[r]] <- result
  }
  
  results <- list()
  for (method in c("CV", "U_G")){
    result_list_method <- lapply(result_list, function(lst) lst[[method]])
    result_all <-  do.call("rbind", result_list_method) %>% as.data.frame()
    result_all <- merge(as.data.frame(psi0_pnt), result_all, by=c("a"))
    
    if(bootstrap){
      result_summary <- result_all %>% 
        filter(SE != 0) %>% 
        mutate(bias = abs(y_hat - psi0),
               bias_se_ratio = bias / SE,
               # bias_se_ratio_bt = bias / SE_bt,
               # cover_rate_bt = as.numeric(ci_lwr_bt <= psi0 & psi0 <= ci_upr_bt),
               cover_rate = as.numeric(ci_lwr <= psi0 & psi0 <= ci_upr)) %>% 
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
    } else {
      result_summary <- result_all %>% 
        filter(SE != 0) %>% 
        mutate(bias = abs(y_hat - psi0),
               bias_se_ratio = bias / SE,
               cover_rate = as.numeric(ci_lwr <= psi0 & psi0 <= ci_upr)) %>% 
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
    }
    
    
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
run_simu_1round_scalers <- function(gen_data_func, eval_points, y_type, n, lambda_scalers){
  
  obs <- gen_data_func(n)

  y_name = "Y"
  x_names = names(obs)[names(obs) != 'Y']
  
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
    
    num_basis <- sum(hal_fit$coefs[-1] != 0)
    
    #--------------------------estimations---------------------------------
    coef <- hal_fit$coefs
    basis_mat <- cbind(1, as.matrix(hal_fit$x_basis))
    
    nonzero_idx <- which(coef != 0)
    
    if(length(nonzero_idx) > 0) {
      
      psi_hat <- sapply(eval_points, function(a){ X_new <- X
                                                  X_new$A = a
                                                  mean(predict(hal_fit, new_data = X_new)) } )
      
      psi_hat_pnt_scaled <- cbind(eval_points, matrix(psi_hat, ncol=1), lambda_CV, lambda_scaler, hal_scaled_fit_time, num_basis)
      colnames(psi_hat_pnt_scaled) <- c("a", "y_hat", "lambda", "lambda_scaler", "hal_fit_time", "n_basis")
      
      # IC-based inference
      psi_hat_pnt_scaled_se <- IC_based_se(X, Y, hal_fit, eval_points)
      psi_hat_pnt_scaled <- as.data.frame(psi_hat_pnt_scaled) %>% 
        mutate(SE = psi_hat_pnt_scaled_se,
               ci_lwr = y_hat - 1.96 * SE,
               ci_upr = y_hat + 1.96 * SE)
      
      if(y_type == "binomial"){
        bounds <- c(0, 1)
      } else {
        bounds <- c(min(Y), max(Y))
      }
      psi_hat_pnt_scaled[,"y_hat"] <- pmax(bounds[1], psi_hat_pnt_scaled[,"y_hat"])
      psi_hat_pnt_scaled[,"y_hat"] <- pmin(psi_hat_pnt_scaled[,"y_hat"], bounds[2])
      
      psi_hat_pnt_scaled[,"ci_lwr"] <- pmax(bounds[1], psi_hat_pnt_scaled[,"ci_lwr"])
      psi_hat_pnt_scaled[,"ci_lwr"] <- pmin(psi_hat_pnt_scaled[,"ci_lwr"], bounds[2])
      
      psi_hat_pnt_scaled[,"ci_upr"] = pmax(bounds[1], psi_hat_pnt_scaled[,"ci_upr"])
      psi_hat_pnt_scaled[,"ci_upr"] <- pmin(psi_hat_pnt_scaled[,"ci_upr"], bounds[2])
      
      # bootstrap-based inference
      # psi_hat_pnt_scaled_bt_bds <- bootstrap_inference(X, Y, eval_points, hal_fit, y_type)
      # psi_hat_pnt_scaled$ci_lwr_bt <- psi_hat_pnt_scaled_bt_bds$lower_bd
      # psi_hat_pnt_scaled$ci_upr_bt <- psi_hat_pnt_scaled_bt_bds$upper_bd
      # psi_hat_pnt_scaled$SE_bt <- psi_hat_pnt_scaled_bt_bds$SE
      
    } else {
      psi_hat_pnt_scaled = as.data.frame(cbind(eval_points, NA, lambda_CV, lambda_scaler, hal_scaled_fit_time, NA, NA, NA, NA, NA, NA, NA))
      names(psi_hat_pnt_scaled) = c("a", "y_hat", "lambda", "lambda_scaler", "hal_fit_time", "n_basis", "SE", "ci_lwr", "ci_upr", "ci_lwr_bt", "ci_upr_bt", "SE_bt")
    }
    results[[i]] = psi_hat_pnt_scaled
  }
  
  names(results) = paste0("scale=", round(lambda_scalers, 4))
  
  return(results)
}


##############################################################

run_simu_scaled_rep <- function(gen_data_func, eval_points, y_type, n, rounds, return_all_rslts = F){
  lambda_scalers = c(1.2, 1.1, 10^seq(from=0, to=-3, length=20))
  result_list <- list()
  for(r in 1:rounds){
    print(paste0("round ", r))
    result <- tryCatch({
      run_simu_1round_scalers(gen_data_func, eval_points, y_type, n=n, lambda_scalers=lambda_scalers)
    }, error = function(e) {
      print(paste0("Error: ", e$message))
      NULL
    })
    
    while(is.null(result)) {
      print('retry with a new generated data')
      result <- tryCatch({
        run_simu_1round_scalers(gen_data_func, eval_points, y_type, n=n, lambda_scalers=lambda_scalers)
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

    result_list_scale <- lapply(result_list, function(lst) lst[[i]][c('a', 'y_hat', 'lambda', 'lambda_scaler', 'hal_fit_time', 'n_basis', 'SE', 'ci_lwr', 'ci_upr')])
    no_empirical_CI_proportion[i] <- mean(sapply(result_list_scale, function(rlt) any(is.na(rlt[,colnames(rlt) == 'SE']))))
    result_all <-  do.call("rbind", result_list_scale) %>% as.data.frame()
    result_all <- merge(as.data.frame(psi0_pnt), result_all, by=c("a"))
    
    result_summary <- result_all %>%
      filter((SE != 0) | (is.na(SE))) %>%
      mutate(bias = abs(y_hat - psi0),
             bias_se_ratio = bias / SE,
             # bias_se_ratio_bt = bias / SE_bt,
             # cover_rate_bt = as.numeric(ci_lwr_bt <= psi0 & psi0 <= ci_upr_bt) ,
             cover_rate = as.numeric(ci_lwr <= psi0 & psi0 <= ci_upr)) %>%
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
  for (i in 2:length(lambda_scalers)) {
    result_summary <- rbind(result_summary, results[[i]]$result_summary)
  }
  results$result_summary <- result_summary

  results$no_empirical_CI_proportion <- no_empirical_CI_proportion

  return(results)
}

