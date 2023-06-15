
fit_hal_all_criteria <- function(X, Y, y_type){
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
  
  hal_fit_time = end - start
  if(units(hal_fit_time) == 'secs'){
    hal_fit_time = as.numeric(hal_fit_time)
  } else if (units(hal_fit_time) == 'mins'){
    hal_fit_time = as.numeric(hal_fit_time)  * 60
  } else if (units(hal_fit_time) == 'hours'){
    hal_fit_time = as.numeric(hal_fit_time)  * 60 * 24
  }
  
  hal_u_g_fit_time = hal_fit_time + hal_cv_fit_time
  
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
  
  hal_fit_time = end - start
  if(units(hal_fit_time) == 'secs'){
    hal_fit_time = as.numeric(hal_fit_time)
  } else if (units(hal_fit_time) == 'mins'){
    hal_fit_time = as.numeric(hal_fit_time)  * 60
  } else if (units(hal_fit_time) == 'hours'){
    hal_fit_time = as.numeric(hal_fit_time)  * 60 * 24
  }
  
  hal_u_l_fit_time = hal_fit_time + hal_cv_fit_time
  
  #================================return================================
  hal_fit_list <- list(hal_CV = hal_CV, hal_u_g = hal_u_g, hal_u_l = hal_u_l)
  lambda_list <- list(lambda_CV = lambda_CV, lambda_u_g = lambda_u_g, lambda_u_l = lambda_u_l)
  hal_fit_time_list <- list(hal_cv_fit_time = hal_cv_fit_time, hal_u_g_fit_time = hal_u_g_fit_time, hal_u_l_fit_time = hal_u_l_fit_time)
  
  return(list(hal_fit_list = hal_fit_list, lambda_list = lambda_list, hal_fit_time_list = hal_fit_time_list, lambda_u_l_idx = lambda_u_l_idx))
}


IC_based_se <- function(hal_fit, eval_points){
  
  coef <- hal_fit$coefs
  basis_mat <- cbind(1, as.matrix(hal_fit$x_basis))
  
  nonzero_idx <- which(coef != 0)
  
  if(length(nonzero_idx) > 0) {
    coef_nonzero <- coef[nonzero_idx]
    basis_mat_nonzero <- as.matrix(basis_mat[, nonzero_idx])
    
    IC_beta <- cal_IC_for_beta(X = basis_mat_nonzero, 
                               Y = Y, 
                               Y_hat = predict(hal_fit, new_data = X, type = "response"),
                               beta_n = coef_nonzero
    )
    
    se <- c()
    
    for (i in 1:length(eval_points)) {
      X_new <- X
      X_new$A = eval_points[i]
      
      # efficient influence curve
      x_basis_a <- make_design_matrix(as.matrix(X_new), hal_fit$basis_list, p_reserve = 0.75)
      x_basis_a_nonzero <- as.matrix(cbind(1, x_basis_a)[, nonzero_idx])
      
      IC_EY <- cal_IC_for_EY(X_new = x_basis_a_nonzero, 
                             beta_n = coef_nonzero, IC_beta = IC_beta)
      
      # empirical SE
      se[i] <- sqrt(var(IC_EY)/n)
      #ci_lwr <- Y_hat - 1.96 * SE
      #ci_upr <- Y_hat + 1.96 * SE

    }
  } else {
    se <- NA
  }

  return(se)
}

IC_based_se_u_l <- function(hal_fit, single_eval_point){
  
  coef <- hal_fit$coefs
  basis_mat <- cbind(1, as.matrix(hal_fit$x_basis))
  
  nonzero_idx <- which(coef != 0)
  
  if(length(nonzero_idx) > 0) {
    coef_nonzero <- coef[nonzero_idx]
    basis_mat_nonzero <- as.matrix(basis_mat[, nonzero_idx])
    
    IC_beta <- cal_IC_for_beta(X = basis_mat_nonzero, 
                               Y = Y, 
                               Y_hat = predict(hal_fit, new_data = X, type = "response"),
                               beta_n = coef_nonzero
    )
    
    X_new <- X
    X_new$A = eval_points[i]
      
    # efficient influence curve
    x_basis_a <- make_design_matrix(as.matrix(X_new), hal_fit$basis_list, p_reserve = 0.75)
    x_basis_a_nonzero <- as.matrix(cbind(1, x_basis_a)[, nonzero_idx])
      
    IC_EY <- cal_IC_for_EY(X_new = x_basis_a_nonzero, 
                           beta_n = coef_nonzero, IC_beta = IC_beta)
      
    # empirical SE
    se <- sqrt(var(IC_EY)/n)
      
  } else {
    se <- NA
  }
  
  return(se)
}

bootstrap_inference <- function(X, Y, eval_points, hal_fit, y_type, B = 200){
  
  basis_list <- hal_fit$basis_list

  y_hat_B <- matrix(ncol = length(eval_points))
  
  for (b in 1:B) {
    #--------------data--------------
    idx <- sample(1:length(Y), length(Y), replace = T)
    Xb <- X[idx,]
    Yb <- Y[idx]
    
    if (length(basis_list) > 0) {
      # generate basis matrix
      x_basis <- hal9001::make_design_matrix(as.matrix(Xb), basis_list)
      
      #--------------fit--------------
      lasso_fit <- tryCatch({
        lasso_fit <- glmnet::glmnet(x = x_basis, y = Yb, 
                                    family = y_type, 
                                    lambda = hal_fit$lambda_star,
                                    intercept = FALSE, standardize = FALSE)
      },
      error = function(){
        lasso_fit <- NA
      })
      
      y_hat_b <- c()
      for (i in 1:length(eval_points)) {
        X_new <- Xb
        X_new$A = eval_points[i]
  
        x_basis_a <- hal9001::make_design_matrix(as.matrix(X_new), hal_fit$basis_list)
        
        #--------------prediction--------------
        if (y_type == "binomial") {
          preds <- predict(lasso_fit, x_basis_a, type = "response")
        } else {
          preds <- predict(lasso_fit, x_basis_a)
        }
        
        y_hat_b[i] = mean(preds)
      }
    } else {
      y_hat_b = rep(NA, length(eval_points))
    }

    y_hat_B <- rbind(y_hat_B, y_hat_b)
  }
  
  #--------------confidence bounds--------------
  y_hat_B <- y_hat_B[-1,]
  lower_bd <- apply(y_hat_B, 2, quantile, probs = 0.05/2, na.rm = T)
  upper_bd <- apply(y_hat_B, 2, quantile, probs = 1-0.05/2, na.rm = T)
  
  return(list(lower_bd=lower_bd, upper_bd=upper_bd))
  #   
  #   basis_list <- hal_fit$basis_list
  #   # copy_map <- hal_fit$copy_map
  #   
  #   y_hat_B <- matrix(ncol = length(eval_points))
  #   for (b in 1:B) {
  #     #--------------data--------------
  #     idx <- sample(1:length(Y), length(Y), replace = T)
  #     Xb <- X[idx,]
  #     Yb <- Y[idx]
  #     
  #     # generate basis matrix
  #     if (length(basis_list) > 0) {
  #       x_basis <- hal9001::make_design_matrix(as.matrix(Xb), basis_list)
  #       # unique_columns <- as.numeric(names(copy_map))
  #       # x_basis <- x_basis[, unique_columns]
  #     } else {
  #       x_basis <- matrix(1, ncol = 2, nrow = nrow(Xb))
  #     }
  #     x_basis <- as.matrix(x_basis)
  #     
  #     #--------------fit--------------
  #     is_glmnet = T
  #     if (dim(x_basis)[2] <= 1) {
  #       # dim of X_basis < 2. make it larger
  #       x_basis <- cbind(matrix(1, ncol = 1, nrow = nrow(X)), x_basis)
  #       x_basis <- cbind(matrix(0, ncol = 1, nrow = nrow(X)), x_basis)
  #       lasso_fit <- glm(Yb ~ x_basis, x = FALSE, y = FALSE, family = y_type)
  #       is_glmnet = F
  #     } else {
  #       lasso_fit <- tryCatch({
  #         lasso_fit <- glmnet::glmnet(x = x_basis, y = Yb, 
  #                                     family = y_type, 
  #                                     lambda = hal_fit$lambda_star,
  #                                     intercept = FALSE, standardize = FALSE)
  #       },
  #       error = function(){
  #         lasso_fit <- glm.fit(x = x_basis, y = Yb, family = y_type)
  #         lasso_fit <- glm(Yb ~ x_basis, x = FALSE, y = FALSE, family = y_type)
  #         is_glmnet = F
  #       })
  #     }
  #     
  #     #--------------predictions, 95% lower and upper bounds--------------
  #     y_hat_b <- c()
  #     for (i in 1:length(eval_points)) {
  #       X_new <- Xb
  #       X_new$A = eval_points[i]
  #       
  #       # generate basis matrix
  #       if (length(basis_list) > 0){
  #         x_basis_a <- hal9001::make_design_matrix(as.matrix(X_new), hal_fit$basis_list)
  #         # x_basis_a <- hal9001::apply_copy_map(x_basis_a, hal_fit$copy_map)
  #         
  #         if(dim(x_basis)[2] <= 1){
  #           x_basis_a <- cbind(matrix(1, ncol = 1, nrow = nrow(X)), x_basis_a)
  #           x_basis_a <- cbind(matrix(0, ncol = 1, nrow = nrow(X)), x_basis_a)
  #         }
  #       } else {
  #         x_basis_a <- matrix(1, ncol = 2, nrow = nrow(X_new))
  #       }
  # 
  #       # prediction
  #       #-----
  #       if (y_type == "binomial") {
  #         preds <- predict(lasso_fit, x_basis_a, type = "response")
  #       } else {
  #         preds <- predict(lasso_fit, x_basis_a)
  #       }
  #       #-----
  #       # beta_hat <- stats::coef(lasso_fit)
  #       # beta_hat[is.na(beta_hat)] <- 0
  #       # beta_hat <- as.matrix(beta_hat)
  #       # preds <- as.vector(
  #       #   Matrix::tcrossprod(x = x_basis_a, y = beta_hat[-1]) + beta_hat[1]
  #       # )
  #       # if (y_type == "binomial") preds <- stats::plogis(preds)
  #       #-----
  #       y_hat_b[i] = mean(preds)
  #     }
  #     y_hat_B <- rbind(y_hat_B, y_hat_b)
  #   }
  #   
  # y_hat_B <- y_hat_B[-1,]
  # lower_bd <- apply(y_hat_B, 2, quantile, probs = 0.05/2)
  # upper_bd <- apply(y_hat_B, 2, quantile, probs = 1-0.05/2)
  # 
  # return(list(lower_bd=lower_bd, upper_bd=upper_bd))
  
}

bootstrap_inference_u_l <- function(X, Y, eval_points, hal_fit, y_type, lambda_u_l_idx, B = 200){
  
  y_hat_B <- matrix(ncol = length(eval_points))
  
  for (b in 1:B) {
    #--------------data--------------
    idx <- sample(1:length(Y), length(Y), replace = T)
    Xb <- X[idx,]
    Yb <- Y[idx]
    
    
    x_basis <- list()
    lasso_fit <- list()
    
    for (i in 1:length(hal_fit)) {
      
      # generate basis matrix
      basis_list <- hal_fit[[i]]$basis_list
      
      # --------------fit glmnet--------------
      if (length(basis_list) > 0) {
        
        x_basis_l <- hal9001::make_design_matrix(as.matrix(Xb), basis_list)
        
        lasso_fit_l <- tryCatch({
          lasso_fit_l <- glmnet::glmnet(x = x_basis_l, 
                                        y = Yb, 
                                        family = y_type, 
                                        lambda = hal_fit[[i]]$lambda_star,
                                        intercept = FALSE, standardize = FALSE)
        },
        error = function(){
          lasso_fit_l <- NA
        })
        
      } else {
        
        x_basis_l <- NA
        lasso_fit_l <- NA
        
      } 
      
      x_basis[[i]] = x_basis_l
      lasso_fit[[i]] = lasso_fit_l
      
    }
    
    #--------------predictions, 95% lower and upper bounds--------------
    y_hat_b <- c()
    for (i in 1:length(eval_points)) {
      X_new <- Xb
      X_new$A = eval_points[i]
      u_l_idx = lambda_u_l_idx[i]
      
      if (length(basis_list) > 0){
        
        # generate basis matrix
        x_basis_a <- hal9001::make_design_matrix(as.matrix(X_new), hal_fit[[u_l_idx]]$basis_list)
        
        # prediction
        if(any(!is.na(lasso_fit[[u_l_idx]]))){
          
          if (y_type == "binomial") {
            preds <- predict(lasso_fit[[u_l_idx]], x_basis_a, type = "response")
          } else {
            preds <- predict(lasso_fit[[u_l_idx]], x_basis_a)
          }
          
          y_hat_b[i] = mean(preds)
          
        } else {
          y_hat_b[i] = NA
        }
        
      } else {
        y_hat_b[i] = NA
      }
      
    }
    y_hat_B <- rbind(y_hat_B, y_hat_b)
  }
  
  y_hat_B <- y_hat_B[-1,]
  lower_bd <- apply(y_hat_B, 2, quantile, probs = 0.05/2, na.rm = T)
  upper_bd <- apply(y_hat_B, 2, quantile, probs = 1-0.05/2, na.rm = T)
  
  return(list(lower_bd=lower_bd, upper_bd=upper_bd))
  #   
  #   y_hat_B <- matrix(ncol = length(eval_points))
  #   for (b in 1:B) {
  #     #--------------data--------------
  #     idx <- sample(1:length(Y), length(Y), replace = T)
  #     Xb <- X[idx,]
  #     Yb <- Y[idx]
  #     
  #     
  #     x_basis <- list()
  #     is_glmnet <- list()
  #     lasso_fit <- list()
  #     
  #     for (i in 1:length(hal_fit)) {
  #       
  #       # generate basis matrix
  #       basis_list <- hal_fit[[i]]$basis_list
  #       if (length(basis_list) > 0) {
  #         x_basis_l <- hal9001::make_design_matrix(as.matrix(Xb), basis_list)
  #       } else {
  #         x_basis_l <- matrix(1, ncol = 2, nrow = nrow(Xb))
  #       }
  #       x_basis[[i]] <- as.matrix(x_basis_l)
  #       
  #       #--------------fit--------------
  #       is_glmnet_l = T
  #       
  #       if (dim(x_basis_l)[2] <= 1) {
  #         # dim of X_basis < 2. make it larger
  #         x_basis_l <- cbind(matrix(1, ncol = 1, nrow = nrow(X)), x_basis_l)
  #         x_basis_l <- cbind(matrix(0, ncol = 1, nrow = nrow(X)), x_basis_l)
  #         lasso_fit_l <- glm(Yb ~ x_basis_l, x = FALSE, y = FALSE, family = y_type)
  #         is_glmnet_l = F
  #       } else {
  #         lasso_fit_l <- tryCatch({
  #           lasso_fit_l <- glmnet::glmnet(x = x_basis_l, y = Yb, 
  #                                       family = y_type, 
  #                                       lambda = hal_fit[[i]]$lambda_star,
  #                                       intercept = FALSE, standardize = FALSE)
  #         },
  #         error = function(){
  #           lasso_fit_l <- glm.fit(x = x_basis_l, y = Yb, family = y_type)
  #           lasso_fit_l <- glm(Yb ~ x_basis_l, x = FALSE, y = FALSE, family = y_type)
  #           is_glmnet_l = F
  #         })
  #       }
  #       is_glmnet[[i]] = is_glmnet_l
  #       lasso_fit[[i]] = lasso_fit_l
  #       
  #     }
  #     
  #     
  #     #--------------predictions, 95% lower and upper bounds--------------
  #     y_hat_b <- c()
  #     for (i in 1:length(eval_points)) {
  #       X_new <- Xb
  #       X_new$A = eval_points[i]
  #       u_l_idx = lambda_u_l_idx[i]
  #       
  #       # generate basis matrix
  #       if (length(basis_list) > 0){
  #         x_basis_a <- hal9001::make_design_matrix(as.matrix(X_new), hal_fit[[u_l_idx]]$basis_list)
  # 
  #         if(dim(x_basis[[u_l_idx]])[2] <= 1){
  #           x_basis_a <- cbind(matrix(1, ncol = 1, nrow = nrow(X)), x_basis_a)
  #           x_basis_a <- cbind(matrix(0, ncol = 1, nrow = nrow(X)), x_basis_a)
  #         }
  #       } else {
  #         x_basis_a <- matrix(1, ncol = 2, nrow = nrow(X_new))
  #       }
  #       
  #       # prediction
  #       if (y_type == "binomial") {
  #         preds <- predict(lasso_fit[[u_l_idx]], x_basis_a, type = "response")
  #       } else {
  #         preds <- predict(lasso_fit[[u_l_idx]], x_basis_a)
  #       }
  # 
  #       
  #       y_hat_b[i] = mean(preds)
  #     }
  #     y_hat_B <- rbind(y_hat_B, y_hat_b)
  #   }
  #   
  #   y_hat_B <- y_hat_B[-1,]
  #   lower_bd <- apply(y_hat_B, 2, quantile, probs = 0.05/2)
  #   upper_bd <- apply(y_hat_B, 2, quantile, probs = 1-0.05/2)
  #   
  #   return(list(lower_bd=lower_bd, upper_bd=upper_bd))
  #   
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
  x_names = names(obs)[names(obs) != 'Y']
  y_type = "binomial"

  Y <- as.numeric(as.matrix(obs %>% select(all_of(y_name))))
  X <- obs %>% 
    select(all_of(x_names)) %>% 
    mutate_if(sapply(., is.factor), as.numeric)
  
  
  fit_hal_all_criteria_rslts <- fit_hal_all_criteria(X, Y, y_type)
  
  #================================CV-HAL================================
  lambda_CV <- fit_hal_all_criteria_rslts$lambda_list$lambda_CV
  psi_hat <- sapply(eval_points, function(a){ X_new <- X
                                              X_new$A = a
                                              mean(predict(fit_hal_all_criteria_rslts$hal_fit_list$hal_CV, new_data = X_new)) } )

  psi_hat_pnt_cv <- cbind(eval_points, matrix(psi_hat, ncol=1), lambda_CV, 1, fit_hal_all_criteria_rslts$hal_fit_time_list$hal_cv_fit_time)
  colnames(psi_hat_pnt_cv) <- c("a", "Y_hat", "lambda", "lambda_scaler", "hal_fit_time")
  
  # IC-based inference
  psi_hat_pnt_cv_se <- IC_based_se(fit_hal_all_criteria_rslts$hal_fit_list$hal_CV, eval_points)
  psi_hat_pnt_cv <- as.data.frame(psi_hat_pnt_cv) %>% 
    mutate(SE = psi_hat_pnt_cv_se,
           ci_lwr = Y_hat - 1.96 * SE,
           ci_upr = Y_hat + 1.96 * SE)

  
  # bootstrap-based inference
  psi_hat_pnt_cv_bt_bds <- bootstrap_inference(X, Y, eval_points, fit_hal_all_criteria_rslts$hal_fit_list$hal_CV, y_type)
  psi_hat_pnt_cv$ci_lwr_bt <- psi_hat_pnt_cv_bt_bds$lower_bd
  psi_hat_pnt_cv$ci_upr_bt <- psi_hat_pnt_cv_bt_bds$upper_bd
  
  
  #================================global undersmoothing================================
  
  if(any(fit_hal_all_criteria_rslts$hal_fit_list$hal_u_g$coefs[-1] != 0)){
    
    psi_hat <- sapply(eval_points, function(a){ X_new <- X
    X_new$A = a
    mean(predict(fit_hal_all_criteria_rslts$hal_fit_list$hal_u_g, new_data = X_new)) } )
    
    lambda_scaler = fit_hal_all_criteria_rslts$lambda_list$lambda_u_g / lambda_CV
    psi_hat_pnt_u_g <- cbind(eval_points, matrix(psi_hat, ncol=1), fit_hal_all_criteria_rslts$lambda_list$lambda_u_g, lambda_scaler, fit_hal_all_criteria_rslts$hal_fit_time_list$hal_u_g_fit_time)
    colnames(psi_hat_pnt_u_g) <- c("a", "Y_hat", "lambda", "lambda_scaler", "hal_fit_time")
    
    # IC-based inference
    psi_hat_pnt_u_g_se <- IC_based_se(fit_hal_all_criteria_rslts$hal_fit_list$hal_CV, eval_points)
    psi_hat_pnt_u_g <- as.data.frame(psi_hat_pnt_u_g) %>% 
      mutate(SE = psi_hat_pnt_u_g_se,
             ci_lwr = Y_hat - 1.96 * SE,
             ci_upr = Y_hat + 1.96 * SE)
    
    # bootstrap-based inference
    psi_hat_pnt_u_g_bt_bds <- bootstrap_inference(X, Y, eval_points, fit_hal_all_criteria_rslts$hal_fit_list$hal_u_g, y_type)
    psi_hat_pnt_u_g$ci_lwr_bt <- psi_hat_pnt_u_g_bt_bds$lower_bd
    psi_hat_pnt_u_g$ci_upr_bt <- psi_hat_pnt_u_g_bt_bds$upper_bd
    
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
      SE <- IC_based_se_u_l(hal_fit, a)
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
    
    result_list[[b]] <- result
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
      psi_hat_pnt_scaled_se <- IC_based_se(hal_fit, eval_points)
      psi_hat_pnt_scaled <- as.data.frame(psi_hat_pnt_scaled) %>% 
        mutate(SE = psi_hat_pnt_scaled_se,
               ci_lwr = Y_hat - 1.96 * SE,
               ci_upr = Y_hat + 1.96 * SE)
      
      
      # bootstrap-based inference
      psi_hat_pnt_scaled_bt_bds <- bootstrap_inference(X, Y, eval_points, hal_fit, y_type)
      psi_hat_pnt_scaled$ci_lwr_bt <- psi_hat_pnt_scaled_bt_bds$lower_bd
      psi_hat_pnt_scaled$ci_upr_bt <- psi_hat_pnt_scaled_bt_bds$upper_bd
      
    } else {
      psi_hat_pnt_scaled = as.data.frame(cbind(eval_points, NA, lambda_CV, lambda_scaler, hal_scaled_fit_time, NA, NA, NA, NA, NA))
      names(psi_hat_pnt_scaled) = c("a", "Y_hat", "lambda", "lambda_scaler", "hal_fit_time", "SE", "ci_lwr", "ci_upr", "ci_lwr_bt", "ci_upr_bt")
    }
    results[[i]] = psi_hat_pnt_scaled
  }
  
  names(results) = paste0("scale=", round(lambda_scalers, 4))
  
  return(results)
}


run_simu_scaled_rep <- function(gen_data_functions, n, B, return_all_rslts = F){
  lambda_scalers = c(1.2, 1.1, 10^seq(from=0, to=-3, length=20))
  result_list <- list()
  for(b in 1:B){
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
    
    result_list[[b]] <- result
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

