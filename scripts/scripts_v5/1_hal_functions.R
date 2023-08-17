

###############################################################################
#'  implement undersmoothed HAL
#'
#' @details fit lasso with a sequence of candidates lambdas using \code{\link[glmnet]{glmnet}})
#' check a global criterion (\eqn{P_n(\phi_{s,i}(Y-\bar{Q}_{n,\lambda}\leq \freq{\sigma_n}{\sqrt{n}log(n)}))})
#' and select the largest lambda which satisfies the criterion.
#' @param X An input \code{matrix} with dimensions number of observations -by-
#'  number of covariates that will be used to derive the design matrix of basis
#'  functions.
#' @param Y A \code{numeric} vector of observations of the outcome variable.
#' @param fit_init The initial HAL fit object from the output list of \code{undersmooth_init}.
#' @param basis_mat The selected basis matrix from initial fit for undersmoothing,
#'  obtained from the output list of \code{undersmooth_init}.
#' @param Nlam Number of lambda candidates. The sequence ranges from \code{fit_init$lambda_star} to
#' \code{fit_init$lambda_star*10^(-3)}.
#' @param family A \code{character} or a \code{\link[stats]{family}} object
#'  (supported by \code{\link[glmnet]{glmnet}}) specifying the error/link
#'  family for a generalized linear model. See more details in \code{fit_hal}.

undersmooth_hal <- function(X,
                            Y,
                            fit_init,
                            basis_mat,
                            a = NULL,
                            Nlam = 20,
                            family = "gaussian"){
  
  if (family != "binomial"){
    preds_init <- predict(fit_init, new_data = X)
  }else {
    preds_init <- predict(fit_init, new_data = X, type = "response")
  }
  
  # estimates of sd in each direction using initial fit
  resid_init <- preds_init - Y
  sd_est  <- apply(basis_mat, 2, function(u) sd(resid_init*u))
  
  # refit on new lambda sequence
  us_lambda <- fit_init$lambda_star*10^seq(from=0, to=-3, length=Nlam)
  us_fit <- glmnet(fit_init$x_basis, Y, lambda=us_lambda, family = family, standardize = FALSE)

  
  if(identical(us_fit$df, 0)){
    res <- list("lambda_init" = fit_init$lambda_star,
                "lambda_under" = fit_init$lambda_star,
                "spec_under" = NA)
    return(res)
  }
  

  if (family != "binomial"){
    pred_mat <- predict(us_fit, fit_init$x_basis)
  }else {
    pred_mat <- predict(us_fit, fit_init$x_basis, type = "response")
  }
  resid_mat <- pred_mat - Y
  
  max_score <- get_maxscore(basis_mat = basis_mat,
                            resid_mat = resid_mat,
                            sd_est = sd_est,
                            Nlam = Nlam, us_fit = us_fit)
  
  # get the first lambda that satisfies the criteria
  lambda_under <- us_lambda[max_score <= 1/(sqrt(n)*log(n))][1] # over under-smoothing 
  
  
  # collect results
  coef_mat <- as.matrix(us_fit$beta)
  
  spec_under <- list("lambda" = us_lambda,
                     "l1_norm" = NA,
                     "n_coef" = NA)
  
  spec_under$l1_norm <- apply(coef_mat, 2, function(x){sum(abs(x))})
  spec_under$n_coef <- apply(coef_mat, 2, function(x){sum(x != 0)})
  
  res <- list("lambda_init" = fit_init$lambda_star,
              "lambda_under" = lambda_under,
              "spec_under" = spec_under)
  return(res)
}

###############################################################################
#'  undersoomthed HAL helper function for global criterion
#'
#' @details For each candidate lambda, do:
#'     1). standardize the score formed by each basis.
#'     2). calculate the mean of the standardized scores for each basis.
#' Select the max of the mean.
#' @param basis_mat The selected basis matrix from initial fit for undersmoothing,
#'  obtained from the output list of \code{undersmooth_init}.
#' @param resid_mat The residual matrix with each column the residuals correspongding to a lambda.
#' @param sd_est A numeric vector containing the sd of each column of \code{basis_mat}.
#' @param Nlam Number of lambda candidates.
#' @param us_fit The \code{glmnet} fit of the sequence of candidate lambdas.

get_maxscore <- function(basis_mat, resid_mat, sd_est, Nlam, us_fit){
  
  basis_mat_sd <- sweep(basis_mat, 2, sd_est, FUN = '/')
  score_all <- apply(basis_mat_sd, 2, function(u) {
    score_mat <- resid_mat * u
    score_mean <- apply(score_mat, 2, mean)
  })
  # absolute value
  max_score <- apply(abs(score_all), 1, max, na.rm=T)
  return(max_score)
}



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
                                       basis_mat = CV_basis_mat,
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
                                               basis_mat = CV_basis_mat_local,
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
                                       basis_mat = CV_basis_mat,
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

###############################################################################
#'  With given fitted HAL object and evaluation points, return the empirical SE

IC_based_se <- function(X, Y, hal_fit, eval_points){
  
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
    
    if(any(! is.na(IC_beta))){
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
    
  } else {
    se <- NA
  }
  
  return(se)
}

###############################################################################

IC_based_se_u_l <- function(X, Y, hal_fit, single_eval_point){
  
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
    X_new$A = single_eval_point
    
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



###############################################################################
# calculating efficient influence curves

cal_IC_for_beta <- function(X, Y, Y_hat, beta_n, family = 'binomial'){
  n <- dim(X)[1] 
  p <- length(beta_n)
  if (!is.matrix(X)) X <- as.matrix(X)
  
  # 1. calculate score: X'(Y - phi(X))
  res <- Y-Y_hat
  score <- sweep(t(X), 2, res, `*`)
  
  # 2. calculate the derivative of phi:
  if(family == 'binomial'){
    d_phi_scaler <- as.vector(exp(- beta_n %*% t(X)) / ((1 + exp(- beta_n %*% t(X)))^2))
    d_phi <- sweep(X, 1, d_phi_scaler, `*`)
  } else {
    d_phi = - X
  }
  
  # 3. -E_{P_n}(X d_phi)^(-1)
  tmat <- t(X) %*% d_phi / n
  if(! is.matrix(try(solve(tmat), silent = TRUE))){
    return(NA)
  }
  tmat <- -solve(tmat)
  
  # 4. calculate influence curves
  IC <- tmat %*% score
  
  return(IC)
}


cal_IC_for_EY <- function(X_new, beta_n, IC_beta, family = 'binomial'){
  if (!is.matrix(X_new)) X_new <- as.matrix(X_new)
  
  if(family == 'binomial'){
    d_phi_scaler_new <- as.vector(exp(- beta_n %*% t(X_new)) / ((1 + exp(- beta_n %*% t(X_new)))^2))
    d_phi_new <- sweep(X_new, 1, d_phi_scaler_new, `*`)
  } else {
    d_phi_new = X_new
  }
  
  IC = diag(d_phi_new %*% IC_beta)
  
  return(IC)
}


cal_IC_for_ATE <- function(X_new_a, X_new_0, beta_n, IC_beta, family = 'binomial'){
  
  if (!is.matrix(X_new_a)) X_new_a <- as.matrix(X_new_a)
  if (!is.matrix(X_new_0)) X_new_0 <- as.matrix(X_new_0)
  
  if (family == 'binomial') {
    d_phi_scaler_new_a <- as.vector(exp(- beta_n %*% t(X_new_a)) / ((1 + exp(- beta_n %*% t(X_new_a)))^2))
    d_phi_new_a <- sweep(X_new_a, 1, d_phi_scaler_new_a, `*`)
    
    d_phi_scaler_new_0 <- as.vector(exp(- beta_n %*% t(X_new_0)) / ((1 + exp(- beta_n %*% t(X_new_0)))^2))
    d_phi_new_0 <- sweep(X_new_0, 1, d_phi_scaler_new_0, `*`)
    
    d_phi_new <- d_phi_new_a - d_phi_new_0
  } else {
    d_phi_new <- X_new_a - X_new_0
  }
  
  IC = diag(d_phi_new %*% IC_beta)
  
  return(IC)
}


###############################################################################

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
      
      # maybe
      
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
  SE <- sqrt(apply(y_hat_B, 2, var, na.rm = T))
  
  return(list(lower_bd=lower_bd, upper_bd=upper_bd, SE=SE))
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


###############################################################################

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
  SE <- sqrt(apply(y_hat_B, 2, var, na.rm = T))
  
  return(list(lower_bd=lower_bd, upper_bd=upper_bd, SE=SE))
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

