options(scipen = 999)

########################
# fit undersmoothed HAL with initial lambda * 10^k, 
########################
undersmooth_hal_fit <- function(n, X, Y, y_type, lambda_scaler){
  
  hal_init <- undersmooth_init(X, Y, family = y_type)
  init_lambda <- hal_init$fit_init$lambda_star
  sclaed_lambda <- init_lambda*lambda_scaler

  print(paste0("=> Fitting initial HAL with lambda:       ", round(init_lambda,5)))
  dgd_hal_fit <- fit_hal(X = X,
                         Y = Y,
                         family = y_type,
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
                         lambda = init_lambda
  )
  y_preds <- predict(dgd_hal_fit, new_data = X)
  mse <- sum((y_preds - Y)^2)
  auc <- auc(Y, y_preds)
  print(paste0("    MSE: ", round(mse, 4), ", AUC: ", round(auc, 4)))
  
  print(paste0("=> Fitting HAL with scaled lambda: ", round(sclaed_lambda,5)))
  dgd_under_hal_fit <- fit_hal(X = X,
                               Y = Y,
                               family = y_type,
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
                               lambda = sclaed_lambda
  )
  y_preds <- predict(dgd_under_hal_fit, new_data = X)
  mse <- sum((y_preds - Y)^2)
  auc <- auc(Y, y_preds)
  print(paste0("    MSE: ", round(mse, 4), ", AUC: ", round(auc, 4)))
  
  hal_fit <- list(dgd_hal_fit, dgd_under_hal_fit)
  names(hal_fit) <- c("dgd_hal_fit", "dgd_under_hal_fit")
  
  lambdas <- list(sclaed_lambda,
                  init_lambda)
  names(lambdas) <- c("lambda_under","lambda_init")
  
  output <- list(hal_fit, lambdas)
  names(output) <- c("hal_fit", "lambdas")
  
  return(output)
}


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
  
  # 3. E_{P_n}(X d_phi)^(-1)
  tmat <- solve(t(X) %*% d_phi / n)
  
  # 4. calculate influence curves
  IC <- t(tmat %*% t(score))
  
  return(IC)
}

cal_IC_for_phi <- function(X_new, beta_n, IC_beta){
  
  if (!is.matrix(X_new)) X_new <- as.matrix(X_new)
  
  d_phi_scaler_new <- as.vector(exp(- beta_n %*% t(X_new)) / ((1 + exp(- beta_n %*% t(X_new)))^2))
  d_phi_new <- sweep(X_new, 1, d_phi_scaler_new, `*`)
  
  IC = diag(d_phi_new %*% t(IC_beta))
  
  return(IC)
}


########################
# calculating ate
########################
cal_ate_hal <- function(dgd_hal_fit, dgd_under_hal_fit,
                        n, X, Y,
                        a1, a0, z){
  
  X1 <- X
  X1$A = a1
  X1$Z = z
  
  X0 <- X
  X0$A = a0
  X0$Z = z
  
  y1_preds <- predict(dgd_hal_fit, new_data = X1)
  y0_preds <- predict(dgd_hal_fit, new_data = X0)
  psi_ss <- mean(y1_preds - y0_preds)
  
  y1_preds_under <- predict(dgd_under_hal_fit, new_data = X1)
  y0_preds_under <- predict(dgd_under_hal_fit, new_data = X0)
  psi_ss_under <- mean(y1_preds_under - y0_preds_under)

  return(list("psi_ss_under" = psi_ss_under,
              "psi_ss_init" = psi_ss))
}

# calculating the ate with given a values 
# z =0 or 1
ate_hal_avec0_z01 <- function(a_vec, obs_df, lambda_scaler){
  
  ate_a_0_ss_under_hal <- rep(NA, length(a_vec))
  ate_a_1_ss_under_hal <- rep(NA, length(a_vec))
  ate_a_0_ss_hal <- rep(NA, length(a_vec))
  ate_a_1_ss_hal <- rep(NA, length(a_vec))
  
  y_name = "Y"
  x_names = c("W", "A", "Z")
  n <- nrow(obs_df)
  Y <- as.numeric(as.matrix(obs_df %>% select(all_of(y_name))))
  X <- obs_df %>% 
    select(all_of(x_names)) %>% 
    mutate_if(sapply(., is.factor), as.numeric)
  
  undersmooth_hal_fit_rlts <-  undersmooth_hal_fit(n, X, Y, y_type = "binomial", lambda_scaler=lambda_scaler)
  hal_fit <- undersmooth_hal_fit_rlts$hal_fit
  dgd_hal_fit <- hal_fit$dgd_hal_fit
  dgd_under_hal_fit <- hal_fit$dgd_under_hal_fit

  for (i in 1:length(a_vec)) {
    a <- a_vec[i]
    
    ate_hal_a_0_relts <- cal_ate_hal(
      dgd_hal_fit = dgd_hal_fit,
      dgd_under_hal_fit = dgd_under_hal_fit,
      n, X, Y,
      a1 = a, a0 = 0, z = 0
    )
    
    ate_hal_a_1_relts <- cal_ate_hal(
      dgd_hal_fit = dgd_hal_fit,
      dgd_under_hal_fit = dgd_under_hal_fit,
      n, X, Y,
      a1 = a, a0 = 0, z = 1
    )  
    
    ate_a_1_ss_under_hal[i] = ate_hal_a_1_relts$psi_ss_under
    ate_a_0_ss_under_hal[i] = ate_hal_a_0_relts$psi_ss_under
    
    ate_a_1_ss_hal[i] = ate_hal_a_1_relts$psi_ss_init
    ate_a_0_ss_hal[i] = ate_hal_a_0_relts$psi_ss_init  
  }
  
  ate_hal_results <- list(ate_a_1_ss_under_hal,
                          ate_a_1_ss_hal,
                          ate_a_0_ss_under_hal,
                          ate_a_0_ss_hal,
                          undersmooth_hal_fit_rlts$lambdas$lambda_under,
                          undersmooth_hal_fit_rlts$lambdas$lambda_init)
  names(ate_hal_results) <- c("ate_a_1_ss_under_hal",
                              "ate_a_1_ss_hal",
                              "ate_a_0_ss_under_hal",
                              "ate_a_0_ss_hal",
                              "lambda_under",
                              "lambda_init")
  return(ate_hal_results)
}


# calculating the ate with given a values using glm
# z =0 or 1
ate_glm_avec0_z01 <- function(a_vec, obs_df){
  
  # fit glm with the true formula
  print("  Fitting glm:       ")
  glm_fit <- glm(formula = Y ~ W + A + W*A + Z, 
                 family = binomial,
                 data = obs_df)
  
  y_preds <- predict(glm_fit)
  mse <- sum((y_preds - obs_df$Y)^2)
  auc <- auc(obs_df$Y, y_preds)
  print(paste0("    MSE: ", round(mse, 4), ", AUC: ", round(auc, 4)))
  
  # calculate the ATE with given a_vec
  ate_glm_a_0 <- rep(NA, length(a_vec))
  ate_glm_a_1 <- rep(NA, length(a_vec))
  
  for (i in 1:length(a_vec)) {
    a <- a_vec[i]
    
    X1 <- obs_df
    X1$A = a
    X0 <- obs_df
    X0$A = 0
    
    # for Z ==1
    X1$Z = 1
    X0$Z = 1
    
    y1_preds <- predict(glm_fit, new_data = X1)
    y0_preds <- predict(glm_fit, new_data = X0)
    psi_ss_z1 <- mean(y1_preds - y0_preds)

    # for Z = 0
    X1$Z = 0
    X0$Z = 0
    
    y1_preds <- predict(glm_fit, new_data = X1)
    y0_preds <- predict(glm_fit, new_data = X0)
    psi_ss_z0 <- mean(y1_preds - y0_preds)

    
    ate_glm_a_1[i] = psi_ss_z1
    ate_glm_a_0[i] = psi_ss_z0
  }
  
  ate_glm_results <- list(ate_glm_a_1,
                          ate_glm_a_0)
  names(ate_glm_results) <- c("ate_glm_a_1",
                              "ate_glm_a_0")
  return(ate_glm_results)
}

ci_cov_rate <- function(vals, sd, true_val){
  cov_rate <- rep(NA, length(a_vec))
  for (i in 1:length(a_vec)) {
    ci_l <- vals[,i] - 1.96 * sd[i]
    ci_u <- vals[,i] + 1.96 * sd[i]
    cov_rate[i] <- mean(ci_l <= true_val[i] & true_val[i] <= ci_u)
  }
  return(cov_rate)
}

##############################################################
# This function runs the simulation for B rounds.
# For each round, it will generate a new sample with size n from the provided 
# generate_data function.

# Output: 
# A list with two data frames containing the average of estimated values of 
# the B rounds, using HAL (initial and undersmoothed). 
# One with z=1, and one with z=0.
##############################################################
run_simu <- function(generate_data, n, B, get_estimates = FALSE){
  
  results_list <- list()
  
  for(b in 1:B){
    print(paste0("---------------------starting round ", b, "---------------------"))
    Obs1_n <- generate_data(n)
    
    # CV HAL
    CV-hal <- undersmooth_init(X, Y, family = y_type)
    CV_lambda <- hal_init$fit_init$lambda_star
    
    print(paste0("=> Fitting initial HAL with lambda:       ", round(init_lambda,5)))
    dgd_hal_fit <- fit_hal(X = X,
                           Y = Y,
                           family = y_type,
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
                           lambda = init_lambda
    )
    y_preds <- predict(dgd_hal_fit, new_data = X)
    mse <- sum((y_preds - Y)^2)
    auc <- auc(Y, y_preds)
    print(paste0("    MSE: ", round(mse, 4), ", AUC: ", round(auc, 4)))
    
    
    
  }
  for(i in 1:length(lambda_scalers)){
    scale_i <- lambda_scalers[i]
    print(paste0("=============================initial lambda * ", scale_i, "============================="))
    
    ate_a_1_ss_under_hal_simu <- list()
    ate_a_1_ss_hal_simu <- list()
    ate_a_0_ss_under_hal_simu <- list()
    ate_a_0_ss_hal_simu <- list()
    
    lambda_under <- rep(NA, B)
    lambda_init <- rep(NA, B)
    
    for (b in 1:B) {
      
      
      
      ate_hal_a_z <- ate_hal_avec0_z01(a_vec, obs_df=Obs1_n, lambda_scaler = scale_i)
      # ate_glm_a_z <- ate_glm_avec0_z01(a_vec, obs_df=Obs1_n)
      
      ate_a_1_ss_under_hal_simu[[b]] = ate_hal_a_z$ate_a_1_ss_under_hal
      ate_a_0_ss_under_hal_simu[[b]] = ate_hal_a_z$ate_a_0_ss_under_hal
      ate_a_1_ss_hal_simu[[b]] = ate_hal_a_z$ate_a_1_ss_hal
      ate_a_0_ss_hal_simu[[b]] = ate_hal_a_z$ate_a_0_ss_hal
      # ate_glm_a_1[[b]] = ate_glm_a_z$ate_glm_a_1
      # ate_glm_a_0[[b]] = ate_glm_a_z$ate_glm_a_0
      
      lambda_under[b] = ate_hal_a_z$lambda_under
      lambda_init[b] = ate_hal_a_z$lambda_init
    }
    
    ate_a_1_ss_under_hal_simu <- do.call("rbind", ate_a_1_ss_under_hal_simu)
    ate_a_1_ss_under_hal_simu_means <- colMeans(ate_a_1_ss_under_hal_simu)
    ate_a_1_ss_under_hal_simu_sd <- apply(ate_a_1_ss_under_hal_simu, 2, sd)
    ate_a_1_ss_under_hal_simu_coverage_rate <- ci_cov_rate(vals=ate_a_1_ss_under_hal_simu, 
                                                           sd = ate_a_1_ss_under_hal_simu_sd, 
                                                           true_val = psi0_a_1)
    
    ate_a_0_ss_under_hal_simu <- do.call("rbind", ate_a_0_ss_under_hal_simu)
    ate_a_0_ss_under_hal_simu_means <- colMeans(ate_a_0_ss_under_hal_simu)
    ate_a_0_ss_under_hal_simu_sd <- apply(ate_a_0_ss_under_hal_simu, 2, sd)
    ate_a_0_ss_under_hal_simu_coverage_rate <- ci_cov_rate(vals=ate_a_0_ss_under_hal_simu, 
                                                           sd = ate_a_0_ss_under_hal_simu_sd, 
                                                           true_val = psi0_a_0)
    
    ate_a_1_ss_hal_simu <- do.call("rbind", ate_a_1_ss_hal_simu)
    ate_a_1_ss_hal_simu_means <- colMeans(ate_a_1_ss_hal_simu)
    ate_a_1_ss_hal_simu_sd <- apply(ate_a_1_ss_hal_simu, 2, sd)
    ate_a_1_ss_hal_simu_coverage_rate <- ci_cov_rate(vals=ate_a_1_ss_hal_simu, 
                                                     sd = ate_a_1_ss_hal_simu_sd, 
                                                     true_val = psi0_a_1)
    
    ate_a_0_ss_hal_simu <- do.call("rbind", ate_a_0_ss_hal_simu)
    ate_a_0_ss_hal_simu_means <- colMeans(ate_a_0_ss_hal_simu)
    ate_a_0_ss_hal_simu_sd <- apply(ate_a_0_ss_hal_simu, 2, sd)
    ate_a_0_ss_hal_simu_coverage_rate <- ci_cov_rate(vals=ate_a_0_ss_hal_simu, 
                                                     sd = ate_a_0_ss_hal_simu_sd, 
                                                     true_val = psi0_a_0)
    
    results_z1 <- data.frame(targ_par = paste0("ATE(",a_vec,",0)"), 
                             psi0 = psi0_a_1,
                             # ss_glm_mean = ate_glm_a_1_relts_means,
                             # ss_glm_sd = ate_glm_a_1_relts_sd,
                             mean_under = ate_a_1_ss_under_hal_simu_means,
                             sd_u = ate_a_1_ss_under_hal_simu_sd,
                             cov_rate_u = ate_a_1_ss_under_hal_simu_coverage_rate,
                             mean_init = ate_a_1_ss_hal_simu_means,
                             sd_i = ate_a_1_ss_hal_simu_sd,
                             cov_rate_i = ate_a_1_ss_hal_simu_coverage_rate)
    results_z1 <- results_z1 %>%
      mutate(bias_u = mean_under - psi0, 
             bias_i = mean_init - psi0,
             mse_u = bias_u^2 + sd_u^2,
             mse_i = bias_i^2 + sd_i^2) %>% 
      select(targ_par, psi0, mean_under, mean_init, bias_u, bias_i, sd_u, sd_i, mse_u, mse_i, cov_rate_u, cov_rate_i) 
    
    results_z0 <- data.frame(targ_par = paste0("ATE(",a_vec,",0)"), 
                             psi0 = psi0_a_0,
                             # ss_glm_mean = ate_glm_a_0_relts_means,
                             # ss_glm_sd = ate_glm_a_0_relts_sd,
                             mean_under = ate_a_0_ss_under_hal_simu_means,
                             sd_u = ate_a_0_ss_under_hal_simu_sd,
                             cov_rate_u = ate_a_1_ss_hal_simu_coverage_rate,
                             mean_init = ate_a_0_ss_hal_simu_means,
                             sd_i = ate_a_0_ss_hal_simu_sd,
                             cov_rate_i = ate_a_0_ss_hal_simu_coverage_rate)
    results_z0 <- results_z0 %>%
      mutate(bias_u = mean_under - psi0, 
             bias_i = mean_init - psi0,
             mse_u = bias_u^2 + sd_u^2,
             mse_i = bias_i^2 + sd_i^2) %>% 
      select(targ_par, psi0, mean_under, mean_init, bias_u, bias_i, sd_u, sd_i, mse_u, mse_i, cov_rate_u, cov_rate_i) 
    
    lambda_smry <- rbind(c(mean(lambda_under), sd(lambda_under)),
                         c(mean(lambda_init), sd(lambda_init)))
    colnames(lambda_smry) <- c("mean", "sd")
    row.names(lambda_smry) <- c("undersmoothed", "initial")
    
    ate_ests <- list(ate_a_1_ss_under_hal_simu,
                     ate_a_0_ss_under_hal_simu,
                     ate_a_1_ss_hal_simu,
                     ate_a_0_ss_hal_simu)
    names(ate_ests) <- c("ate_a_1_ss_under_hal_simu",
                         "ate_a_0_ss_under_hal_simu",
                         "ate_a_1_ss_hal_simu",
                         "ate_a_0_ss_hal_simu")
    results <- list("ATE_z=1"=results_z1, "ATE_z=0"=results_z0, "lambda_summary" = lambda_smry, "ate_ests" = ate_ests)
    
    results_list[[i]] <- results
  }
  
  names(results_list) <- paste0("lambda_scaler_", round(lambda_scalers, 4))
  return(results_list)
}







plot_perforences <- function(results_list, z_para, target_para){
  est_avg <- c()
  bias <- c()
  sd <- c()
  cr <- c()
  target <- c()
  z <- c()

  if(target_para == "ATE(0.5, 0)"){
    j = 1
  } else if (target_para == "ATE(5, 0)"){
    j = 10
  } else if (target_para == "ATE(2.5, 0)") {
    j = 5
  }
  psi0 <- ifelse(z_para==1, psi0_a_1[j], psi0_a_0[j])
  for(i in 1:length(lambda_scalers)){
    if(z_para == 1){
      est_avg <- c(est_avg, results_list[[i]]$`ATE_z=1`$mean_under[j])
      bias <- c(bias, results_list[[i]]$`ATE_z=1`$bias_u[j])
      sd <- c(sd, results_list[[i]]$`ATE_z=1`$sd_u[j])
      cr <- c(cr, results_list[[i]]$`ATE_z=1`$cov_rate_u[j])
      target <- c(target, target_para)
      z <- c(z, 1)
    } else {
      est_avg <- c(est_avg, results_list[[i]]$`ATE_z=0`$mean_under[j])
      bias <- c(bias, results_list[[i]]$`ATE_z=0`$bias_u[j])
      sd <- c(sd, results_list[[i]]$`ATE_z=0`$sd_u[j])
      cr <- c(cr, results_list[[i]]$`ATE_z=0`$cov_rate_u[j])
      target <- c(target, target_para)
      z <- c(z, 0)
    }
  }
  
  perform_df <- data.frame(lambda_scalers = lambda_scalers,
                           target = target,
                           z = z,
                           est_avg = est_avg,
                           bias = bias,
                           sd = sd,
                           cr = cr) %>%
    mutate(bias_d_df = abs(bias)/sd)
  
  
  p_est_avg <- ggplot(perform_df) +  
    geom_point(aes(x = lambda_scalers, y = est_avg)) + 
    geom_hline(aes(yintercept=psi0)) +
    labs(title="Estimation average") +
    theme()
  
  p_bias <- ggplot(perform_df, 
                   aes(x = lambda_scalers, y = bias)) +  
    geom_point() + 
    labs(title="Bias") 
  
  p_sd <- ggplot(perform_df, 
                 aes(x = lambda_scalers, y = sd)) +  
    geom_point() + 
    labs(title="Standard deviation") 
  
  p_bias_d_df <- ggplot(perform_df, 
                        aes(x = lambda_scalers, y = bias_d_df)) +  
    geom_point() + 
    labs(title="|Bias| / Standard deviation") 
  
  p_cr <- ggplot(perform_df, 
                 aes(x = lambda_scalers, y = cr)) +  
    geom_point() + 
    labs(title="Coverage rate") 
  
  p <- grid.arrange(p_est_avg, p_bias, p_sd, p_bias_d_df, p_cr,
                    layout_matrix = rbind(c(NA,1,1,NA),
                                          c(NA,1,1,NA),
                                          c(2,2,3,3),
                                          c(2,2,3,3),
                                          c(4,4,5,5),
                                          c(4,4,5,5)),
                    top = textGrob(paste0("HAL-based plug in estimator performence for ", target_para, ", z=", z_para), 
                                   gp=gpar(fontsize=11, fontface = 'bold')))
  return(p)
}
