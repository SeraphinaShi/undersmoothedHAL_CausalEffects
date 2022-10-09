
########################
# fit undersmoothed HAL
########################
undersmooth_hal_fit <- function(n, X, Y, y_type){
  
  hal_init <- undersmooth_init(X, Y, family = y_type)
  hal_undersmooth <- undersmooth_hal(X, Y,
                                     fit_init = hal_init$fit_init,
                                     basis_mat = hal_init$basis_mat,
                                     family = y_type)

  printf("  Fitting initial HAL with lambda:       %f \n", hal_undersmooth$lambda_init)
  dgd_hal_fit <- fit_hal(X = X,
                         Y = Y,
                         family = y_type,
                         num_knots = hal9001:::num_knots_generator(
                           max_degree = ifelse(ncol(X) >= 20, 2, 3),
                           smoothness_orders = 1,
                           base_num_knots_0 = 500,
                           base_num_knots_1 = 100
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
                         lambda = hal_undersmooth$lambda_init
  )
  
  printf("  Fitting undersmoothed HAL with lambda: %f \n", hal_undersmooth$lambda_under)
  dgd_under_hal_fit <- fit_hal(X = X,
                               Y = Y,
                               family = y_type,
                               num_knots = hal9001:::num_knots_generator(
                                 max_degree = ifelse(ncol(X) >= 20, 2, 3),
                                 smoothness_orders = 1,
                                 base_num_knots_0 = 500,
                                 base_num_knots_1 = 100
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
                               lambda = hal_undersmooth$lambda_under
  )
  
  hal_fit <- list(dgd_hal_fit, dgd_under_hal_fit)
  names(hal_fit) <- c("dgd_hal_fit", "dgd_under_hal_fit")
  
  return(hal_fit)
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
  
  y_preds <- predict(dgd_hal_fit, new_data = X)
  y1_preds <- predict(dgd_hal_fit, new_data = X1)
  y0_preds <- predict(dgd_hal_fit, new_data = X0)
  psi_ss <- mean(y1_preds - y0_preds)
  mse <- sum((y_preds - Y)^2)
  
  y_preds_under <- predict(dgd_under_hal_fit, new_data = X)
  y1_preds_under <- predict(dgd_under_hal_fit, new_data = X1)
  y0_preds_under <- predict(dgd_under_hal_fit, new_data = X0)
  psi_ss_under <- mean(y1_preds_under - y0_preds_under)
  mse_under <- sum((y_preds_under - Y)^2)
  
  return(list("psi_ss_under" = psi_ss_under,
              "mse_under" = mse_under,
              "psi_ss_init" = psi_ss,
              "mse_init" = mse))
}

# calculating the ate with given a vlues 
# z =0 or 1
ate_hal_avec0_z01 <- function(a_vec, obs_df){
  
  ate_hal_a_0_relts_ss_under_hal <- rep(NA, length(a_vec))
  ate_hal_a_1_relts_ss_under_hal <- rep(NA, length(a_vec))
  ate_hal_a_0_relts_ss_hal <- rep(NA, length(a_vec))
  ate_hal_a_1_relts_ss_hal <- rep(NA, length(a_vec))
  
  y_name = "Y"
  x_names = c("W", "A", "Z")
  n <- nrow(obs_df)
  Y <- as.numeric(as.matrix(obs_df %>% select(all_of(y_name))))
  X <- obs_df %>% 
    select(all_of(x_names)) %>% 
    mutate_if(sapply(., is.factor), as.numeric)
  
  hal_fit <- undersmooth_hal_fit(n, X, Y, y_type = "binomial")
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
    
    ate_hal_a_1_relts_ss_under_hal[i] = ate_hal_a_1_relts$psi_ss_under
    ate_hal_a_0_relts_ss_under_hal[i] = ate_hal_a_0_relts$psi_ss_under
    
    ate_hal_a_1_relts_ss_hal[i] = ate_hal_a_1_relts$psi_ss_init
    ate_hal_a_0_relts_ss_hal[i] = ate_hal_a_0_relts$psi_ss_init  
  }
  
  ate_hal_results <- list(ate_hal_a_1_relts_ss_under_hal,
                          ate_hal_a_1_relts_ss_hal,
                          ate_hal_a_0_relts_ss_under_hal,
                          ate_hal_a_0_relts_ss_hal)
  names(ate_hal_results) <- c("ate_hal_a_1_relts_ss_under_hal",
                              "ate_hal_a_1_relts_ss_hal",
                              "ate_hal_a_0_relts_ss_under_hal",
                              "ate_hal_a_0_relts_ss_hal")
  return(ate_hal_results)
}





run_simu <- function(generate_data, n, B){
  
  ate_hal_a_1_relts_ss_under_hal_simu <- list()
  ate_hal_a_1_relts_ss_hal_simu <- list()
  ate_hal_a_0_relts_ss_under_hal_simu <- list()
  ate_hal_a_0_relts_ss_hal_simu <- list()
  
  for (b in 1:B) {
    print(paste0("-----starting round ", b, "-----"))
    
    Obs1_n <- generate_data(n)
    ate_hal_a_z_n500 <- ate_hal_avec0_z01(a_vec, obs_df=Obs1_n)
    
    ate_hal_a_1_relts_ss_under_hal_simu[[b]] = ate_hal_a_z_n500$ate_hal_a_1_relts_ss_under_hal
    ate_hal_a_0_relts_ss_under_hal_simu[[b]] = ate_hal_a_z_n500$ate_hal_a_0_relts_ss_under_hal
    ate_hal_a_1_relts_ss_hal_simu[[b]] = ate_hal_a_z_n500$ate_hal_a_1_relts_ss_hal
    ate_hal_a_0_relts_ss_hal_simu[[b]] = ate_hal_a_z_n500$ate_hal_a_0_relts_ss_hal
  }
  
  ate_hal_a_1_relts_ss_under_hal_simu <- do.call("rbind", ate_hal_a_1_relts_ss_under_hal_simu)
  ate_hal_a_1_relts_ss_under_hal_simu_means <- colMeans(ate_hal_a_1_relts_ss_under_hal_simu)
  ate_hal_a_1_relts_ss_under_hal_simu_sd <- apply(ate_hal_a_1_relts_ss_under_hal_simu, 2, sd)
  
  ate_hal_a_0_relts_ss_under_hal_simu <- do.call("rbind", ate_hal_a_0_relts_ss_under_hal_simu)
  ate_hal_a_0_relts_ss_under_hal_simu_means <- colMeans(ate_hal_a_0_relts_ss_under_hal_simu)
  ate_hal_a_0_relts_ss_under_hal_simu_sd <- apply(ate_hal_a_0_relts_ss_under_hal_simu, 2, sd)
  
  ate_hal_a_1_relts_ss_hal_simu <- do.call("rbind", ate_hal_a_1_relts_ss_hal_simu)
  ate_hal_a_1_relts_ss_hal_simu_means <- colMeans(ate_hal_a_1_relts_ss_hal_simu)
  ate_hal_a_1_relts_ss_hal_simu_sd <- apply(ate_hal_a_1_relts_ss_hal_simu, 2, sd)
  
  ate_hal_a_0_relts_ss_hal_simu <- do.call("rbind", ate_hal_a_0_relts_ss_hal_simu)
  ate_hal_a_0_relts_ss_hal_simu_means <- colMeans(ate_hal_a_0_relts_ss_hal_simu)
  ate_hal_a_0_relts_ss_hal_simu_sd <- apply(ate_hal_a_0_relts_ss_hal_simu, 2, sd)
  
  results_z1 <- data.frame(targ_par = paste0("ATE(",a_vec,",0)"), 
                           psi0 = psi0_a_1,
                           ss_under_hal_mean = ate_hal_a_1_relts_ss_under_hal_simu_means,
                           ss_under_hal_sd = ate_hal_a_1_relts_ss_under_hal_simu_sd,
                           ss_hal_mean = ate_hal_a_1_relts_ss_hal_simu_means,
                           ss_hal_sd = ate_hal_a_1_relts_ss_hal_simu_sd)
  results_z0 <- data.frame(targ_par = paste0("ATE(",a_vec,",0)"), 
                           psi0 = psi0_a_0,
                           ss_under_hal_mean = ate_hal_a_0_relts_ss_under_hal_simu_means,
                           ss_under_hal_sd = ate_hal_a_0_relts_ss_under_hal_simu_sd,
                           ss_hal_mean = ate_hal_a_0_relts_ss_hal_simu_means,
                           ss_hal_sd = ate_hal_a_0_relts_ss_hal_simu_sd)
  results <- list("ATE_z=1"=results_z1, "ATE_z=0"=results_z0)
  return(results)
}
