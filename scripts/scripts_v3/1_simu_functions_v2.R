options(scipen = 999)

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
  
  y_name = "Y"
  x_names = c("W", "A", "Z")
  y_type = "binomial"
  
  results_list <- list()
  
  for(b in 1:B){
    print(paste0("---------------------starting round ", b, "---------------------"))
    obs_df <- generate_data(n)

    Y <- as.numeric(as.matrix(obs_df %>% select(all_of(y_name))))
    X <- obs_df %>% 
      select(all_of(x_names)) %>% 
      mutate_if(sapply(., is.factor), as.numeric)
    
    # ========================Lambdas========================
    # ------------------------CV HAL------------------------
    hal_CV <- undersmooth_init(X, Y, family = y_type)
    lambda_CV <- hal_CV$fit_init$lambda_star
    
    # ------------------------Globally undersmoothed HAL------------------------
    hal_undersmooth <- undersmooth_hal(X, Y,
                                       fit_init = hal_CV$fit_init,
                                       basis_mat = hal_CV$basis_mat,
                                       family = y_type)
    lambda_global_u <- hal_undersmooth$lambda_under
    
    # ------------------------Locally undersmoothed HAL------------------------
    # lambda_local_u <- rep(NA, 2*length(a_vec))
    # for (z in c(0,1)) {
    #   for (a in a_vec) {
    #     X_new = X
    #     X_new$A = a
    #     X_new$Z = z 
    #     
    #   }
    # }
    
    # ------------------------Over a grid of penalties------------------------
    lambda_grid <- lambda_CV*lambda_scalers
    
    # ===============Fit HAL working model with given penalties=================
    # ===============Get estimates from each HAL=================
    lambdas_names <- c("CV", "global_u", paste0("CV*", round(lambda_scalers, 4)))
    lambdas <- c(lambda_CV, lambda_global_u, lambda_grid)
    
    estimates <- list()
    for (i in 1:length(lambdas)) {
      HAL_fit <- fit_hal(X = X,
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
                         lambda = lambda[i]
      )
      ate_ests[[i]] <- cal_ate_hal(HAL_fit)
    }
    
    
  }
  
  return(results_list)
}





cal_ate_hal <- function(dgd_hal_fit){
  
  ate = c()
  for (z in c(0,1)) {
    for(a in a_vec){
      X1 <- X
      X1$A = a
      X1$Z = z
      
      X0 <- X
      X0$A = 0
      X0$Z = z
      
      y1_preds <- predict(dgd_hal_fit, new_data = X1)
      y0_preds <- predict(dgd_hal_fit, new_data = X0)
      psi <- mean(y1_preds - y0_preds)
      
      ate <- append(ate, psi)
    }
  }
  
  results <- list("ATE_ests_z=0"=ate[1:length(a_vec)], "ATE_ests_z=1"=ate[(length(a_vec)+1):length(ate)])
  return(results)
}


