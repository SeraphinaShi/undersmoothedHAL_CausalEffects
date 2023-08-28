
library(mgcv)


##############################################################
# This function runs the simulation with given:
#       - data generating function
#       - sample size n
# returns the estimated ATE and empirical 95% CI
##############################################################
run_simu_1round_gam_poly <- function(gen_data_functions, eval_points, y_type, n, method){
  
  obs <- gen_data_functions(n)
  
  y_name <- "Y"
  x_names <- names(obs)[names(obs) != y_name]
  
  if (method == 'GAM'){
    
    form <- as.formula(paste(y_name, "~", paste(paste0("s(", x_names, ")"), collapse = " + ")))
    model_fit <- gam(form, data = obs, family = binomial())
    
  } else if (method == 'POLY'){
    
    # Fit a polynomial logistic regression model
    form <- as.formula(paste(y_name, "~", paste(paste0("poly(", x_names, ", degree = 2, raw = TRUE)"), collapse = " + ")))
    model_fit <- glm(form, data = obs, family = binomial)
    
  }

  
  psi_hat_pnt <- matrix(NA, nrow = length(eval_points), ncol = 5)
  colnames(psi_hat_pnt) = c("a", "y_hat", "SE", "ci_lwr", "ci_upr")
  for (i in 1:length(eval_points)) {
    psi_hat_pnt[i,1] = eval_points[i]
    
    X_new <- obs[x_names]
    X_new[, colnames(X_new)=='A'] = eval_points[i]
    
    predictions <- predict(model_fit, newdata = X_new, type = "response", se.fit = TRUE)
    
    y_hat <-  mean(predictions$fit)
    SE <-  mean(predictions$se.fit)
    
    psi_hat_pnt[i,2] = y_hat
    psi_hat_pnt[i,3] = SE
    psi_hat_pnt[i,4] = y_hat - 1.96 * SE
    psi_hat_pnt[i,5] = y_hat + 1.96 * SE
  } 
  
  # bootstrap-based inference
  # psi_hat_pnt_bt <- bootstrap_inference_gam_poly(obs, eval_points, form, method)
  # psi_hat_pnt <- merge(psi_hat_pnt, psi_hat_pnt_bt, by = "a")

  return(psi_hat_pnt)
}


bootstrap_inference_gam_poly <- function(obs, eval_points, form, method){
  
  n_bootstrap <- 500  # Number of bootstrap iterations
  
  # Create an empty matrix to store bootstrap estimates
  bootstrap_estimates <- matrix(NA, nrow = length(eval_points), ncol = n_bootstrap)
  
  for (b in 1:n_bootstrap) {
    # Generate bootstrap sample by sampling with replacement from the original data
    bootstrap_data <- obs[sample(nrow(obs), replace = TRUE), ]
    
    if (method == 'GAM'){
      # Fit the GAM model on the bootstrap sample
      bootstrap_model <- gam(form, data = bootstrap_data, family = binomial())
    } else if (method == 'POLY'){
      # Fit a polynomial logistic regression model
      bootstrap_model <- glm(form, data = obs, family = binomial)
    }
    
    
    # Obtain predictions for the evaluation points
    for (i in 1:length(eval_points)) {
      X_new <- bootstrap_data[names(obs)[names(obs) != "Y"]]
      X_new[, colnames(X_new) == 'A'] <- eval_points[i]
      
      predictions <- predict(bootstrap_model, newdata = X_new, type = "response")
      bootstrap_estimates[i, b] <- mean(predictions)
    }
  }
  
  # Calculate the bootstrap point estimates and confidence intervals
  psi_hat_pnt_bt <- matrix(NA, nrow = length(eval_points), ncol = 4)
  colnames(psi_hat_pnt_bt) <- c("a", "SE_bt", "ci_lwr_bt", "ci_upr_bt")
  
  for (i in 1:length(eval_points)) {
    psi_hat_pnt_bt[i, 1] <- eval_points[i]
    bootstrap_sample <- bootstrap_estimates[i, ]
    
    psi_hat_pnt_bt[i, 2] <- sd(bootstrap_sample)
    psi_hat_pnt_bt[i, 3] <- quantile(bootstrap_sample, 0.025)
    psi_hat_pnt_bt[i, 4] <- quantile(bootstrap_sample, 0.975)
  }
  
  return(psi_hat_pnt_bt)
}



##############################################################
# This function runs the simulation for B rounds with given:
#       - data generating function
#       - sample size n
#       - number of simulations: B

# returns the estimated ATE, empirical 95% CI, and coverage rates
##############################################################
run_simu_gam_poly_rep <- function(gen_data_functions, eval_points, y_type, n, rounds, return_all_rslts = F, method = "GAM"){
  
  result_list <- list()
  
  for(r in 1:rounds){
    print(paste0("round ", r))
    result <- tryCatch({
      run_simu_1round_gam_poly(gen_data_functions, eval_points, y_type, n=n, method=method)
    }, error = function(e) {
      print(paste0("Error: ", e$message))
      NULL
    })
    
    while(is.null(result)) {
      print('retry with a new generated data')
      result <- tryCatch({
        run_simu_1round_gam_poly(gen_data_functions, eval_points, y_type, n=n, method=method)
      }, error = function(e) {
        print(paste0("Error: ", e$message))
        NULL
      })
    }
    
    result_list[[r]] <- result
  }
  
  result_all <-  do.call("rbind", result_list) %>% as.data.frame()
  result_all <- merge(as.data.frame(psi0_pnt), result_all, by=c("a"))
  
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
    mutate(method = method)
  
  if(return_all_rslts){
    results <- list(result_summary = result_summary,
                    all_results = result_list)
  } else {
    results <- list(result_summary = result_summary)
  }
  
  return(results)
}
