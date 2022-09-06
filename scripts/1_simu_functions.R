###############################################################################
## Undersoomthed HAL function
## ----------------------------------------------------------------------------
## inputs: 
##       df:            a dataframe
##       yname:         a string  
##       xname:         a vector of string
##       Nlam:          a scalar, number of candidates lambda
##       family_type:   a string, type of y, ("binomial", "gaussian")
###############################################################################

undermoothed_HAL <- function(df, y_name, x_name, Nlam, family_type = "binomial"){
  
  # variables
  y <- as.numeric(as.matrix(df %>% select(all_of(yname))))
  x <- df %>% 
    select(all_of(xname)) %>% 
    mutate_if(sapply(., is.factor), as.numeric)
  n <- nrow(df)
  
  # get initial fit
  print("---initial fitting---")
  tic()
  fit <- fit_hal(X = x,
                 Y = y, 
                 return_x_basis = TRUE,
                 family = family_type,
                 num_knots = num_knots_generator(
                   max_degree = ifelse(ncol(x) >= 20, 2, 3),
                   smoothness_orders = 1,
                   base_num_knots_0 = 500,
                   base_num_knots_1 = max(100, ceiling(sqrt(n)))
                 )
  )
  toc()
  
  # only non-zero direction
  init_coef <-fit$coefs[-1]
  nonzero_col <- which(init_coef != 0)
  
  # if all coefs are zero, skip undersmooth and use the initial fit
  if (length(nonzero_col) == 0){
    res <- list("lambda_init" = fit$lambda_star,
                "lambda_under"= fit$lambda_star)
  }else{
    # refit on new lambda sequence
    print("---candidates fitting---")
    tic()
    us_lambda <- fit$lambda_star*10^seq(from=0, to=-3, length=Nlam)
    us_fit <- glmnet(fit$x_basis, y, lambda=us_lambda, family = type, standardize = FALSE)
    toc()
    
    # evaluate refits
    if (type == "gaussian"){
      pred_mat <- predict(us_fit, fit$x_basis)
      preds_init <- predict(fit, new_data = x)
    }else if (type == "binomial"){
      pred_mat <- predict(us_fit, fit$x_basis, type = "response")
      preds_init <- predict(fit, new_data = x) 
    }
    
    resid_mat <- pred_mat-y
    basis_mat <- as.matrix(fit$x_basis)
    basis_mat <- as.matrix(basis_mat[, nonzero_col])
    
    # estimates of sd in each direction using initial fit
    resid_init <- preds_init-y
    sd_est <- rep(NA, ncol(basis_mat))
    
    for (i in 1:ncol(basis_mat)){
      u <- basis_mat[,i]
      score_init <- resid_init * u
      sd_est[i] <- sd(score_init)
    }
    
    # check the criterion 
    print("---checking criterion---")
    tic()
    max_score <- get_maxscore(basis_mat = basis_mat, 
                              resid_mat = resid_mat,
                              sd_est = sd_est, 
                              Nlam = Nlam, us_fit = us_fit)
    toc()
    
    # get the first lambda that satisfies the criteria
    lambda_under <- us_lambda[max_score <= 1/(sqrt(n)*log(n))][1]
    
    # collect and save results
    coef_mat <- as.matrix(us_fit$beta)
    
    df_under <- data.frame("lambda" = NA,
                           "l1_norm" = NA,
                           "n_coef" = NA)
    
    for (j in 1:Nlam){
      df_under[j,1] = us_lambda[j]
      df_under[j,2] = sum(abs(coef_mat[,j]))
      df_under[j,3] = sum(coef_mat[,j] != 0)
    }
    
    res <- list("lambda_init" = fit$lambda_star,
                "lambda_under" = lambda_under,
                "df_under" = df_under)
  }
  return(res)
}