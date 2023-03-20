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
run_simu_1round <- function(gen_data_functions, n, lambda_scaler = 1, undersmooth=F){
  obs <- gen_data_functions(n)
  
  y_name = "Y"
  x_names = c("W", "A", "Z")
  y_type = "binomial"
  
  Y <- as.numeric(as.matrix(obs %>% select(all_of(y_name))))
  X <- obs %>% 
    select(all_of(x_names)) %>% 
    mutate_if(sapply(., is.factor), as.numeric)
  
  # fitting HAL
  CV_hal <- fit_hal(X = X, Y = Y, family = y_type,
                    return_x_basis = TRUE,
                    num_knots = hal9001:::num_knots_generator(
                      max_degree = ifelse(ncol(X) >= 20, 2, 3),
                      smoothness_orders = 1,
                      base_num_knots_0 = 20,
                      base_num_knots_1 = 20 # max(100, ceiling(sqrt(n)))
                      )
  )
  CV_lambda <- CV_hal$lambda_star
  
  if((!undersmooth) & lambda_scaler == 1){
    hal_fit <- CV_hal
    lambda = CV_lambda
  } else {
    if(undersmooth){
      
      CV_nonzero_col <- which(CV_hal$coefs[-1] != 0)
      CV_basis_mat <- as.matrix(CV_hal$x_basis)
      CV_basis_mat <- as.matrix(CV_basis_mat[, CV_nonzero_col])
      
      hal_undersmooth <- undersmooth_hal(X, Y,
                                         fit_init = CV_hal,
                                         basis_mat = CV_basis_mat,
                                         family = y_type)
      lambda = hal_undersmooth$lambda_under
      
      while(is.na(lambda)){
        print("  Since the learned undersmoothed lambda is NA, refitting lambdas.")
        
        CV_hal <- fit_hal(X = X, Y = Y, family = y_type,
                          return_x_basis = TRUE,
                          num_knots = hal9001:::num_knots_generator(
                            max_degree = ifelse(ncol(X) >= 20, 2, 3),
                            smoothness_orders = 1,
                            base_num_knots_0 = 20,
                            base_num_knots_1 = 20 # max(100, ceiling(sqrt(n)))
                          )
        )
        CV_lambda <- CV_hal$lambda_star
        
        CV_nonzero_col <- which(CV_hal$coefs[-1] != 0)
        CV_basis_mat <- as.matrix(CV_hal$x_basis)
        CV_basis_mat <- as.matrix(CV_basis_mat[, CV_nonzero_col])
        
        hal_undersmooth <- undersmooth_hal(X, Y,
                                           fit_init = CV_hal,
                                           basis_mat = CV_basis_mat,
                                           family = y_type)
        lambda = hal_undersmooth$lambda_under
      }
      
      
    } else {
      lambda = lambda_scaler * CV_lambda
    }
    
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
                       lambda = lambda
    )

  } 
  
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
  
  if(any(is.na(IC_beta))){
    # print("failed to calculate IC_beta because of the given basis matrix is not invertible")
    return(NA)
  }
  
  # calculate ATE
  psi_hat_10pnt = matrix(ncol = 6)
  for (z in c(1,0)) {
    for (a in seq(0.5,5,0.5)) {
      X_new <- X
      X_new$A = a
      X_new$Z = z
      
      X_new_0 <- X_new
      X_new_0$A = 0
      
      Ya_hat <- predict(hal_fit, new_data = X_new)
      Y0_hat <- predict(hal_fit, new_data = X_new_0)
      psi_hat <- mean(Ya_hat - Y0_hat)
      
      # efficient influence curve
      x_basis_a <- make_design_matrix(as.matrix(X_new), hal_fit$basis_list, p_reserve = 0.75)
      x_basis_a_nonzero <- as.matrix(cbind(1, x_basis_a)[, nonzero_idx])
      
      x_basis_0 <- make_design_matrix(as.matrix(X_new_0), hal_fit$basis_list, p_reserve = 0.75)
      x_basis_0_nonzero <- as.matrix(cbind(1, x_basis_0)[, nonzero_idx])
      
      IC_ATE <- cal_IC_for_ATE(X_new_a = x_basis_a_nonzero, 
                               X_new_0 = x_basis_0_nonzero, 
                               beta_n = coef_nonzero, IC_beta = IC_beta)
      
      # empirical SE and 95% confidence interval
      
      SE <- sqrt(var(IC_ATE)/n)
      ci_lwr <- psi_hat - 1.96 * SE
      ci_upr <- psi_hat + 1.96 * SE
      
      psi_hat_10pnt <- rbind(psi_hat_10pnt, matrix(c(a, z, psi_hat, SE, ci_lwr, ci_upr), nrow = 1))
    }   
  }
  psi_hat_10pnt <- psi_hat_10pnt[-1, ]
  colnames(psi_hat_10pnt) <- c("a", "z", "psi_hat", "SE", "ci_lwr", "ci_upr")
  
  psi_hat_10pnt <- cbind(psi_hat_10pnt, lambda, CV_lambda, lambda_scaler)
  
  return(psi_hat_10pnt)
}


##############################################################
# This function runs the simulation for B rounds with given:
#       - data generating function
#       - sample size n
#       - number of simulations: B

# returns the estimated ATE, empirical 95% CI, and coverage rates
##############################################################
run_simu_rep <- function(gen_data_functions, n, B, lambda_scaler=1, undersmooth=F, return_all_rslts = F){
  result_list <- list()
  for(b in 1:B){
    result <- run_simu_1round(gen_data_functions, n=n, lambda_scaler, undersmooth)
    while(any(is.na(result))){
      result <- run_simu_1round(gen_data_functions, n=n, lambda_scaler, undersmooth)
    }
    result_list[[b]] <- result
  }
  result_all <-  do.call("rbind", result_list) %>% as.data.frame()
  result_all <- merge(as.data.frame(psi0_10pnt), result_all, by=c("a", "z"))
  
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
    summarise(across(everything(), mean)) %>% 
    ungroup()
  
  if(return_all_rslts){
    return(list(result_summary = result_summary,
                all_results = result_list))
  } else {
    return(result_summary)
  }
}



##############################################################
# This function runs the simulation for B rounds with given:
#       - data generating function
#       - sample size n
#       - number of simulations: B

# returns the estimated ATE, empirical 95% CI, and coverage rates
##############################################################
plot_perforences_1lambda_alla <- function(df, z_para=1, est_plot_only=F, plot_list=F, add_oracal=F){
  
  df <- df %>% filter(z==z_para)
  
  color = ifelse(z_para==1, "#F8766D", "#00BFC4")
  oracal_color = ifelse(z_para==1,"#00BA38", "#F564E3")

  p_est_avg <- ggplot(data=df, aes(x=a)) +
    geom_line(aes(y=psi0), alpha = 0.5, color="darkgrey") +
    geom_errorbar(aes(ymin=ci_lwr, ymax=ci_upr, color='Empirical'), width=0.7) +
    geom_point(aes(y=psi0), color = "black") +
    geom_point(aes(y=psi_hat, color='Empirical'), shape=17, size=2) +
    labs(x="a", y="ATE", title = "Estimation") +
    scale_color_manual(name='Method',
      breaks=c('Empirical', 'Oracal'),
      values=c('Empirical'=color, 'Oracal'=oracal_color)) +
    scale_x_continuous(limits = c(0, 5), breaks = 0:5)

  if(add_oracal){
    p_est_avg <- p_est_avg + 
      geom_errorbar(aes(ymin=oracal_ci_lwr, ymax=oracal_ci_upr, color='Oracal'), width=0.7) 
  } 
  
  legend <- get_legend(p_est_avg)
  p_est_avg <- p_est_avg + theme(legend.position='none')
  
  if(est_plot_only){
    return(p_est_avg + labs(title=""))
  } else {
    p_bias <- ggplot(df, aes(x = a, y = bias)) +  
      geom_line(color=color) +
      geom_point(color=color) + 
      labs(title="|Bias|") 
    
    if(!add_oracal){
      p_se <- ggplot(df, aes(x = a, y = SE)) +  
        geom_line(color=color) +
        geom_point(color=color) + 
        labs(title="Standard Error") 
      
      p_bias_d_df <- ggplot(df, aes(x = a, y = bias_se_ratio)) +  
        geom_line(color=color) + 
        geom_point(color=color) + 
        labs(title="|Bias| / Standard Error") 
      
      p_cr <- ggplot(df, aes(x = a, y = cover_rate)) +  
        geom_line(color=color) + 
        geom_point(color=color) +
        labs(title="95% CI Coverage Rate")
    } else {
      p_se <- ggplot(df, aes(x = a)) +  
        geom_line(aes(y = SE, color='Empirical'),alpha=0.7) +
        geom_point(aes(y = SE, color='Empirical'),alpha=0.7) + 
        geom_line(aes(y = oracal_SE, color='Oracal'),alpha=0.7) +
        geom_point(aes(y = oracal_SE, color='Oracal'),alpha=0.7) + 
        scale_color_manual(name='Method',
          breaks=c('Empirical', 'Oracal'),
          values=c('Empirical'=color, 'Oracal'=oracal_color)) + 
        theme(legend.position='none') +
        labs(title="Standard Error") 
      
      
      p_bias_d_df <- ggplot(df, aes(x = a)) +  
        geom_line(aes(y = bias_se_ratio, color='Empirical'), alpha=0.7) +
        geom_point(aes(y = bias_se_ratio, color='Empirical'), alpha=0.7) + 
        geom_line(aes(y = oracal_bias_se_ratio, color='Oracal'), alpha=0.7) +
        geom_point(aes(y = oracal_bias_se_ratio, color='Oracal'), alpha=0.7) + 
        scale_color_manual(name='Method',
          breaks=c('Empirical', 'Oracal'),
          values=c('Empirical'=color, 'Oracal'=oracal_color)) + 
        theme(legend.position='none') +
        labs(title="|Bias| / Standard Error") 
      
      p_cr <- ggplot(df, aes(x = a)) +  
        geom_line(aes(y = cover_rate, color='Empirical'), alpha=0.7) +
        geom_point(aes(y = cover_rate, color='Empirical'), alpha=0.7) + 
        geom_line(aes(y = oracal_cover_rate, color='Oracal'), alpha=0.7) +
        geom_point(aes(y = oracal_cover_rate, color='Oracal'), alpha=0.7) + 
        scale_color_manual(name='Method',
          breaks=c('Empirical', 'Oracal'),
          values=c('Empirical'=color, 'Oracal'=oracal_color)) + 
        theme(legend.position='none') +
        labs(title="95% CI Coverage Rate")
    }
    
    
    if(plot_list){
      return(list(p_est_avg, p_bias, p_se, p_bias_d_df, p_cr))
    } else {
      p <- grid.arrange(p_est_avg, p_bias, p_se, p_bias_d_df, p_cr, legend,
                          layout_matrix = rbind(c(NA,1,1,NA),
                                                c(NA,1,1,6),
                                                c(2,2,3,3),
                                                c(2,2,3,3),
                                                c(4,4,5,5),
                                                c(4,4,5,5)),
                          top = textGrob(paste0("HAL-based plug in estimator performence for ATE"), 
                                         gp=gpar(fontsize=11, fontface = 'bold')))
        return(p)
    }
  }

}


estimation_qqplot <- function(results_list, z_para){
  df <- do.call("rbind", results_list) %>% 
    as.data.frame() %>% filter(z == z_para)
  
  p <- ggplot(df, aes(sample = psi_hat)) + 
    stat_qq() + stat_qq_line() +
    # theme(axis.title.x=element_blank(), 
    #       axis.title.y=element_blank()) +
    labs(x = "Theoretical Quantiles",
         y = "Sample Quantiles",
         title = paste0("Normal Q-Q Plot for Estimated ATE when z=",z_para)) +
    facet_wrap(~a, nrow=2) 
  
  return(p)
}



plot_perforences_alllambda_1a <- function(df, a_para, z_para, add_oracal=F){
  df <- df %>% filter(a == a_para, z == z_para)
  
  p_est_avg <- ggplot(df) +  
    geom_line(aes(x = lambda_scaler, y = psi_hat), color = "grey") + 
    geom_point(aes(x = lambda_scaler, y = psi_hat)) + 
    geom_hline(aes(yintercept=psi0)) +
    labs(title="Estimation average") +
    theme()
  
  p_bias <- ggplot(df, aes(x = lambda_scaler, y = bias)) +  
    geom_line(color = "grey") + 
    geom_point() + 
    labs(title="Bias") 
  
  p_se <- ggplot(df, aes(x = lambda_scaler)) +  
    geom_line(aes(y = SE), color = "grey") + 
    geom_point(aes(y = SE, color = "Empirical")) + 
    labs(title="Standard Error") + 
    scale_color_manual(name='method',
                       breaks=c('Empirical', 'Oracal'),
                       values=c('Empirical'='black', 'Oracal'='#8B6508'))

  
  p_bias_d_df <- ggplot(df, aes(x = lambda_scaler)) + 
    geom_line(aes(y = bias_se_ratio), color = "grey") + 
    geom_point(aes(y = bias_se_ratio, color = "Empirical")) + 
    labs(title="|Bias| / Standard Error") + 
    scale_color_manual(name='method',
                       breaks=c('Empirical', 'Oracal'),
                       values=c('Empirical'='black', 'Oracal'='#8B6508')) + 
    theme(legend.position='none') 
  
  p_cr <- ggplot(df, aes(x = lambda_scaler)) +  
    geom_line(aes(y = cover_rate), color = "grey") + 
    geom_point(aes(y = cover_rate, color = "Empirical")) + 
    labs(title="Coverage rate") + 
    scale_color_manual(name='method',
                       breaks=c('Empirical', 'Oracal'),
                       values=c('Empirical'='black', 'Oracal'='#8B6508')) + 
    theme(legend.position='none')  
  
  if(!add_oracal){
    legend <- get_legend(p_se)
    p_se <- p_se + theme(legend.position='none')
  } else {
    p_se <- p_se + 
      geom_line(aes(y = oracal_SE), color = "#EEC591") + 
      geom_point(aes(y = oracal_SE, color = "Oracal")) 
    
    legend <- get_legend(p_se)
    p_se <- p_se + theme(legend.position='none')
    
    p_bias_d_df <- p_bias_d_df +  
      geom_line(aes(y = oracal_bias_se_ratio), color = "#EEC591") + 
      geom_point(aes(y = oracal_bias_se_ratio, color = "Oracal")) 
    
    p_cr <- p_cr +  
      geom_line(aes(y = oracal_cover_rate), color = "#EEC591") + 
      geom_point(aes(y = oracal_cover_rate, color = "Oracal")) 
  }
  
  p <- grid.arrange(p_est_avg, p_bias, p_se, p_bias_d_df, p_cr, legend, 
                    layout_matrix = rbind(c(NA,1,1,NA),
                                          c(NA,1,1,6),
                                          c(2,2,3,3),
                                          c(2,2,3,3),
                                          c(4,4,5,5),
                                          c(4,4,5,5)),
                    top = textGrob(paste0("HAL-based plug in estimator performence for a=", a_para, ", z=", z_para), 
                                   gp=gpar(fontsize=11, fontface = 'bold')))
  return(p)
}

