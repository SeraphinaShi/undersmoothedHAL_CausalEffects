#' Undersmoothed HAL
#'
#' Note: Current undersmoothed HAL use a global criterion.
#' Future work need to be done for user-specified criterion driven by the target parameter.
#'
#'
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
                            criterion = NULL,
                            a = NULL,
                            Nlam = 20,
                            family = "gaussian"){
  
  preds_init <- predict(fit_init, new_data = X)
  # estimates of sd in each direction using initial fit
  resid_init <- preds_init - Y
  sd_est  <- apply(basis_mat, 2, function(u) sd(resid_init*u))
  
  # refit on new lambda sequence
  us_lambda <- fit_init$lambda_star*10^seq(from=0, to=-3, length=Nlam)
  us_fit <- glmnet(fit_init$x_basis, Y, lambda=us_lambda, family = family, standardize = FALSE)
  # >   us_fit <- glmnet(fit_init$x_basis, Y, lambda=us_lambda, family = family, standardize = FALSE)
  # Warning messages:
  #   1: from glmnet C++ code (error code -1); Convergence for 1th lambda value not reached after maxit=100000 iterations; solutions for larger lambdas returned 
  # 2: In getcoef(fit, nvars, nx, vnames) :
  #   an empty model has been returned; probably a convergence issue

  
  if(identical(us_fit$df, 0)){
    res <- list("lambda_init" = fit_init$lambda_star,
                "lambda_under" = fit_init$lambda_star,
                "spec_under" = NA)
    return(res)
  }
  
  # check the criterion (global)
  # TBD user-specified criterion (e.g. target parameter driven)
  if (is.null(criterion)){
    
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
  } 
  
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




########################
# calculating efficient influence curves
########################
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


########################
# Auto-undersmoothing HAL
########################

AU_HAL_fit <- function(X,Y,
                       family='binomial',
                       smoothness_orders = 1){
  
  #================================CV-HAL================================
  hal_fit <- fit_hal(X = X, Y = Y, family = y_type,
                     return_x_basis = TRUE,
                     num_knots = hal9001:::num_knots_generator(
                       max_degree = ifelse(ncol(X) >= 20, 2, 3),
                       smoothness_orders = 1,
                       base_num_knots_0 = 20,
                       base_num_knots_1 = 20 # max(100, ceiling(sqrt(n)))
                     )
  )
  lambda_CV <- hal_fit$lambda_star
  
  #================================GU-HAL================================
  CV_nonzero_col <- which(hal_fit$coefs[-1] != 0)
  
  if (length(CV_nonzero_col) != 0){
    
    CV_basis_mat <- as.matrix(hal_fit$x_basis)
    CV_basis_mat_nonzero <- as.matrix(CV_basis_mat[, CV_nonzero_col])
    
    hal_undersmooth <- undersmooth_hal(X, Y,
                                       fit_init = hal_fit,
                                       basis_mat = CV_basis_mat,
                                       family = y_type)
    lambda_u_g = hal_undersmooth$lambda_under
    
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
    if(any(hal_u_g$coefs[-1] != 0)){
      hal_fit = hal_u_g
    }
  }
  
  return(hal_fit)
}



