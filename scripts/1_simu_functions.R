#' Undersmoothed HAL
#'
#' Note: Current undersmoothed HAL use a global criterion.
#' Future work need to be done for user-specified criterion driven by the target parameter.
#'
#'

###############################################################################
#' initialize undersmoothed HAL
#'
#' @details fit a regular 0-order HAL(\code{smoothness_orders = 0}) with binning
#'  (\code{num_knots}), select the set of basis functions with non-zero coefs
#'  which will be used in following procedure. Note this is not the only valid initializing
#'  procedure, the user can define their own initialization.
#' @param X An input \code{matrix} with dimensions number of observations -by-
#'  number of covariates that will be used to derive the design matrix of basis
#'  functions.
#' @param Y A \code{numeric} vector of observations of the outcome variable.
#' @param family A \code{character} or a \code{\link[stats]{family}} object
#'  (supported by \code{\link[glmnet]{glmnet}}) specifying the error/link
#'  family for a generalized linear model. See more details in \code{fit_hal}.

undersmooth_init <- function(X, Y, family = "gaussian"){
  
  n <- length(Y)
  
  # get initial fit
  fit_init <- fit_hal(X = X,
                      Y = Y,
                      smoothness_orders = 0,
                      return_x_basis = TRUE,
                      family = family,
                      num_knots = hal9001:::num_knots_generator(
                        max_degree = ifelse(ncol(X) >= 20, 2, 3),
                        smoothness_orders = 0,
                        base_num_knots_0 = max(100, ceiling(sqrt(n)))
                      )
  )
  
  # select the non-zero directions/basis
  init_coef <-fit_init$coefs[-1]
  nonzero_col <- which(init_coef != 0)
  basis_mat <- as.matrix(fit_init$x_basis)
  basis_mat <- as.matrix(basis_mat[, nonzero_col])
  
  res <- list("fit_init" = fit_init,
              "basis_mat" = basis_mat)
  return(res)
}

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
                            Nlam = 20,
                            family = "gaussian"){
  
  preds_init <- predict(fit_init, new_data = X)
  # estimates of sd in each direction using initial fit
  resid_init <- preds_init - Y
  sd_est  <- apply(basis_mat, 2, function(u) sd(resid_init*u))
  
  # refit on new lambda sequence
  us_lambda <- fit_init$lambda_star*10^seq(from=0, to=-3, length=Nlam)
  us_fit <- glmnet(fit_init$x_basis, Y, lambda=us_lambda, family = family, standardize = FALSE)
  
  # evaluate refits
  if (family != "binomial"){
    pred_mat <- predict(us_fit, fit_init$x_basis)
  }else {
    pred_mat <- predict(us_fit, fit_init$x_basis, type = "response")
  }
  resid_mat <- pred_mat - Y
  
  # check the criterion (global)
  # TBD user-specified criterion (e.g. target parameter driven)
  if (is.null(criterion)){
    max_score <- get_maxscore(basis_mat = basis_mat,
                              resid_mat = resid_mat,
                              sd_est = sd_est,
                              Nlam = Nlam, us_fit = us_fit)
    
    # get the first lambda that satisfies the criteria
    lambda_under <- us_lambda[max_score <= 1/(sqrt(n)*log(n))][1]
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
  max_score <- apply(abs(score_all), 1, max)
  return(max_score)
}




########################
# calculating ate
########################
cal_ate_hal <- function(df, y_name, x_names, y_type, a1, a0, z){
  n <- nrow(df)
  Y <- as.numeric(as.matrix(df %>% select(all_of(y_name))))
  X <- df %>% 
    select(all_of(x_names)) %>% 
    mutate_if(sapply(., is.factor), as.numeric)
  
  hal_init <- undersmooth_init(X, Y, family = y_type)
  hal_undersmooth <- undersmooth_hal(X, Y,
                                     fit_init = hal_init$fit_init,
                                     basis_mat = hal_init$basis_mat,
                                     family = y_type)

  
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

