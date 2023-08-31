library(origami)
library(hal9001)

X <- cbind(rnorm(n = 20), rnorm(n = 20), rnorm(n = 20))
Y <- sample(c(0,1), size = 20, replace = T)

hal_fit = fit_hal(X, Y, family = "binomial")
summary(hal_fit)

predict(hal_fit, new_data = X)
predict(hal_fit, new_data = X, type = "response")
predict(hal_fit, new_data = X, type = "link")



# In the predict function: 
#   # return predictions if link function scale is acceptable
#   if (type == "link") {
#     # output predictions on the link function scale
#     return(preds)
#   }
# 
# # apply inverse family (link function) transformations
# if (inherits(object$family, "family")) {
#   inverse_link_fun <- object$family$linkinv
#   preds <- inverse_link_fun(preds)
# } else if (object$family == "binomial") {
#   preds <- stats::plogis(preds)
# } else if (object$family %in% c("poisson", "cox")) {
#   preds <- exp(preds)
# }
# 
# # bound predictions within observed outcome bounds if on response scale
# bounds <- object$prediction_bounds
# if (!is.null(bounds)) {
#   bounds <- sort(bounds)
#   if (is.matrix(preds)) {
#     preds <- apply(preds, 2, pmax, bounds[1])
#     preds <- apply(preds, 2, pmin, bounds[2])
#   } else {
#     preds <- pmax(bounds[1], preds)
#     preds <- pmin(preds, bounds[2])
#   }
# }
# 
# # output predictions on the response scale
# return(preds)