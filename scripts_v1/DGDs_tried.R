###############################################################################
# generate_data1
# positivity violations between W and A, try generate_data2
###############################################################################

generate_data1 <- function(n, a=NA, z=NA){
  # exogenous variables
  U_W <- runif(n, 0, 1)
  U_A <- rnorm(n,0,1)
  U_Z <- runif(n, 0, 1)
  U_Y <- runif(n, 0, 1)
  
  
  # endogenous variables
  W <- U_W
  
  if(is.na(a)){
    A <-  10 - 10*W + U_A
    A[A<=0] = 0
    A[A>=10] = 10
  } else {
    A <- rep(a, n)
  }
  
  if(is.na(z)){
    Z <- as.numeric(U_Z < plogis(0.2*(W+A)))
  } else {
    Z <- rep(z, n)
  }
  
  Y <- as.numeric(U_Z < plogis(0.05 * (-W + 5*A - 3 * W * A - Z)))
  
  
  # data frame
  O <- data.frame(W, A, Z, Y)
  return(O)
}


###############################################################################
# generate_data2
# no positivity violations, but glm cannot fit this one
# not sure why, need to dig deeper
###############################################################################
generate_data2 <- function(n, a=NA, z=NA){
  # exogenous variables
  U_W <- runif(n, 0, 1)
  U_A <- rnorm(n,0,4)
  U_Z <- runif(n, 0, 1)
  U_Y <- runif(n, 0, 1)
  
  
  # endogenous variables
  W <- U_W
  
  if(is.na(a)){
    A <-  10 - 10*W  + U_A
    A[A<=0] = 0
    A[A>=10] = 10
  } else {
    A <- rep(a, n)
  }
  
  if(is.na(z)){
    Z <- as.numeric(U_Z < plogis(0.1*(W+A)))
  } else {
    Z <- rep(z, n)
  }
  
  Y <- as.numeric(U_Y < plogis(0.01 * (-W + 5*A - 3 * W * A - Z)))
  
  
  # data frame
  O <- data.frame(W, A, Z, Y)
  return(O)
}


###############################################################################
# generate_data3
###############################################################################
generate_data3 <- function(n, a=NA, z=NA){
  # exogenous variables
  U_W <- rnorm(n, 0, 1)
  U_A <- rnorm(n, 0, 2)
  U_Z <- runif(n, 0, 1)
  U_Y <- runif(n, 0, 1)
  
  # endogenous variables
  W <- U_W
  
  if(is.na(a)){
    A <-  2 - 0.5*W + U_A
    A[A<=0] = 0
    A[A>=5] = 5
  } else {
    A <- rep(a, n)
  }
  
  if(is.na(z)){
    Z <- as.numeric(U_Z < plogis(2-W-A))
  } else {
    Z <- rep(z, n)
  }
  
  Y <- as.numeric(U_Y < plogis(W + 2*A + Z - 0.5 * W * A))
  
  # data frame
  O <- data.frame(W, A, Z, Y)
  return(O)
}

