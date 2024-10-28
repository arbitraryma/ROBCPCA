# This function standardizes each multivariate time series (MTS) object by removing the mean of each variable (column-wise) 
# and possibly scaling by some measure of spread (e.g., standard deviation). 
standardize_mts <- function(mts_data) {
  # mts_data is a list where each element is an MTS object (a matrix)
  lapply(mts_data, function(ts) {
    scale(ts, center = TRUE, scale = FALSE) # Subtracts the column means
  })
}


# This function calculates the cross-covariance matrices for each MTS object in a list. 
# It constructs a block matrix including both the lag-0 covariance matrix and the lag-1 cross-covariance matrix:

cross_covariance_lag1 <- function(mts_data) {
  # mts_data is a multivariate time series (matrix or data frame)
  
  n <- nrow(mts_data) # number of time points
  p <- ncol(mts_data) # number of variables
  
  # Check if there are enough time points
  if (n < 2) {
    stop("The multivariate time series should have at least 2 time points.")
  }
  
  # Compute Gamma_0 (lag-0 covariance matrix), biased
  Gamma_0 <- cov(mts_data)
  
  # Compute Gamma_1 (lag-1 cross-covariance matrix), biased
  X_t <- mts_data[1:(n-1), ]    # Data from time t
  X_t_plus_1 <- mts_data[2:n, ]  # Data from time t+1
  
  # Cross-covariance between X_t and X_{t+1}, biased
  Gamma_1 <- cov(X_t, X_t_plus_1)
  
  # Construct the block matrix
  block_matrix <- rbind(
    cbind(Gamma_0, Gamma_1),   # First row: [Gamma(0), Gamma(1)]
    cbind(t(Gamma_1), Gamma_0) # Second row: [Gamma(-1), Gamma(0)]
  )
  
  output <- list(block_matrix, X_t, X_t_plus_1)
  return(output)
}



# this is the function when we consider lag-0 and lag-2 

cross_covariance_lag2 <- function(mts_data) {
  # mts_data is a multivariate time series (matrix or data frame)
  
  n <- nrow(mts_data) # number of time points
  p <- ncol(mts_data) # number of variables
  
  # Check if there are enough time points
  if (n < 3) { # Need at least 3 time points for lag 2
    stop("The multivariate time series should have at least 3 time points.")
  }
  
  # Compute Gamma_0 (lag-0 covariance matrix), biased
  Gamma_0 <- cov(mts_data)
  
  # Compute Gamma_2 (lag-2 cross-covariance matrix), biased
  X_t <- mts_data[1:(n-2), ]    # Data from time t
  X_t_plus_2 <- mts_data[3:n, ]  # Data from time t+2
  
  # Cross-covariance between X_t and X_{t+2}, biased
  Gamma_2 <- cov(X_t, X_t_plus_2)
  
  # Construct the block matrix with lag 2
  block_matrix <- rbind(
    cbind(Gamma_0, Gamma_2),   # First row: [Gamma(0), Gamma(2)]
    cbind(t(Gamma_2), Gamma_0) # Second row: [Gamma(-2), Gamma(0)]
  )
  
  output <- list(block_matrix, X_t, X_t_plus_2)
  return(output)
}



# this is the function when we consider lag-0 lag-1, lag-2 

cross_covariance_lag12 <- function(mts_data) {
  # mts_data is a multivariate time series (matrix or data frame)
  
  n <- nrow(mts_data) # number of time points
  p <- ncol(mts_data) # number of variables
  
  # Check if there are enough time points
  if (n < 3) { # Need at least 3 time points for lag 2
    stop("The multivariate time series should have at least 3 time points.")
  }
  

  # Compute Sigma(0) (lag-0 covariance matrix), biased
  Sigma_0 <- cov(mts_data)
  
  # Compute Sigma(1) (lag-1 cross-covariance matrix), biased
  X_t_lag1 <- mts_data[1:(n-1), ]    # Data from time t
  X_t_plus_1 <- mts_data[2:n, ]      # Data from time t+1
  Sigma_1 <- cov(X_t_lag1, X_t_plus_1)
  
  # Compute Sigma(2) (lag-2 cross-covariance matrix), biased
  X_t_lag2 <- mts_data[1:(n-2), ]    # Data from time t
  X_t_plus_2 <- mts_data[3:n, ]      # Data from time t+2
  Sigma_2 <- cov(X_t_lag2, X_t_plus_2)
  
  # Construct the 2p*2p block matrix
  combined_matrix <- 
    rbind(cbind(Sigma_0, Sigma_1),
    cbind(t(Sigma_1), Sigma_0)) +
    rbind(cbind(Sigma_0, Sigma_2),   # First row: [Gamma(0), Gamma(1)] + [Gamma(0), Gamma(2)]
    cbind(t(Sigma_2), Sigma_0))      # Second row: [Gamma(-1), Gamma(0)] + [Gamma(-2), Gamma(0)]
    
  
   X_t_plus_1 <- mts_data[2:(n-1), ] # extract the same length (n-1), so we can combine with when doing projection
  
  
  
  output <- list(combined_matrix, X_t_lag2, X_t_plus_1, X_t_plus_2)
  return(output)
}





combine_xt_xt1 <- function(tsxt, tsxt1) {
  # Combine tsxt and tsxt1 by columns
  combined_matrix <- cbind(tsxt, tsxt1)
  return(combined_matrix)
}




combine_lag12 <- function(tsxt,tsxt1,tsxt2){
  # Combine tsxt and tsxt1 by columns
  combined_matrix <- cbind(tsxt, tsxt1,tsxt2)
  return(combined_matrix)
}

