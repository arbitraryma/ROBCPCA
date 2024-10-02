# This code gives the algorithm 2 in the paper 

library(stats)
library(r2spss)
library(dvmisc)
library(cellWise)
library(rrcov)

# the input ts is your MTS objects, it should be a list!
# and the k here is the number of clusters you want 

initi_label <- function(ts,k){
  n <- length(ts)
  # store the info of centroids
  max_iter <- 100
  mu <- list()
  sigma <- list()
  # start with a random centroid
  
  diff <- Inf
  idealprop <- rep(1 / k, k)
  prop <- rep(0,k)
  # t <- ifelse(k <= 5, 0.2, 0.10)
  t <- ifelse(k <= 5, ifelse(k == 2, 0.3, 0.2), 0.2)
  # robustly estimate the location and covariance matrix of the centroid
  iter <- 1
  while (diff > t | any(prop==0) && iter < max_iter ) {
    a <- sample(n,1)
  est <- CovMrcd(ts[[a]],alpha = 0.75)
  mu[[1]] <- est@center
  sigma[[1]] <- est@cov
  
  # with the fist centroid, we compute the mahalanobis distance
  # notice, for a multivariate time series, it's mahalanobis distance to the chosen centroid is a vector
  # hence, we use a metric mean to represent the vector, somehow, it happens that a MTS object can contain
  # outliers, hence, we use trimmed mean to make the result more resistant to outliers
  
  dist <- matrix(0, n, k)

  # use the same strategy     
    for (m in 2:k) {
      dist[, m-1] <- sapply(1:n, function(i) {
        min(sapply(1:(m-1), function(j) {
          mean(trim(mahalanobis(ts[[i]], mu[[j]], sigma[[j]]), 0.2, tails = 'upper'))
          # mean(trim(mahalanobis(ts[[i]], mu[[j]], sigma[[j]]), 0.1, tails = 'both'))
        }))
      })
      
      sqdist <- dist[, m - 1]^2
      prob <- sqdist / sum(sqdist)
      next_centroid_idx <- sample(1:n, 1, prob = prob)
      est <- CovMrcd(ts[[next_centroid_idx]],alpha = 0.75)
      mu[[m]] <- est@center
      sigma[[m]] <- est@cov
    }
    
    # Ensure the last column is also filled
    dist[, k] <- sapply(1:n, function(i) {
      min(sapply(1:k, function(j) {
       
        mean(trim(mahalanobis(ts[[i]], mu[[j]], sigma[[j]]), 0.2, tails = 'upper'))
      }))
    })
    
    # Assign points to the closest centroid
    assignment <- sapply(1:n, function(i) {
      which.min(sapply(1:k, function(j) {
        mean(trim(mahalanobis(ts[[i]], mu[[j]], sigma[[j]]), 0.2, tails = 'upper'))
       
      }))
    })
    
    # Check proportion criteria
    
    actual_props <- table(factor(assignment, levels = 1:k)) / n  # Ensure all groups are represented
    
    prop <- as.numeric(actual_props)  # Convert to numeric vector
    diff <- max(abs(prop - idealprop))
    iter <- iter +1
  }
  
  if (iter>= max_iter ) {
    # this functions is done only if it surpasses 100 times but still cannot have a balanced grouping, so just do 1 more time,
    # and make sure each group has at least 1 object
    while (any(prop==0) ) {
      
    a <- sample(n,1)
    est <- CovMrcd(ts[[a]],alpha = 0.75)
    mu[[1]] <- est@center
    sigma[[1]] <- est@cov
    
    # with the fist centroid, we compute the mahalanobis distance
    # notice, for a multivariate time series, it's mahalanobis distance to the chosen centroid is a vector
    # hence, we use a metric mean to represent the vector, somehow, it happens that a MTS object can contain
    # outliers, hence, we use trimmed mean to make the result more resistant to outliers
    
    dist <- matrix(0, n, k)
    
    # use the same strategy     
    for (m in 2:k) {
      dist[, m-1] <- sapply(1:n, function(i) {
        min(sapply(1:(m-1), function(j) {
          mean(trim(mahalanobis(ts[[i]], mu[[j]], sigma[[j]]), 0.2, tails = 'upper'))
          # mean(trim(mahalanobis(ts[[i]], mu[[j]], sigma[[j]]), 0.1, tails = 'both'))
        }))
      })
      
      sqdist <- dist[, m - 1]^2
      prob <- sqdist / sum(sqdist)
      next_centroid_idx <- sample(1:n, 1, prob = prob)
      est <- CovMrcd(ts[[next_centroid_idx]],alpha = 0.75)
      mu[[m]] <- est@center
      sigma[[m]] <- est@cov
    }
    
    # Ensure the last column is also filled
    dist[, k] <- sapply(1:n, function(i) {
      min(sapply(1:k, function(j) {
        mean(trim(mahalanobis(ts[[i]], mu[[j]], sigma[[j]]), 0.2, tails = 'upper'))
      
      }))
    })
    
    # Assign points to the closest centroid
    assignment <- sapply(1:n, function(i) {
      which.min(sapply(1:k, function(j) {
         mean(trim(mahalanobis(ts[[i]], mu[[j]], sigma[[j]]), 0.2, tails = 'upper'))
       
      }))
    })
    
    # Check proportion criteria
    
    actual_props <- table(factor(assignment, levels = 1:k)) / n  # Ensure all groups are represented
    
    prop <- as.numeric(actual_props)  # Convert to numeric vector
  }
  }
    return(assignment)
}    
    
  





