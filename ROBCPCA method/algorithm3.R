library(pracma)
library(base)
library(fungible)


# This function is Algorithm 3, this is with lag 0

robcpca <-  function(ts,k,p=NULL)
{
  
  # the initial step 
  n <- length(ts) # number of MTS  observations
  k <- k   # number of clusters 
  # l <- c(rep(1/k,k))  # initial probablity in each cluster to have the common projection axes, evenly distributed
  dim <- ncol(ts[[1]])  # for cpca, we assume all ts are with equal dimension
  ts1 <- list()
  sigma <- list()
  m <- list()
  # first adaption, using the robust covariance estimator
  for (i in 1:n) {
    
    m[[i]] <- colMeans(ts[[i]])
    ts1[[i]] <- sweep(ts[[i]], 2, m[[i]], `-`)
    sigma[[i]] <- cov(ts1[[i]])
    sigma[[i]] <- sigma[[i]] + 0.01* eye(dim)
  }
  
  group <- list()
  s <- list()
  
  indx <- initi_label(ts,k)
  
  
  for (j in 1:k) {
    # in the original paper, they initially average the number l <- n/k, somehow they may cannot
    # be divided evenly. Here we just consider to devide them nearly equally
    group[[j]] <- sigma[indx==j]
    s[[j]] <- comaxe(group[[j]])
  }
  
  # the refinement step
  t <- 0
  tmax <- 100
  toterror <- rep(0,tmax+1) # this is, the total projection error
  error <-  matrix(rep(0,n*k), ncol = k) # this is, the projection error for each observation to each cluster
  
  while (t <= tmax | length(unique(indx))!=k) {
    t <- t+1 
    # indx <- NULL
    for (v in 1:n) { # for each observation
      for (l in 1:k) { # for each cluster
        y <- ts1[[v]]%*%s[[l]][[1]] %*% t(s[[l]][[1]]) # the projection 
        df <- ts1[[v]] - y
        error[v, l] <- normF(df)
      }
    }
    
    # a contains the minimum error  
    a <- apply(error, 1, function(x) {minerror = min(x)})
    # b contains the the group index for which it has the minim error
    indx <- apply(error, 1, function(x) {cluster_index = which.min(x) })
    
    toterror[t+1] <- sum(a)
    
    if (abs(toterror[t+1]-toterror[t])<.Machine$double.eps ) {
      break
    }
  } 
  # this following while is set to guarantee that for the initilization step, each cluster can have at least one MTS object
    while(length(unique(indx))!=k) {
      indx <- sample.int(k,n,replace = TRUE,prob = rep(1/k,k))
    }
    # renew the cluster 
    for (r in 1:k) {
      group[[r]] <- sigma[indx==r]
      s[[r]] <- comaxe(group[[r]])
    }    
    
 
  return(indx)
}




# this function considers with only lag one, and we would like to test whether adding one lag can improve the accuracy 
robcpca_lag1 <- function(ts,k,p=NULL)
{
  
  # the initial step 
  n <- length(ts) # number of MTS  observations
  p <- ncol(ts[[1]])
  k <- k   # number of clusters 
  l <- c(rep(1/k,k))  # initial probobility in each cluster to have the common projection axes, evenly distributed
  dim <- ncol(ts[[1]])  # for cpca, we assume all ts are with equal dimension
  ts1 <- list()
  sigma <- list()
  
  # for the robust version, let's also test whether this part we can use the estlocscale to do 
  ts1 <- standardize_mts(ts)
  
  # Apply the cross_covariance function to each MTS object
  results <- lapply(ts1, cross_covariance_lag1)
  
  # Extract the block matrices (sigma), X_t (tsxt), and X_t_plus_1 (tsxt1) into separate lists
  
  sigma <- lapply(results, function(x) x[[1]]) # with dimension 2p * 2p
  sigma <- lapply(sigma, function(mat) {
    mat + diag(0.01, nrow = nrow(mat), ncol = ncol(mat))  # Add 0.01 to the diagonal
  })
  
  tsxt <- lapply(results, function(x) x[[2]])  # with dimension (n-1) * p
  tsxt1 <- lapply(results, function(x) x[[3]])
  
  combined_list <- mapply(combine_xt_xt1, tsxt, tsxt1, SIMPLIFY = FALSE) # this gives dimension (n-1)*2p
  
  group <- list()
  s <- list()
  indx <- initi_label(ts,k)
  
  for (j in 1:k) {
    # in the original paper, they initially average the number l <- n/k, somehow they may cannot
    # be divided evenly. Here we just consider to devide them nearly equally
    group[[j]] <- sigma[indx==j]
    s[[j]] <- comaxe(group[[j]])
  }
  
  # the refinement step
  t <- 0
  tmax <- 100
  toterror <- rep(0,tmax+1) # this is, the total projection error
  error <-  matrix(rep(0,n*k), ncol = k) # this is, the projection error for each observation to each cluster
  
  while (t <= tmax) {
    t <- t+1 
    # indx <- NULL
    for (m in 1:n) { # for each observation

      for (l in 1:k) { # for each cluster
        y <- combined_list[[m]]%*%s[[l]][[1]] %*% t(s[[l]][[1]]) # the projection 
        df <- combined_list[[m]] - y
        error[m, l] <- normF(df)
      }
    }
    
    # a contains the minimum error  
    a <- apply(error, 1, function(x) {minerror = min(x)})
    # b contains the the group index for which it has the minim error
    indx <- apply(error, 1, function(x) {cluster_index = which.min(x) })
    
    toterror[t+1] <- sum(a)
    
    if (abs(toterror[t+1]-toterror[t])<.Machine$double.eps ) {
      break
    }
    
    # this helps to prevent the super bad initialization where we may
    # encounter a group without any MTS object, thus meet an error later
    while(length(unique(indx))!=k) {
      indx <- sample.int(k,n,replace = TRUE,prob = rep(1/k,k))
    }
    # renew the cluster 
    for (r in 1:k) {
      group[[r]] <- sigma[indx==r]
      s[[r]] <- comaxe(group[[r]])
    }    
    
  } 

  return(indx)
 
}





# this function considers with only lag zero and two
# and we would like to test whether adding one lag can improve the accuracy 
robcpca_lag_2 <- function(ts,k,p=NULL)
{
  
  # the initial step 
  n <- length(ts) # number of MTS  observations
  p <- ncol(ts[[1]])
  k <- k   # number of clusters 
  l <- c(rep(1/k,k))  # initial probability in each cluster to have the common projection axes, evenly distributed
  dim <- ncol(ts[[1]])  # for cpca, we assume all ts are with equal dimension
  ts1 <- list()
  sigma <- list()
  
  # for the robust version, let's also test whether this part we can use the estlocscale to do 
  ts1 <- standardize_mts(ts)
  
  # Apply the cross_covariance function to each MTS object
  results <- lapply(ts1, cross_covariance_lag2)
  
  # Extract the block matrices (sigma), X_t (tsxt), and X_t_plus_1 (tsxt1) into separate lists
  
  sigma <- lapply(results, function(x) x[[1]]) # with dimension 2p * 2p
  sigma <- lapply(sigma, function(mat) {
    mat + diag(0.01, nrow = nrow(mat), ncol = ncol(mat))  # Add 0.01 to the diagonal
  })
  
  tsxt <- lapply(results, function(x) x[[2]])  # with dimension (n-2) * p
  tsxt1 <- lapply(results, function(x) x[[3]])
  
  combined_list <- mapply(combine_xt_xt1, tsxt, tsxt1, SIMPLIFY = FALSE) # this gives dimension (n-1)*2p
  
  group <- list()
  s <- list()
  indx <- initi_label(ts,k)
  
  for (j in 1:k) {
    # in the original paper, they initially average the number l <- n/k, somehow they may cannot
    # be divided evenly. Here we just consider to devide them nearly equally
    group[[j]] <- sigma[indx==j]
    s[[j]] <- comaxe(group[[j]])
  }
  
  # the refinement step
  t <- 0
  tmax <- 100
  toterror <- rep(0,tmax+1) # this is, the total projection error
  error <-  matrix(rep(0,n*k), ncol = k) # this is, the projection error for each observation to each cluster
  
  while (t <= tmax) {
    t <- t+1 
    
    for (m in 1:n) { # for each observation
      
      for (l in 1:k) { # for each cluster
        y <- combined_list[[m]]%*%s[[l]][[1]] %*% t(s[[l]][[1]]) # the projection 
        df <-combined_list[[m]] - y
        error[m, l] <- normF(df)
      }
    }
    
    # a contains the minimum error  
    a <- apply(error, 1, function(x) {minerror = min(x)})
    # b contains the the group index for which it has the minim error
    indx <- apply(error, 1, function(x) {cluster_index = which.min(x) })
    
    toterror[t+1] <- sum(a)
    
    if (abs(toterror[t+1]-toterror[t])<.Machine$double.eps ) {
      break
    }
    
    # this helps to prevent the super bad initialization where we may
    # encounter a group without any MTS object, thus meet an error later
    while(length(unique(indx))!=k) {
      indx <- sample.int(k,n,replace = TRUE,prob = rep(1/k,k))
    }
    # renew the cluster 
    for (r in 1:k) {
      group[[r]] <- sigma[indx==r]
      s[[r]] <- comaxe(group[[r]])
    }    
    
  } 
  return(indx)
}








# this is a trail for lag1 + lag 2, sum error 
# this function considers with lag 0, 1 and two, and we would like to test whether adding one lag can improve the accuracy 
robcpca_lag_12 <- function(ts,k,p=NULL)
{
  
  # the initial step 
  n <- length(ts) # number of MTS  observations
  p <- ncol(ts[[1]])
  k <- k   # number of clusters 
  l <- c(rep(1/k,k))  # initial probobility in each cluster to have the common projection axes, evenly distributed
  dim <- ncol(ts[[1]])  # for cpca, we assume all ts are with equal dimension
  ts1 <- list()
  sigma <- list()
  
  # for the robust version, let's also test whether this part we can use the estlocscale to do 
  ts1 <- standardize_mts(ts)
  
  # Apply the cross_covariance function to each MTS object
  results <- lapply(ts1, cross_covariance_lag1)
  results2 <-  lapply(ts1, cross_covariance_lag2)
  
  # Extract the block matrices (sigma), X_t (tsxt), and X_t_plus_1 (tsxt1) into separate lists
  
  sigma <- lapply(results, function(x) x[[1]]) # with dimension 2p * 2p
  sigma <- lapply(sigma, function(mat) {
    mat + diag(0.01, nrow = nrow(mat), ncol = ncol(mat))  # Add 0.01 to the diagonal
  })
  
  tsxt <- lapply(results, function(x) x[[2]])  # with dimension (n-1) * p
  tsxt1 <- lapply(results, function(x) x[[3]])
  
  combined_list <- mapply(combine_xt_xt1, tsxt, tsxt1, SIMPLIFY = FALSE) # this gives dimension (n-1)*2p
  
  
  sigma2 <- lapply(results2, function(x) x[[1]]) # with dimension 2p * 2p
  sigma2 <- lapply(sigma2, function(mat) {
    mat + diag(0.01, nrow = nrow(mat), ncol = ncol(mat))  # Add 0.01 to the diagonal
  })
  
  tsxt2 <- lapply(results2, function(x) x[[2]])  # with dimension (n-2) * p
  tsxt12 <- lapply(results2, function(x) x[[3]])
  combined_list2 <- mapply(combine_xt_xt1, tsxt2, tsxt12, SIMPLIFY = FALSE) 
  
  
  group <- list()
  group2 <- list()
  s <- list()
  s2 <- list()
  indx <- initi_label(ts,k)
  
  for (j in 1:k) {
    # in the original paper, they initially average the number l <- n/k, somehow they may cannot
    # be divided evenly. Here we just consider to devide them nearly equally
    group[[j]] <- sigma[indx==j]
    s[[j]] <- comaxe(group[[j]])
    group2[[j]] <- sigma2[indx == j]
    s2[[j]] <- comaxe(group2[[j]])
  }
  
  # the refinement step
  t <- 0
  tmax <- 100
  toterror <- rep(0,tmax+1) # this is, the total projection error
  error <-  matrix(rep(0,n*k), ncol = k) # this is, the projection error for each observation to each cluster
  
  while (t <= tmax) {
    t <- t+1 
    # indx <- NULL
    for (m in 1:n) { # for each observation
     
      for (l in 1:k) { # for each cluster
        y1 <- combined_list[[m]]%*%s[[l]][[1]] %*% t(s[[l]][[1]]) # the projection 
        y2 <- combined_list2[[m]]%*%s2[[l]][[1]] %*% t(s2[[l]][[1]]) 

        df <- normF(combined_list[[m]] - y1) + normF(combined_list2[[m]] - y2)
       error[m, l] <- df
      }
    }
    
    # a contains the minimum error  
    a <- apply(error, 1, function(x) {minerror = min(x)})
    # b contains the the group index for which it has the minim error
    indx <- apply(error, 1, function(x) {cluster_index = which.min(x) })
    
    toterror[t+1] <- sum(a)
    
    if (abs(toterror[t+1]-toterror[t])<.Machine$double.eps ) {
      break
    }
    
    # this helps to prevent the super bad initialization where we may
    # encounter a group without any MTS object, thus meet an error later
    while(length(unique(indx))!=k) {
      indx <- sample.int(k,n,replace = TRUE,prob = rep(1/k,k))
    }
   for (j in 1:k) {
    # in the original paper, they initially average the number l <- n/k, somehow they may cannot
    # be divided evenly. Here we just consider to devide them nearly equally
    group[[j]] <- sigma[indx==j]
    s[[j]] <- comaxe(group[[j]])
    group2[[j]] <- sigma2[indx == j]
    s2[[j]] <- comaxe(group2[[j]])
  }
     
    
  } 

  return(indx)

}


