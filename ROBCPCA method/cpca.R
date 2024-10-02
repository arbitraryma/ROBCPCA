library(base)
library(stats)
library(fungible)

cpca <- function(ts,k,p=NULL)
{
  
  # the initial step 
  n <- length(ts) # number of MTS  observations
  k <- k   # number of clusters 
  l <- c(rep(1/k,k))  # initial probobility in each cluster to have the common projection axes, evenly distributed
  dim <- ncol(ts[[1]])  # for cpca, we assume all ts are with equal dimension
  ts1 <- list()
  sigma <- list()
 
  for (i in 1:n) {
    
    # robustly estimate the location and normalization 
    ts1[[i]] <- sweep(ts[[i]], 2, colMeans(ts[[i]]), `-`)
    sigma[[i]] <- cov(ts1[[i]])
  }
  
  group <- list()
  s <- list()
  # repidx <- matrix(rep(0,n*3),ncol = 3)
  # miner <- rep(0,3)
  # for(ran in 1:3 ){ # the result can be significantly influence by the initiation,
  # we do three times of initiation and choose the best one to be the final clustering result
  # initialization: random allocate data into different clusters 
  indx <- sample.int(k,n,replace = TRUE,prob = rep(1/k,k))
  
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
        y <- ts1[[m]]%*%s[[l]][[1]] %*% t(s[[l]][[1]]) # the projection 
        df <- ts1[[m]] - y
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
  #   repidx[,ran] <- indx
  #   miner[ran] <- toterror[t+1] 
  # }
  # r is the iterate process time, break until the error converges 
  #fin <- which.min(miner)
  return(indx)
  # return(repidx[,fin])
}
