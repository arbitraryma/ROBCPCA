library(base)


# This function is  Algorithm 1 in the paper

# This function can help to select p (number of components to retain) if it is not 
# specified by the user, it automatically retain p that account for more than 90% of variability 
# and we put a cap at 10 according to Hubert

# The output is the common projection axes and p

# Default is null input of p, one can also test the result when specifying the number of p,
# with a specified p, this function is then same as Li
comaxe <- function(sigma,p=NULL){
  # sigma <- group[[1]]
  sigma_com <- 0
  n <- length(sigma)
  for (i in 1:n) {
    sigma_com <- sigma_com + sigma[[i]]
  }
  sigma_com <- 1/n * sigma_com
  a <- svd(sigma_com)
  b <- a$d
  c <- a$u
  prop <- b/sum(b)
  cum <- cumsum(prop)
  # Calculate the cumulative variance explained
  cumulative_variance <- cum
  # Determine how many eigenvalues are needed to exceed 90% of the total variance
  if (is.null(p)) {
   p <-  which(cumulative_variance >= 0.9)[1]
  }else
  {p<-p }
  S <- c[,1:p]
  # common projection 
  com <- list(S,p)
  return(com)
}

