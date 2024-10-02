##### Loading needed packages
library(signal)
#install.packages("remotes")   # if you don't have the package already
#remotes::install_github("nickpoison/astsa/astsa_build")
library(astsa)
library(plot.matrix)

##### Extracting oscillations from a signal
extract_rhythms <- function(x,fs,rhythms = c("delta","theta","alpha","beta","gamma"),
                            method = "butter",twosided = FALSE){
  y <- vector("list",length(rhythms))
  for(i in 1:length(rhythms)){
    if(rhythms[i] == "delta"){
      if(method == "butter"){
        filter_coeffs <- signal::butter(4, c(0.5,4)/(0.5 * fs), type="pass")
      } else if(method == "fir"){
        filter_coeffs <- signal::fir1(100, c(0.5,4)/(0.5 * fs), type="pass")
      }
      if(twosided){
        y[[i]] <- as.numeric(signal::filtfilt(filter_coeffs, x))
      } else {
        y[[i]] <- as.numeric(signal::filter(filter_coeffs, x))
      }
    } else if(rhythms[i] == "theta"){
      if(method == "butter"){
        filter_coeffs <- signal::butter(4, c(4,8)/(0.5 * fs), type="pass")
      } else if(method == "fir"){
        filter_coeffs <- signal::fir1(100, c(4,8)/(0.5 * fs), type="pass")
      }
      if(twosided){
        y[[i]] <- as.numeric(signal::filtfilt(filter_coeffs, x))
      } else {
        y[[i]] <- as.numeric(signal::filter(filter_coeffs, x))
      }
    } else if(rhythms[i] == "alpha"){
      if(method == "butter"){
        filter_coeffs <- signal::butter(4, c(8,12)/(0.5 * fs), type="pass")
      } else if(method == "fir"){
        filter_coeffs <- signal::fir1(100, c(8,12)/(0.5 * fs), type="pass")
      }
      if(twosided){
        y[[i]] <- as.numeric(signal::filtfilt(filter_coeffs, x))
      } else {
        y[[i]] <- as.numeric(signal::filter(filter_coeffs, x))
      }
    } else if(rhythms[i] == "beta"){
      if(method == "butter"){
        filter_coeffs <- signal::butter(4, c(12,30)/(0.5 * fs), type="pass")
      } else if(method == "fir"){
        filter_coeffs <- signal::fir1(100, c(12,30)/(0.5 * fs), type="pass")
      }
      if(twosided){
        y[[i]] <- as.numeric(signal::filtfilt(filter_coeffs, x))
      } else {
        y[[i]] <- as.numeric(signal::filter(filter_coeffs, x))
      }
    } else if(rhythms[i] == "gamma"){
      if(method == "butter"){
        filter_coeffs <- signal::butter(4, c(30,45)/(0.5 * fs), type="pass")
      } else if(method == "fir"){
        filter_coeffs <- signal::fir1(100, c(30,45)/(0.5 * fs), type="pass")
      }
      if(twosided){
        y[[i]] <- as.numeric(signal::filtfilt(filter_coeffs, x))
      } else {
        y[[i]] <- as.numeric(signal::filter(filter_coeffs, x))
      }
    }
    
  }
  return(y)
}

##### Generate AR(2) coefficients for specific frequency band
fband.arcoefs <- function(fband,samp.rate = 100){
  if(fband == "delta"){
    peaklocation <- 2
    startband <- 0.5
    endband <- 4
    sharpness <- 0.03
  } else if(fband == "theta"){
    peaklocation <- 6
    startband <- 4
    endband <- 8
    sharpness <- 0.03
  } else if(fband == "alpha"){
    peaklocation <- 10
    startband <- 8
    endband <- 12
    sharpness <- 0.03
  } else if(fband == "beta"){
    peaklocation <- 22.5
    startband <- 12
    endband <- 30
    sharpness <- 0.05
  } else if(fband == "gamma"){
    peaklocation <- 37.5
    startband <- 30
    endband <- 45
    sharpness <- 0.075
  }
  psi<- peaklocation/samp.rate
  par <- c(exp(-sharpness)*2*cos(2*pi*psi), -exp(-2*sharpness))
  return(par)
}



# generate weights for contribution from the latent processes 
gen.coef <- function(idx, p) {   
  #idx is a vector created with c(), contains the index of the bands that contribute the most
  #1 - delta, 2 - theta, 3 - alpha, 4 - beta, 5 - gamma
  #p is the sum of the coefficients corresponding to the bands on idx
  
  coeffs <- numeric(5)
  w1 <- p/length(idx)
  eps1 <- 1e-1
  
  
  for (i in idx) {
    coeffs[i] <- ifelse(i != max(idx), runif(1, w1 - eps1, w1 + eps1), 
                        p - sum(coeffs[idx[-length(idx)]]))
  }
  if (p < 1) {
    idxc <- (1:5)[-idx]
    w2 <- (1-p)/(5 - length(idx))
    eps2 <- 1e-2
    for (i in idxc) {
      coeffs[i] <- ifelse(i != max(idxc), runif(1, w2 - eps2, w2 + eps2), 
                          1 - p - sum(coeffs[idxc]))
    }
  }
  
  coeffs   
}

#generates the latent process from a specific band 
#bands are delta, theta, alpha, beta and gamma
gen.latent <- function(fband, lenT, samp.rate = 100) {
  #fband is the frequency band
  
  #lenT is the length of the time series
  
  #samp.rate is the sampling rate
  
  phis <- fband.arcoefs(fband, samp.rate)
  Z <- arima.sim(lenT, model = list(ar = c(phis[1], phis[2])))
  Z <- Z - mean(Z)
  Z <- extract_rhythms(Z,fs = samp.rate,rhythms = fband,
                       method = "butter",twosided = TRUE)[[1]]
  Z
}

#generates multivariate time series as mixing of latent processes
gen.ts <- function(d, p, channels, clusters, lenT, samp.rate = 100) {
  #d is the dimension of the resulting time series
  
  #p is a c() vector with probabilities as in gen.coef, one p per cluster
  
  #channels is a list with c() vectors, the i-th c() vector of the list must
  #contain the indices of which channels will contribute the most for the i-th cluster
  
  #list of c() vectors, the c() vector on the i-th position 
  #contains which components of the time series belong to cluster i 
  
  #lenT is the length of the time series
  
  #samp.rate is the sampling rate
  
  A <- t(gen.A(d, p, channels, clusters))
  bands <- c("delta", "theta", "alpha", "beta", "gamma")
  Z <- apply(do.call(rbind, lapply(bands, function(fband) gen.latent(fband, lenT))), 
             1, scale)
  Z %*% A
}


### Generate AR(2) coefficients manually for delta, theta, alpha, beta and gamma processes
### INPUT: (1) peaklocation, (2) sampling rate, (3) sharpness/bandwidth 
### OUTPUT: AR(2) coefficients phi1 and phi2

T = 100;
timeT = (1:T);
samprate = 100 #Hertz

### Delta parameters 
peaklocation_delta = 2.00
sharpness_delta = 0.05 

M_delta = exp(sharpness_delta);
psi_delta = peaklocation_delta/samprate;
phi1_delta = (2/M_delta)*cos(2*pi*psi_delta);
phi2_delta = -(1/(M_delta^2));

### Theta Oscillations 
### Inputs
peaklocation_theta = 6.00
sharpness_theta = 0.05  

M_theta = exp(sharpness_theta);
psiM_theta = peaklocation_theta/samprate;
phi1_theta = (2/M_theta)*cos(2*pi*psiM_theta);
phi2_theta = -(1/(M_theta^2));


### Alpha Oscillations 
### Inputs
peaklocation_alpha = 10.00
sharpness_alpha = 0.05 #range (0, \infty)

M_alpha = exp(sharpness_alpha);
psi_alpha = peaklocation_alpha/samprate;
phi1_alpha = (2/M_alpha)*cos(2*pi*psi_alpha);
phi2_alpha = -(1/(M_alpha^2));


### Beta Oscillations 
### Inputs
peaklocation_beta = 20.00
sharpness_beta = 0.10 #range (0, \infty)

M_beta = exp(sharpness_beta);
psi_beta = peaklocation_beta/samprate;
phi1_beta = (2/M_beta)*cos(2*pi*psi_beta);
phi2_beta = -(1/(M_beta^2));



### Gamma Oscillations 
### Inputs
peaklocation_gamma = 40.00
sharpness_gamma = 0.10 #range (0, \infty)

M_gamma = exp(sharpness_gamma);
psi_gamma = peaklocation_gamma/samprate;
phi1_gamma = (2/M_gamma)*cos(2*pi*psi_gamma);
phi2_gamma = -(1/(M_gamma^2));



bands <- c("delta", "theta", "alpha", "beta", "gamma")

# for each replications, we generate different latent process and different coefficient matrices
# the matrices are generated randomly using the uniform distributions, the dominated frequncies will be given 
# the most proportions and others will be given a very small amount of percentages 

# for cluster 1, we generate processes that are dominated by delta and theta 





# here, we consider the matrix to have length 100, with 32/64/ 256 channels 

n <- 100# this is, the time length, for EEG, it is generally considered to be the same 

# generate the latent process

# construct the coefficient matrix to obtain the final MTS, using function gen.coef



set.seed(1)
f <- c(32,64,128)  # we consider the common choice of # channels 

arob <- matrix(rep(0,length(f)*rep),nrow = length(f))
ac <- matrix(rep(0,length(f)*rep),nrow = length(f))
arob1   <- matrix(rep(0,length(f)*rep),nrow = length(f))
ac1 <- matrix(rep(0,length(f)*rep),nrow = length(f))
autoselectac  <- matrix(rep(0,length(f)*rep),nrow = length(f))
auto_iniac  <- matrix(rep(0,length(f)*rep),nrow = length(f)) 
auto_ini_l1ac  <- matrix(rep(0,length(f)*rep),nrow = length(f)) 
roblagac1 <- matrix(rep(0,length(f)*rep),nrow = length(f)) 
roblagac2   <- matrix(rep(0,length(f)*rep),nrow = length(f)) 
roblagac12 <- matrix(rep(0,length(f)*rep),nrow = length(f)) 

truearobinx <- matrix(rep(0,length(f)*rep),nrow = length(f))
trueacinx <- matrix(rep(0,length(f)*rep),nrow = length(f))
truearobinx1 <- matrix(rep(0,length(f)*rep),nrow = length(f))
trueacinx1 <- matrix(rep(0,length(f)*rep),nrow = length(f))
trueauto <- matrix(rep(0,length(f)*rep),nrow = length(f))
trueauto_ini <- matrix(rep(0,length(f)*rep),nrow = length(f))
trueroglag1  <- matrix(rep(0,length(f)*rep),nrow = length(f)) 
trueroglag2  <- matrix(rep(0,length(f)*rep),nrow = length(f)) 
trueroglag12  <- matrix(rep(0,length(f)*rep),nrow = length(f)) 

labels <- c(rep(1,10),rep(2,10),rep(3,10))
# the MTS for cluster one 
# for each cluster, we only generate 5 MTS 

for (j in 1:length(dimen)) {
  f <- f[j]
  # control the same coefficient for each replication
  coef1 <- sapply(1:dim, function(x) { gen.coef(c(1,2,3), 0.95) })
  coef2 <- sapply(1:dim, function(x) {gen.coef(c(2,4), 0.95)})
  coef3 <- sapply(1:dim, function(x) {gen.coef(c(4,5), 0.95)})
  
  for (l in 1:100) {
 # check result for roblag1 after l=57
    group1_data_list <- list()
    group2_data_list <- list()
    group3_data_list <- list()
# Perform the matrix multiplication and store in the list

  
# for cluster 1, we consider the main mixture are delta and theta
# the proportion of the dominated is 95%,  with 5% comes from other frequency band 
for (i in 1:10) {
 # group1_data_list[[i]] <- apply(do.call(rbind, lapply(bands, function(fband) gen.latent(fband, n))), 
 #                                1, scale) %*% sapply(1:dim, function(x) {
 #                                  gen.coef(c(1,2,3), 0.95)
  #                               })
  group1_data_list[[i]] <- apply(do.call(rbind, lapply(bands, function(fband) gen.latent(fband, n))), 
                                 1, scale)  %*% coef1
}


# for cluster 2, we consider the main mixture are   alpha and beta 

for (i in 1:10) {
 # group2_data_list[[i]] <- apply(do.call(rbind, lapply(bands, function(fband) gen.latent(fband, n))), 
 #                                1, scale) %*% sapply(1:dim, function(x) {
 #                                  gen.coef(c(2,3), 0.95)
 #                                })
  group2_data_list[[i]] <-  apply(do.call(rbind, lapply(bands, function(fband) gen.latent(fband, n))), 
                                  1, scale)  %*% coef2
}


for (i in 1:10) {
 # group3_data_list[[i]] <- apply(do.call(rbind, lapply(bands, function(fband) gen.latent(fband, n))), 
 #                                1, scale) %*% sapply(1:dim, function(x) {
 #                                  gen.coef(c(3,4,5), 0.95)
 #                               })
  
  group3_data_list[[i]] <-  apply(do.call(rbind, lapply(bands, function(fband) gen.latent(fband, n))), 
                                  1, scale)  %*% coef3
}
  
 ts <- c(group1_data_list,group2_data_list,group3_data_list)

 robinx <- cpca(ts,3,p=1)
 cinx <- cpca(ts,3,p=2)
 robinx1 <- cpca(ts,3,p=3)
 cinx1 <- cpca(ts,3,p=4)
 auto <- cpca(ts,3)
 auto_ini <- robcpca(ts,3)
 roblag1 <- robcpca_lag1(ts,3)
 roblag2 <- robcpca_lag_2(ts,3)
 roblag12 <- robcpca_lag_12(ts,3)
 auto_ini_l1 <- cpca_initial_man(ts,2)
 
 arob[j,l] <- adjustedRandIndex(robinx,labels)
 ac[j,l] <- adjustedRandIndex(cinx,labels)
 arob1[j,l] <- adjustedRandIndex(robinx1,labels)
 ac1[j,l] <- adjustedRandIndex(cinx1,labels)
 autoselectac[j,l] <- adjustedRandIndex(auto,labels)
auto_iniac[j,l] <-  adjustedRandIndex(auto_ini,labels) 
 roblagac1[j,l] <-  adjustedRandIndex(roblag1,labels)
 roblagac2[j,l] <-  adjustedRandIndex(roblag2,labels)
  roblagac12[j,l] <-  adjustedRandIndex(roblag12,labels)

 truearobinx[j,l] <-  true_accuracy(labels, robinx)
 
 trueacinx[j,l] <-  true_accuracy(labels, cinx)
 
 truearobinx1[j,l] <-  true_accuracy(labels, robinx1)
 
 trueacinx1[j,l] <- true_accuracy(labels, cinx1)
 
 trueauto[j,l] <-  true_accuracy(labels,auto)
 
 trueauto_ini[j,l] <-  true_accuracy(labels,auto_ini)
 
  trueroglag1[j,l] <-  true_accuracy(labels,roblag1)
 
  trueroglag2[j,l] <-  true_accuracy(labels,roblag2)
 
 trueroglag12[j,l] <-  true_accuracy(labels,roblag12)


}

}



# the ARI (adjusted random index)
round(rowMeans(arob[1:3,1:29]),4)
round(rowMeans(ac[1:3,1:100]),4)
round(rowMeans(arob1[1:3,1:100]),4)
round(rowMeans(ac1[1:3,1:100]),4)
round(rowMeans(autoselectac[1:3,1:50]),4)
round(rowMeans(auto_iniac[1:3,1:100]),4) 
round(rowMeans(roblagac1[1:3,1:100]),4)
round(rowMeans(roblagac2[1:3,1:100]),4)
round(rowMeans(roblagac12[1:3,1:100]),4)

# 0.8291 0.5999 0.7198
# 0.8068 0.6065 0.7366
# 0.7951 0.5746 0.7042
# 0.8081 0.5901 0.7356
# 0.8081 0.5901 0.7356
# 0.8948 0.7756 0.8787
# 1 1 1
# 1 1 1
# 1 1 1



# sd 
round(apply(arob[1:3,1:100], 1, sd),4)
round(apply(ac[1:3,1:100], 1, sd),4)
round(apply(arob1[1:3,1:100] , 1, sd),4)
round(apply(ac1[1:3,1:100] , 1, sd),4)
round(apply(autoselectac[1:3,1:100], 1, sd),4)
round(apply(auto_iniac[1:3,1:100], 1, sd),4) 
round(apply(roblagac1[1:3,1:100], 1, sd),4)
round(apply(roblagac2[1:3,1:100], 1, sd),4)
round(apply(roblagac12[1:3,1:100], 1, sd),4)
# 0.1715 0.1997 0.2400
# 0.1760 0.2157 0.2448
# 0.1908 0.1907 0.2496
# 0.2038 0.2020 0.2617
# 0.2081 0.2101 0.2403
# 0.0830 0.0737 0.1164
# 0 0 0
# 0 0 0
# 0 0 0


# True accuracy 

round(rowMeans(truearobinx[1:3,1:100]),4)  
round(rowMeans(trueacinx[1:3,1:100]),4)
round(rowMeans(truearobinx1[1:3,1:100]),4)
round(rowMeans(trueacinx1[1:3,1:100]),4)
round(rowMeans(trueauto[1:3,1:100]),4)
round(rowMeans(trueauto_ini[1:3,1:100]),4)
round(rowMeans(trueroglag1[1:3,1:100]),4)
round(rowMeans(trueroglag2[1:3,1:100]),4)
round(rowMeans(trueroglag12[1:3,1:100]),4)
#  0.9147 0.7860 0.8370
#  0.9143 0.7927 0.8113
#  0.8953 0.7970 0.8277
#  0.8830 0.7637 0.7990
#  0.8893 0.7993 0.8240
#  0.9640 0.9173 0.9473 
#  1 1 1
#  1 1 1
#  1 1 1

# sd true accuracy
round(apply(truearobinx[1:3,1:100], 1, sd),4)
round(apply(trueacinx[1:3,1:100] , 1, sd),4)
round(apply(truearobinx1[1:3,1:100] , 1, sd),4)
round(apply(trueacinx1[1:3,1:100] , 1, sd),4)
round(apply(trueauto[1:3,1:100] , 1, sd),4)
round(apply(trueauto_ini[1:3,1:100] , 1, sd),4)  
round(apply(trueroglag1[1:3,1:100], 1, sd),4)
round(apply(trueroglag2[1:3,1:100], 1, sd),4)
round(apply(trueroglag12[1:3,1:100], 1, sd),4)

#  0.1266 0.1600 0.1899
#  0.1282 0.1677 0.1980
#  0.1425 0.1500 0.1948
#  0.1574 0.1706 0.2053
#  0.1584 0.1610 0.1942
#  0.0298 0.0356 0.0833
#  0 0 0
#  0 0 0
#  0 0 0 











# for this simulation, we consider some noises generated by t-distribution with df=3

dimen <- c(32,64,128)  # we consider the common choice of # channels 
label<- c(rep(1,10),rep(2,10),rep(3,10))
accuracy <- matrix(rep(0,3*100),nrow = 3)
accuracyrob <- matrix(rep(0,3*100),nrow = 3)
trueac <- matrix(rep(0,3*100),nrow = 3)
trueacrob <-  matrix(rep(0,3*100),nrow = 3)

# the MTS for cluster one 
# for each cluster, we only generate 5 MTS 
for (j in 1:length(dimen)) {
  f <- f[j]
  # control the same coefficient for each replication
  coef1 <- sapply(1:dim, function(x) { gen.coef(c(1,2,3), 0.95) })
  coef2 <- sapply(1:dim, function(x) {gen.coef(c(2,4), 0.95)})
  coef3 <- sapply(1:dim, function(x) {gen.coef(c(4,5), 0.95)})
  for (l in 1:50) {
    
    group1_data_list <- list()
    group2_data_list <- list()
    group3_data_list <- list()
    # Perform the matrix multiplication 5 times and store in the list
    
    
    # Perform the matrix multiplication and store in the list
    
 
    for (i in 1:10) {
  
      group1_data_list[[i]] <- apply(do.call(rbind, lapply(bands, function(fband) gen.latent(fband, n))), 
                                     1, scale)  %*% coef1
    }
    
    
    for (i in 1:10) {

      group2_data_list[[i]] <-  apply(do.call(rbind, lapply(bands, function(fband) gen.latent(fband, n))), 
                                      1, scale)  %*% coef2
    }
    
    
    for (i in 1:10) {
      group3_data_list[[i]] <-  apply(do.call(rbind, lapply(bands, function(fband) gen.latent(fband, n))), 
                                      1, scale)  %*% coef3
    }
    
    group1_data_list <- add_noise_to_group(group1_data_list)
    group2_data_list <- add_noise_to_group(group2_data_list)
    group3_data_list <- add_noise_to_group(group3_data_list)
    
    
    
    ts <- c(group1_data_list,group2_data_list,group3_data_list)
    
    robinx <- cpca(ts,3,p=1)
    cinx <- cpca(ts,3,p=2)
    robinx1 <- cpca(ts,3,p=3)
    cinx1 <- cpca(ts,3,p=4)
    auto <- cpca(ts,3)
    auto_ini <- robcpca(ts,3)
    roblag1 <- robcpca_lag1(ts,3)
    roblag2 <- robcpca_lag_2(ts,3)
    roblag12 <- robcpca_lag_12(ts,3)
    # auto_ini_l1 <- cpca_initial_man(ts,2)
    
    arob[j,l] <- adjustedRandIndex(robinx,labels)
    ac[j,l] <- adjustedRandIndex(cinx,labels)
    arob1[j,l] <- adjustedRandIndex(robinx1,labels)
    ac1[j,l] <- adjustedRandIndex(cinx1,labels)
    autoselectac[j,l] <- adjustedRandIndex(auto,labels)
    auto_iniac[j,l] <-  adjustedRandIndex(auto_ini,labels) 
    roblagac1[j,l] <-  adjustedRandIndex(roblag1,labels)
    roblagac2[j,l] <-  adjustedRandIndex(roblag2,labels)
    roblagac12[j,l] <-  adjustedRandIndex(roblag12,labels)
    
    truearobinx[j,l] <-  true_accuracy(labels, robinx)
    
    trueacinx[j,l] <-  true_accuracy(labels, cinx)
    
    truearobinx1[j,l] <-  true_accuracy(labels, robinx1)
    
    trueacinx1[j,l] <- true_accuracy(labels, cinx1)
    
    trueauto[j,l] <-  true_accuracy(labels,auto)
    
    trueauto_ini[j,l] <-  true_accuracy(labels,auto_ini)
    
    trueroglag1[j,l] <-  true_accuracy(labels,roblag1)
    
    trueroglag2[j,l] <-  true_accuracy(labels,roblag2)
    
    trueroglag12[j,l] <-  true_accuracy(labels,roblag12)
    
    
  }
  
}

