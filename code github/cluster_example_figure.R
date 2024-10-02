source('/Users/maz0b/Desktop/2nd semester/contamination_func.R')
library(graphics)


m1 <- rep(0,2)
s1 <- generateCorMat(2)

a1 <- dataclustercon1(100,m1,muc(s1,m1),s1,0.1)

m2 <- rep(0.5,2)
s2 <- generateCorMat(2,corrType = 'A09')

b2 <- dataclustercon1(100,m2,muc(s2,m2),s2,0.1)

r <- matrix(c(1,0.5,-0.5,1),ncol = 2)
a2 <- mvrnorm(100,rep(-5,2),r)

q <- matrix(c(1,-0.8,0.8,1),ncol = 2)
b1 <- mvrnorm(100,rep(-10,2),q)

plot(a1[[2]][,1],a1[[2]][,2],col='darkorchid1',pch = 2,ylim=c(-30,5),xlim=c(-25,5),
     ylab="",xlab = '')

points(b2[[2]][,1],b2[[2]][,2],col='deepskyblue3',pch = 3)


points(a2[,1],a2[,2],col='salmon2',pch = 5)
points(b1[,1],b1[,2],col='pink',pch = 6)
