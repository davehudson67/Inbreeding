## load library
library(nimble)

#####################################################
######              CHECK R CODE               ######
#####################################################

## load random sampler
source("../Distributions/Exponential/Dist_Expo.R")

## Set parameters
r <- 0.3

## test non-truncated samples
sims <- rExpo(n=10000, r=r)
hist(sims, freq = FALSE)
## add PDF
lines(0:100, dExpo(0:100, r = r), col = "blue", lwd = 2)

## test left-truncated samples
sims <- rExpo(n=10000, r=r, zL = 20)
hist(sims, freq = FALSE)
## check against brute-force approach
sims1 <- NULL
while(length(sims1) < 10000) {
  temp <- rExpo(10000, r = r)
  temp <- temp[temp >= 20]
  sims1 <- c(sims1, temp)
}
sims1 <- sims1[1:10000]
hist(sims1, add = TRUE, col = "red", freq = FALSE)
## add PDF
lines(20:100, dExpo(20:100, r = r, zL = 20), col = "blue", lwd = 2)

## test right-truncated samples
sims <- rExpo(n=10000, r=r, zU = 5)
hist(sims, freq = FALSE)
## check against brute-force approach
sims1 <- NULL
while(length(sims1) < 10000) {
  temp <- rExpo(10000, r = r)
  temp <- temp[temp <= 5]
  sims1 <- c(sims1, temp)
}
sims1 <- sims1[1:10000]
hist(sims1, add = TRUE, col = "red", freq = FALSE)
## add PDF
lines(0:5, dExpo(0:5, r = r, zU = 5), col = "blue", lwd = 2)

## test interval-truncated samples
sims <- rExpo(n=10000, r=r, zL = 3, zU = 5)
hist(sims, freq = FALSE)
## check against brute-force approach
sims1 <- NULL
while(length(sims1) < 10000) {
  temp <- rExpo(10000, r = r)
  temp <- temp[temp >= 3 & temp <= 5]
  sims1 <- c(sims1, temp)
}
sims1 <- sims1[1:10000]
hist(sims1, add = TRUE, col = "red", freq = FALSE)
## add PDF
lines(3:5, dExpo(3:5, r = r, zL = 3, zU = 5), col = "blue", lwd = 2)

## check user input NA in R code
sims <- rExpo(n=10000, r=r, zL = 3, zU = NA)
summary(sims)

## check Survival times =======================================================================================================
age <- 1:80

## test non-truncated samples
sims <- rExpo(n=10000, r=r)
survivors <- NULL
for(i in 1:length(age)){
  survivors[i]<-(sum(sims>age[i])/10000)
}
plot(survivors, xlab = "Age", ylab = "Proportion of survivors")
lines(1:80, pExpo(1:80, r = r, lower.tail = 0), col = "blue")

## test left-truncated samples
sims <- rExpo(n=10000, r=r, zL = 20)
survivors <- NULL
for(i in 1:length(age)){
  survivors[i]<-(sum(sims>age[i])/10000)
}
plot(survivors, xlab = "Age", ylab = "Proportion of survivors")
lines(20:80, pExpo(20:80, r = r, zL = 20, lower.tail = 0), col = "blue")

## test right-truncated samples
sims <- rExpo(n=10000, r=r, zU = 20)
survivors <- NULL
for(i in 1:length(age)){
  survivors[i]<-(sum(sims>age[i])/10000)
}
plot(survivors, xlab = "Age", ylab = "Proportion of survivors")
lines(1:20, pExpo(1:20, r = r, zU = 20, lower.tail = 0), col = "blue")

## test interval-truncated samples
sims <- rExpo(n=10000, r=r, zL = 5, zU = 20)
survivors <- NULL
for(i in 1:length(age)){
  survivors[i]<-(sum(sims>age[i])/10000)
}
plot(survivors, xlab = "Age", ylab = "Proportion of survivors")
lines(5:20, pExpo(5:20, r = r, zL = 5, zU = 20, lower.tail = 0), col = "blue")

#####################################################
######          CHECK MODELS COMPILE           ######
#####################################################


##  can use exp distribution in Nimble.


#####################################################
######         CHECK USER INPUT ERRORS        #######
#####################################################

# Different dimensions
r<-c(r,r)
x<- 12
zL <- NA
zU <- NA
dExpo(x,r,zL,zU)
x<- c(4, 13)
dExpo(x,r,zL,zU)
pExpo(x,r,zL,zU)

#zL greater than zU
r <- exp(-2)
zL<- 10
zU<- 3
dExpo(x,r,zL,zU)
rExpo(1,r,zL,zU)
pExpo(0.4,r,zL,zU)
qExpo(5,r,zL,zU)

#zL greater than x
zU<-NA
x <- 4
dExpo(x,r,zL,zU)
pExpo(x,r,zL,zU)

#zU less than x
zU <- 3
zL <- NA
dExpo(x,r,zL,zU)
pExpo(x,r,zL,zU)

#Negative parameter
r<--r
dExpo(2,r,zL,zU)
dExpo(c(1, 2),r,zL,zU)
rExpo(4,r,zL,zU)

## length argument to rExpo
r <- -r
rExpo(4,r,zL,zU)
rExpo(c(2, 4),r,zL,zU)
