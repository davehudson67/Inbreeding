## load library
library(nimble)

#####################################################
######              CHECK R CODE               ######
#####################################################

## load random sampler
source("../Distributions/Gompertz/Dist_Gompertz.R")
source("../Distributions/Gompertz/Dist_GompertzNim.R")

## Set parameters
a <- 0.05
b <- 0.02

## test non-truncated samples
sims <- rGompertz(n = 10000, a = a, b = b)
hist(sims, freq = FALSE)
## add PDF
lines(0:100, dGompertz(0:100, a = a, b = b), col = "blue", lwd = 2)

## test left-truncated samples
sims <- rGompertz(n = 10000, a = a, b = b, zL = 20)
hist(sims, freq = FALSE)
## check against brute-force approach
sims1 <- NULL
while(length(sims1) < 10000) {
  temp <- rGompertz(10000, a = a, b = b)
  temp <- temp[temp >= 20]
  sims1 <- c(sims1, temp)
}
sims1 <- sims1[1:10000]
hist(sims1, add = TRUE, col = "red", freq = FALSE)
## add PDF
lines(20:100, dGompertz(20:100, a, b, zL = 20), col = "blue", lwd = 2)

## test right-truncated samples
sims <- rGompertz(n = 10000, a = a, b = b, zU = 5)
hist(sims, freq = FALSE)
## check against brute-force approach
sims1 <- NULL
while(length(sims1) < 10000) {
  temp <- rGompertz(10000, a = a, b = b)
  temp <- temp[temp <= 5]
  sims1 <- c(sims1, temp)
}
sims1 <- sims1[1:10000]
hist(sims1, add = TRUE, col = "red", freq = FALSE)
## add PDF
lines(0:5, dGompertz(0:5, a, b, zU = 5), col = "blue", lwd = 2)

## test interval-truncated samples
sims <-  rGompertz(n = 10000, a = a, b = b, zL = 3, zU = 5)
hist(sims, freq = FALSE)
## check against brute-force approach
sims1 <- NULL
while(length(sims1) < 10000) {
  temp <- rGompertz(10000, a = a, b = b)
  temp <- temp[temp >= 3 & temp <= 5]
  sims1 <- c(sims1, temp)
}
sims1 <- sims1[1:10000]
hist(sims1, add = TRUE, col = "red", freq = FALSE)
 ## add PDF
lines(seq(3, 5, length = 100), dGompertz(seq(3, 5, length = 100), a = a, b = b, zL = 3, zU = 5), col = "blue", lwd = 2)

## check user input NA in R code
sims <- sapply(1:10000, function(x) rGompertz(n = 1, a = a, b = b, zL = 3, zU = NA))
summary(sims)

## check Survival times =======================================================================================================
age <- 1:80

## test non-truncated samples
sims <- rGompertz(n = 10000, a = a, b = b)
survivors <- NULL
for(i in 1:length(age)){
  survivors[i]<-(sum(sims>age[i])/10000)
}
plot(survivors, xlab = "Age", ylab = "Proportion of survivors")
lines(1:80, pGompertz(1:80, a = a, b = b, lower.tail = 0), col = "blue")

## test left-truncated samples
sims <- rGompertz(n = 10000, a = a, b = b, zL = 20)
survivors <- NULL
for(i in 1:length(age)){
  survivors[i]<-(sum(sims>age[i])/10000)
}
plot(survivors, xlab = "Age", ylab = "Proportion of survivors")
lines(20:80, pGompertz(20:80, a = a, b = b, zL = 20, lower.tail = 0), col = "blue")

## test right-truncated samples
sims <- rGompertz(n = 10000, a = a, b = b, zU = 20)
survivors <- NULL
for(i in 1:length(age)){
  survivors[i]<-(sum(sims>age[i])/10000)
}
plot(survivors, xlab = "Age", ylab = "Proportion of survivors")
lines(1:20, pGompertz(1:20, a = a, b = b, zU = 20, lower.tail = 0), col = "blue")

## test interval-truncated samples
sims <- rGompertz(n = 10000, a = a, b = b, zL = 5, zU = 20)
survivors <- NULL
for(i in 1:length(age)){
  survivors[i]<-(sum(sims>age[i])/10000)
}
plot(survivors, xlab = "Age", ylab = "Proportion of survivors")
lines(5:20, pGompertz(5:20, a = a, b = b, zL = 5, zU = 20, lower.tail = 0), col = "blue")

#####################################################
######          CHECK MODELS COMPILE           ######
#####################################################

## check components compile
compGompertzNim <- compileNimble(dgompzNim)
compGompertzNim <- compileNimble(pgompzNim)
compGompertzNim <- compileNimble(qgompzNim)
compGompertzNim <- compileNimble(rgompzNim)

#####################################################
######          CHECK COMPILED MODEL           ######
#####################################################

## check compiled model
compGompertzNim <- compileNimble(rgompzNim)
 
## test non-truncated samples
sims <- sapply(1:10000, function(x) compGompertzNim(n = 1, a = a, b = b))
hist(sims, freq = FALSE)
## check against brute-force approach
lines(seq(0, 100, length = 100), dgompzNim(seq(0, 100, length = 100), a = a, b = b), col = "blue", lwd = 2)

## test survivor function
survivors <- NULL
for(i in 1:length(age)){
  survivors[i]<-(sum(sims>age[i])/10000)
}
plot(survivors, xlab = "Age", ylab = "Proportion of survivors")
lines(1:80, pgompzNim(1:80, a = a, b = b, lower.tail = 0), col = "blue")

#####################################################
######         CHECK USER INPUT ERRORS        #######
#####################################################

# Different dimensions
a<-c(a,a)
x<- 12
zL <- NA
zU <- NA
dGompertz(x,a,b,zL,zU)
x<- c(4, 13)
dGompertz(x,a,b,zL,zU)
pGompertz(x,a,b,zL,zU)

#zL greater than zU
a <- exp(-2)
zL<- 10
zU<- 3
dGompertz(x,a,b,zL,zU)
rGompertz(1,a,b,zL,zU)
pGompertz(0.4,a,b,zL,zU)
qGompertz(5,a,b,zL,zU)

#zL greater than x
zU<-NA
x <- 4
dGompertz(x,a,b,zL,zU)
pGompertz(x,a,b,zL,zU)

#zU less than x
zU <- 3
zL <- NA
dGompertz(x,a,b,zL,zU)
pGompertz(x,a,b,zL,zU)

#Negative parameter
a<--a
dGompertz(2,a,b,zL,zU)
dGompertz(c(1, 2),a,b,zL,zU)
rGompertz(4,a,b,zL,zU)

## length argument to rGompertz
a <- -a
rGompertz(4,a,b,zL,zU)
rGompertz(c(2, 4),a,b,zL,zU)
