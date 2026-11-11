## load library
library(nimble)
library(purrr)

#####################################################
######              CHECK R CODE               ######
#####################################################

## load random sampler
source("../Distributions/Siler/Dist_Siler.R")
## Set parameters
a1 <- exp(-2)
a2 <- exp(-11)
b1 <- 1.2
b2 <- 0.11
c <- exp(-5) 

## test non-truncated samples
sims <- rSiler(n = 10000, a1 = a1, a2 = a2, b1 = b1, b2 = b2, c = c)
hist(sims, freq = FALSE)
## add PDF
lines(0:100, dSiler(0:100, a1 = a1, a2 = a2, b1 = b1, b2 = b2, c = c), col = "blue", lwd = 2)

## test left-truncated samples
sims <- rSiler(n = 10000, a1 = a1, a2 = a2, b1 = b1, b2 = b2, c = c, zL = 20)
hist(sims, freq = FALSE)
## check against brute-force approach
sims1 <- NULL
while(length(sims1) < 10000) {
  temp <- rSiler(10000, a1, a2, b1, b2, c)
  temp <- temp[temp >= 20]
  sims1 <- c(sims1, temp)
}
sims1 <- sims1[1:10000]
hist(sims1, add = TRUE, col = "red", freq = FALSE)
## add PDF
lines(20:100, dSiler(20:100, a1, a2, b1, b2, c, zL = 20), col = "blue", lwd = 2)

## test right-truncated samples
sims <- rSiler(n = 10000, a1 = a1, a2 = a2, b1 = b1, b2 = b2, c = c, zU = 5)
hist(sims, freq = FALSE)
## check against brute-force approach
sims1 <- NULL
while(length(sims1) < 10000) {
  temp <- rSiler(10000, a1 = a1, a2 = a2, b1 = b1, b2 = b2, c = c)
  temp <- temp[temp <= 5]
  sims1 <- c(sims1, temp)
}
sims1 <- sims1[1:10000]
hist(sims1, add = TRUE, col = "red", freq = FALSE)
## add PDF
lines(0:5, dSiler(0:5, a1 = a1, a2 = a2, b1 = b1, b2 = b2, c = c, zU = 5), col = "blue", lwd = 2)

## test interval-truncated samples
sims <- rSiler(n = 10000, a1 = a1, a2 = a2, b1 = b1, b2 = b2, c = c, zL = 3, zU = 5)
hist(sims, freq = FALSE)
## check against brute-force approach
sims1 <- NULL
while(length(sims1) < 10000) {
  temp <- rSiler(10000, a1 = a1, a2 = a2, b1 = b1, b2 = b2, c = c)
  temp <- temp[temp >= 3 & temp <= 5]
  sims1 <- c(sims1, temp)
}
sims1 <- sims1[1:10000]
hist(sims1, add = TRUE, col = "red", freq = FALSE)
## add PDF
lines(seq(3, 5, length = 100), dSiler(seq(3, 5, length = 100), a1 = a1, a2 = a2, b1 = b1, b2 = b2, c = c, zL = 3, zU = 5), col = "blue", lwd = 2)

## check user input NA in R code
sims <- sapply(1:10000, function(x) rSiler(n = 1, a1 = a1, a2 = a2, b1 = b1, b2 = b2, c = c, zL = 3, zU = NA))
summary(sims)

## check speed
system.time(sapply(1:10000, function(x) rSiler(n = 1, a1 = a1, a2 = a2, b1 = b1, b2 = b2, c = c, zL = 3, zU = NA)))
system.time(rSiler(10000, a1 = a1, a2 = a2, b1 = b1, b2 = b2, c = c, zL = 3, zU = NA))

## check Survival times =======================================================================================================
age <- 1:80

## test non-truncated samples
sims <- rSiler(n = 10000, a1 = a1, a2 = a2, b1 = b1, b2 = b2, c = c)
survivors <- NULL
for(i in 1:length(age)){
  survivors[i]<-(sum(sims>age[i])/10000)
}
plot(survivors, xlab = "Age", ylab = "Proportion of survivors")
lines(1:80, pSiler(1:80, a1 = a1, a2 = a2, b1 = b1, b2 = b2, c = c, lower.tail = 0), col = "blue")

## test left-truncated samples
sims <- rSiler(n = 10000, a1 = a1, a2 = a2, b1 = b1, b2 = b2, c = c, zL = 20)
survivors <- NULL
for(i in 1:length(age)){
  survivors[i]<-(sum(sims>age[i])/10000)
}
plot(survivors, xlab = "Age", ylab = "Proportion of survivors")
lines(20:80, pSiler(20:80, a1 = a1, a2 = a2, b1 = b1, b2 = b2, c = c, zL = 20, lower.tail = 0), col = "blue")

## test right-truncated samples
sims <- rSiler(n = 10000, a1 = a1, a2 = a2, b1 = b1, b2 = b2, c = c, zU = 20)
survivors <- NULL
for(i in 1:length(age)){
  survivors[i]<-(sum(sims>age[i])/10000)
}
plot(survivors, xlab = "Age", ylab = "Proportion of survivors")
lines(1:20, pSiler(1:20, a1 = a1, a2 = a2, b1 = b1, b2 = b2, c = c, zU = 20, lower.tail = 0), col = "blue")

## test interval-truncated samples
sims <- rSiler(n = 10000, a1 = a1, a2 = a2, b1 = b1, b2 = b2, c = c, zL = 5, zU = 20)
survivors <- NULL
for(i in 1:length(age)){
  survivors[i]<-(sum(sims>age[i])/10000)
}
plot(survivors, xlab = "Age", ylab = "Proportion of survivors")
lines(5:20, pSiler(5:20, a1 = a1, a2 = a2, b1 = b1, b2 = b2, c = c, zL = 5, zU = 20, lower.tail = 0), col = "blue")

#####################################################
######          CHECK MODELS COMPILE           ######
#####################################################
 
source("../Distributions/Siler/Dist_SilerNim.R")

## check components compile
compSiler <- compileNimble(dsilerNim)
compSiler <- compileNimble(psilerNim)
compSiler <- compileNimble(qsilerNim)
compSiler <- compileNimble(rsilerNim)
 
#####################################################
######          CHECK COMPILED MODEL           ######
#####################################################
 
## check compiled model
compSiler <- compileNimble(rsilerNim)
 
# ## test non-truncated samples
sims <- sapply(1:10000, function(x) compSiler(n = 1, a1 = a1, a2 = a2, b1 = b1, b2 = b2, c = c))
hist(sims, freq = FALSE)
## check against brute-force approach
lines(seq(0, 100, length = 100), dsilerNim(seq(0, 100, length = 100), a1 = a1, a2 = a2, b1 = b1, b2 = b2, c = c), col = "blue", lwd = 2)

## test survivor function
survivors <- NULL
for(i in 1:length(age)){
  survivors[i]<-(sum(sims>age[i])/10000)
}
plot(survivors, xlab = "Age", ylab = "Proportion of survivors")
lines(1:80, psilerNim(1:80, a1 = a1, a2 = a2, b1 = b1, b2 = b2, c = c, lower.tail = 0), col = "blue")

#####################################################
######         CHECK USER INPUT ERRORS        #######
#####################################################

# Different dimensions
a1<-c(a1,a1)
x<- 12
zL <- NA
zU <- NA
dSiler(x,a1,a2,b1,b2,c,zL,zU)
x<- c(4, 13)
dSiler(x,a1,a2,b1,b2,c,zL,zU)
pSiler(x,a1,a2,b1,b2,c,zL,zU)

#zL greater than zU
a1 <- exp(-2)
zL<- 10
zU<- 3
dSiler(x,a1,a2,b1,b2,c,zL,zU)
rSiler(1,a1,a2,b1,b2,c,zL,zU)
pSiler(0.4,a1,a2,b1,b2,c,zL,zU)
qSiler(5,a1,a2,b1,b2,c,zL,zU)

#zL greater than x
zU<-NA
x <- 4
dSiler(x,a1,a2,b1,b2,c,zL,zU)
pSiler(x,a1,a2,b1,b2,c,zL,zU)

#zU less than x
zU <- 3
zL <- NA
dSiler(x,a1,a2,b1,b2,c,zL,zU)
pSiler(x,a1,a2,b1,b2,c,zL,zU)

#Negative parameter
a1<--a1
dSiler(2,a1,a2,b1,b2,c,zL,zU)
dSiler(c(1, 2),a1,a2,b1,b2,c,zL,zU)
rSiler(4,a1,a2,b1,b2,c,zL,zU)

## length argument to rSiler
a1 <- -a1
rSiler(4,a1,a2,b1,b2,c,zL,zU)
rSiler(c(2, 4),a1,a2,b1,b2,c,zL,zU)


