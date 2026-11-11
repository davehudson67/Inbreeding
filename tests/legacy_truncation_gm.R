## load library
library(nimble)
library(purrr)
library(lamW)

#####################################################
######              CHECK R CODE               ######
#####################################################

## load random sampler
source("../Distributions/GompertzMakeham/Dist_GompzMake.R")
source("../Distributions/GompertzMakeham/Dist_GompertzMakehamNim.R")
a <- 0.05
b <- 0.02
c <- exp(-5)

## test non-truncated samples
sims <- rGompzMake(n = 10000, a = a, b = b, c = c)
hist(sims, freq = FALSE)
## add PDF
lines(0:100, dGompzMake(0:100, a = a, b = b, c = c), col = "blue", lwd = 2)

## test left-truncated samples
sims <- rGompzMake(n = 10000, a = a, b = b, c = c, zL = 20)
hist(sims, freq = FALSE)
## check against brute-force approach
sims1 <- NULL
while(length(sims1) < 10000) {
  temp <- rGompzMake(10000, a = a, b = b, c = c)
  temp <- temp[temp >= 20]
  sims1 <- c(sims1, temp)
}
sims1 <- sims1[1:10000]
hist(sims1, add = TRUE, col = "red", freq = FALSE)
## add PDF
lines(20:100, dGompzMake(20:100, a = a, b = b, c = c, zL = 20), col = "blue", lwd = 2)

## test right-truncated samples
sims <- rGompzMake(n = 10000, a = a, b = b, c = c, zU = 5)
hist(sims, freq = FALSE)
## check against brute-force approach
sims1 <- NULL
while(length(sims1) < 10000) {
  temp <- rGompzMake(10000, a = a, b = b, c = c)
  temp <- temp[temp <= 5]
  sims1 <- c(sims1, temp)
}
sims1 <- sims1[1:10000]
hist(sims1, add = TRUE, col = "red", freq = FALSE)
## add PDF
lines(0:5, dGompzMake(0:5, a = a, b = b, c = c, zU = 5), col = "blue", lwd = 2)

## test interval-truncated samples
sims <- rGompzMake(n = 10000, a = a, b = b, c = c, zL = 3, zU = 5)
hist(sims, freq = FALSE)
## check against brute-force approach
sims1 <- NULL
while(length(sims1) < 10000) {
  temp <- rGompzMake(10000, a = a, b = b, c = c)
  temp <- temp[temp >= 3 & temp <= 5]
  sims1 <- c(sims1, temp)
}
sims1 <- sims1[1:10000]
hist(sims1, add = TRUE, col = "red", freq = FALSE)
## add PDF
lines(3:5, dGompzMake(3:5, a = a, b = b, c = c, zL = 3, zU = 5), col = "blue", lwd = 2)

## check user input NA in R code
sims <- rGompzMake(n = 10000, a = a, b = b, c = c, zL = 3, zU = NA)
summary(sims)

## check Survival times =======================================================================================================
age <- 1:80

## test non-truncated samples
sims <- rGompzMake(n = 10000, a = a, b = b, c = c)
survivors <- NULL
for(i in 1:length(age)){
  survivors[i]<-(sum(sims>age[i])/10000)
}
plot(survivors, xlab = "Age", ylab = "Proportion of survivors")
lines(1:80, pGompzMake(1:80, a = a, b = b, c = c, lower.tail = 0), col = "blue")

## test left-truncated samples
sims <- rGompzMake(n = 10000, a = a, b = b, c = c, zL = 20)
survivors <- NULL
for(i in 1:length(age)){
  survivors[i]<-(sum(sims>age[i])/10000)
}
plot(survivors, xlab = "Age", ylab = "Proportion of survivors")
lines(20:80, pGompzMake(20:80, a = a, b = b, c = c, zL = 20, lower.tail = 0), col = "blue")

## test right-truncated samples
sims <- rGompzMake(n = 10000, a = a, b = b, c = c, zU = 20)
survivors <- NULL
for(i in 1:length(age)){
  survivors[i]<-(sum(sims>age[i])/10000)
}
plot(survivors, xlab = "Age", ylab = "Proportion of survivors")
lines(1:20, pGompzMake(1:20, a = a, b = b, c = c, zU = 20, lower.tail = 0), col = "blue")

## test interval-truncated samples
sims <- rGompzMake(n = 10000, a = a, b = b, c = c, zL = 5, zU = 20)
survivors <- NULL
for(i in 1:length(age)){
  survivors[i]<-(sum(sims>age[i])/10000)
}
plot(survivors, xlab = "Age", ylab = "Proportion of survivors")
lines(5:20, pGompzMake(5:20, a = a, b = b, c = c, zL = 5, zU = 20, lower.tail = 0), col = "blue")

#####################################################
######          CHECK MODELS COMPILE           ######
#####################################################

## check components compile
compGompzMake <- compileNimble(dgompzMakeNim)
compGompzMake <- compileNimble(pgompzMakeNim)
compGompzMake <- compileNimble(qgompzMakeNim)
compGompzMake <- compileNimble(rgompzMakeNim)

#####################################################
######          CHECK COMPILED MODEL           ######
#####################################################

## check compiled model
compGompertzMakeNim <- compileNimble(rgompzMakeNim)
 
## test non-truncated samples
sims <- sapply(1:10000, function(x) compGompertzMakeNim(n = 1, a = a, b = b, c = c))
hist(sims, freq = FALSE)
## check against brute-force approach
lines(seq(0, 100, length = 100), dgompzMakeNim(seq(0, 100, length = 100), a = a, b = b, c = c), col = "blue", lwd = 2)

## test survivor function
survivors <- NULL
for(i in 1:length(age)){
  survivors[i]<-(sum(sims>age[i])/10000)
}
plot(survivors, xlab = "Age", ylab = "Proportion of survivors")
lines(1:80, pgompzMakeNim(1:80, a = a, b = b, c = c, lower.tail = 0), col = "blue")

#####################################################
######         CHECK USER INPUT ERRORS        #######
#####################################################

# Different dimensions
a<-c(a,a)
x<- 12
zL <- NA
zU <- NA
dGompzMake(x,a,b,c,zL,zU)
x<- c(4, 13)
dGompzMake(x,a,b,c,zL,zU)
pGompzMake(x,a,b,c,zL,zU)

#zL greater than zU
a <- exp(-2)
zL<- 10
zU<- 3
dGompzMake(x,a,b,c,zL,zU)
rGompzMake(1,a,b,c,zL,zU)
pGompzMake(0.4,a,b,c,zL,zU)
qGompzMake(5,a,b,c,zL,zU)

#zL greater than x
zU<-NA
x <- 4
dGompzMake(x,a,b,c,zL,zU)
pGompzMake(x,a,b,c,zL,zU)

#zU less than x
zU <- 3
zL <- NA
dGompzMake(x,a,b,c,zL,zU)
pGompzMake(x,a,b,c,zL,zU)

#Negative parameter
a<--a
dGompzMake(2,a,b,c,zL,zU)
dGompzMake(c(1, 2),a,b,c,zL,zU)
rGompzMake(4,a,b,c,zL,zU)

## length argument to rGompzMake
a <- -a
rGompzMake(4,a,b,c,zL,zU)
rGompzMake(c(2, 4),a,b,c,zL,zU)

