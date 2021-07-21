library(tidyverse)
library(adegenet)
library(data.table)
library(ape)
library(pegas)
library(seqinr)
library(Rhh)
library(inbreedR)
library(plyr)

rm(list=ls())

#Add data ----------------------------------------------------------------------------------------------------------

CH <- readRDS("Data/BadgersNewREADY300621.rds")
gene.data<- read.table("Data/Gene_data.txt", header = TRUE)

## combine data
CH <- inner_join(gene.data, CH)
gene.data <- select(CH, 1:51)

row.names(gene.data) <- gene.data$ID
info<- gene.data[1:7]
gene.data<- gene.data[,8:51]

## create matrix of indiciators for homozygous = 0 or heterozygous = 1 alleles 
gene.inbr <- convert_raw(gene.data)
gene.inbr$ID <- rownames(gene.data)

## Tests #######################################################################
## create allelic data frame

gene.data$ID <- rownames(gene.data)
gene.data <- select(gene.data, 45, 1:44 )
y <- alleles2loci(gene.data[,-1])
colnames(y)
## rename columns
colnames(y)<- c(paste("ms", seq(1:22), sep = ""))

## create genind object
t <- df2genind(y, sep = "/")
summary(t)

#bartlett test for homogeneity of variances
t2 <- summary(t)
bartlett.test(list(t2$Hexp, t2$Hobs))

#test for Hardy-Weinberg equilibrium
badger.hwt<- hw.test(t)
badger.hwt

##remove 1, 3 and 14
y_red <- y[, -c(1,14)]
t_red <- df2genind(y_red, sep = "/")
badger.hwt_red<- hw.test(t_red)
badger.hwt_red


## calculate measures of inbreeding...
gene.inbr <- gene.inbr[, -c(1,14)]

## calculate MLH using inbreedR ##
gene.inbr <- select(gene.inbr, -ID)
het <- inbreedR::MLH(gene.inbr)
hom <- 1 - het
hist(hom)

?g2_microsats
g2 <- g2_microsats(gene.inbr, 100, 100, CI = 0.95)
print(g2)
plot(g2)

hf <- r2_hf(gene.inbr, nboot = 100)
plot(hf)


## Heterozygosity-heterozygosity correlation
HHC_badgers <- HHC(genes, reps = 1000)
HHC_badgers

## HFC parameters
hf <- r2_hf(genes, nboot = 100, type = "msats")
hf
plot(hf)













## heterozygosity/heterozygosity correlation
h <- h_cor("gene.info", na.string = NULL, n = 100, method = "sh")
mean(h)
quantile(h, probs = c(0.025, 0.975))
hist(h)



## separate into vectors so each can be used as categorical variable in analysis
inb1 <- gene.inbr[, 1]
inb2 <- gene.inbr[, 2]
inb3 <- gene.inbr[, 3]
inb4 <- gene.inbr[, 4]
inb5 <- gene.inbr[, 5]
inb6 <- gene.inbr[, 6]
inb7 <- gene.inbr[, 7]
inb8 <- gene.inbr[, 8]
inb9 <- gene.inbr[, 9]
inb10 <- gene.inbr[, 10]
inb11 <- gene.inbr[, 11]
inb12 <- gene.inbr[, 12]
inb13 <- gene.inbr[, 13]
inb14 <- gene.inbr[, 14]
inb15 <- gene.inbr[, 15]
inb16 <- gene.inbr[, 16]
inb17 <- gene.inbr[, 17]
inb18 <- gene.inbr[, 18]
inb19 <- gene.inbr[, 19]
inb20 <- gene.inbr[, 20]
inb21 <- gene.inbr[, 21]
inb22 <- gene.inbr[, 22]

## calculate measures of inbreeding...

## calculate MLH using inbreedR ##
het <- inbreedR::MLH(gene.inbr)
hom <- 1 - het
hist(hom)

## calculate measures using Rhh ##
## add ID back to dataframe
gene.data$ID <- rownames(gene.data)
gene.data <- select(gene.data, 45, 1:44 )
write_delim(gene.data, "gene.info", delim = " ")

## create the three individual multilocus heterozygosity estimates
mlh <- Rhh::mlh("gene.info", "mlh_out", na.string = NULL, n.digits = 3)
mlh <- mlh[-1, ]
mlh$ID <- as.character(mlh$ID)
mlh$IR <- mlh$IR[,1]
mlh$SH <- mlh$SH[,1]
mlh$HL <- mlh$HL[,1]
summary(mlh)

## calculate using adgenet and pegas package ##

## create allelic data frame
y <- alleles2loci(gene.data[,-1])
colnames(y)
## rename columns
colnames(y)<- c("X1bl",  "X1bm",  "X1bs",  "X1bxl", "X1gl",  "X1gm",  "X1gs",  "X1gxl", "X1ys",  "X1yxl", "X2bs",
                "X2bxl", "X2gl", "X2gs",  "X2gxl", "X2yl",  "X2ys",  "m1",    "m10",  "m12",   "m14",   "m15") 

## create genind object
t <- df2genind(y, sep = "/")
summary(t)

## compute the mean inbreeding for each individual
temp <- adegenet::inbreeding(t, N = 200)

## temp is a list of values sampled from the likelihood distribution of each individual; 
## means values are obtained for all individuals using sapply:
Fbar <- sapply(temp, mean)
summary(Fbar)
hist(Fbar, col="firebrick", main= "Average inbreeding in badgers")

## check correlations between measures of homozygozity
mlh$Fbar <- Fbar
mlh$hom <- hom
mlh[,-1] %>%
  as.data.frame() %>%
  ggpairs()

## save for output analysis ##
#rm(gene.data, gene.inbr, info, t, temp, y, Fbar, het, hom)
#save.image("inbreed_data.RData")

#################################################################################################

## Tests ##

#bartlett test for homogeneity of variances
t2 <- summary(t)
bartlett.test(list(t2$Hexp, t2$Hobs))

#test for Hardy-Weinberg equilibrium
badger.hwt<- hw.test(t)
badger.hwt


## heterozygosity/heterozygosity correlation
h <- h_cor("gene.info", na.string = NULL, n = 100, method = "sh")
mean(h)
quantile(h, probs = c(0.025, 0.975))
hist(h)




## Now using adgenet and pegas package ##

#Create allelic data frame
y <- alleles2loci(gene.data[,-1])

y <- y[,-14]
colnames(y)

#Rename columns
colnames(y)<- c("X1bl",  "X1bm",  "X1bs",  "X1bxl", "X1gl",  "X1gm",  "X1gs",  "X1gxl", "X1ys",  "X1yxl", "X2bs",
                "X2bxl", "X2gl", "X2gs",  "X2gxl", "X2yl",  "X2ys",  "m1",    "m10",  "m12",   "m14",   "m15") 

#create genind object
t <- df2genind(y, sep = "/")
t
summary(t)

#bartlett test for homogeneity of variances
t2 <- summary(t)
bartlett.test(list(t2$Hexp, t2$Hobs))

#test for Hardy-Weinberg equilibrium
badger.hwt<- hw.test(t)
badger.hwt

#Compute the mean inbreeding for each individual...
temp <- inbreeding(t, N = 100)
#class(temp)
#head(names(temp))
#head(temp[[1]],20)

# temp is a list of values sampled from the likelihood distribution of each individual; 
# means values are obtained for all individuals using sapply:
Fbar <- sapply(temp, mean)
summary(Fbar)
hist(Fbar, col="firebrick", main= "Average inbreeding in badgers")
gene.inbr$f <- Fbar

#saveRDS(gene.inbr, file="Data/Gene.info.rds")

## Now using InbreedR package ##

## data.frame with individuals in rows and loci in columns, containing genotypes coded as 0 (homozygote), 1 (heterozygote) and NA (missing)

genes <- gene.inbr[, -c(23,24)]

## calculate g2 - identity disequilibrium





#colnames(gene.inbr) <- colnames(y)
gene.inbr$het <- het

gene.inbr <- join(gene.inbr, mlh, by = "ID")

saveRDS(gene.inbr, "gene_data.rds")


## het vs inbreeding
cor.test(het, mlh$HL)


