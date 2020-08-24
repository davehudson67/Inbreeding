library(tidyverse)
library(adegenet)
library(data.table)
library(ape)
library(pegas)
library(seqinr)

rm(list=ls())

#Add data ----------------------------------------------------------------------------------------------------------

gene.data<- read.table("Gene_data.txt", header = TRUE)
row.names(gene.data) <- gene.data$ID
gene.data.info<- gene.data[1:7]
gene.data<- gene.data[,8:51]

#row.names(gene.data)<-gene.data$ID
#remove extra info
#gene.data<-gene.data[2:37]
#Create allelic data frame

y<- alleles2loci(gene.data)
colnames(y)

#Rename columns
colnames(y)<- c("X1bl",  "X1bm",  "X1bs",  "X1bxl", "X1gl",  "X1gm",  "X1gs",  "X1gxl", "X1ys",  "X1yxl", "X2bs",
                "X2bxl", "X2gl", "X2gs",  "X2gxl", "X2yl",  "X2ys",  "m1",    "m10",  "m12",   "m14",   "m15") 

#create genind object
t<-df2genind(y, sep = "/")
t
summary(t)

#bartlett test for homogeneity of variances
t2 <- summary(t)
bartlett.test(list(t2$Hexp, t2$Hobs))

#test for Hardy-Weinberg equilibrium
badger.hwt<- hw.test(t)
badger.hwt

#Compute the mean inbreeding for each individual...
temp<- inbreeding(t, N = 100)
class(temp)
head(names(temp))
head(temp[[1]],20)

# temp is a list of values sampled from the likelihood distribution of each individual; 
# means values are obtained for all individuals using sapply:

Fbar <- sapply(temp, mean)
hist(Fbar, col="firebrick", main= "Average inbreeding in badgers")
gene.data.info$f_inbreed<-Fbar

saveRDS(gene.data.info, file="Gene.info.rds")
 