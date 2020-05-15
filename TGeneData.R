library(tidyverse)
library(adegenet)
library(data.table)
library(ape)
library(pegas)
library(seqinr)

rm(list=ls())

#Add data ----------------------------------------------------------------------------------------------------------

gene.data<- read.table("Gene_data.txt", header = TRUE)
gene.data.info<- gene.data[1:7]
gene.data<- gene.data[,c(1,16:51)]
row.names(gene.data)<-gene.data$ID
#remove extra info
gene.data<-gene.data[2:37]
#Create allelic data frame
y<- alleles2loci(gene.data)
y

#Rename columns
colnames(y)<- c("X1gl",  "X1gm",  "X1gs",  "X1gxl", "X1ys",  "X1yxl", "X2bs",  "X2bxl", "X2gl", 
                "X2gs",  "X2gxl", "X2yl", "X2ys",  "m1",    "m10",   "m12",   "m14",   "m15") 
#Create genind object
y.gen.ind<- loci2genind(y)

is.genind(y.gen.ind)

y.gen.ind@tab

temp <- inbreeding(y.gen.ind, res.type = )
head(names(temp))
Fbar<- sapply(temp, mean)
hist(Fbar, col="firebrick", main= "Average inbreeding in badgers")

toto<-summary(y.gen.ind)
toto

gene.data.info$Inbreed<-Fbar
rm(temp)
rm(toto)
rm(y)
rm(gene.data)
