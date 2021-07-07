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

gene.data<- read.table("Data/Gene_data.txt", header = TRUE)
row.names(gene.data) <- gene.data$ID
info<- gene.data[1:7]
gene.data<- gene.data[,8:51]

## create matrix of indiciators for homozygous = 0 or heterozygous = 1 alleles 
gene.inbr <- convert_raw(gene.data)
#colnames(gene.inbr)

## calculate measures of inbreeding, leaving one marker out each time...

## calculate MLH using inbreedR ##
hom1 <- 1 - inbreedR::MLH(gene.inbr[,-1])
hom2 <- 1 - inbreedR::MLH(gene.inbr[,-2])
hom3 <- 1 - inbreedR::MLH(gene.inbr[,-3])
hom4 <- 1 - inbreedR::MLH(gene.inbr[,-4])
hom5 <- 1 - inbreedR::MLH(gene.inbr[,-5])
hom6 <- 1 - inbreedR::MLH(gene.inbr[,-6])
hom7 <- 1 - inbreedR::MLH(gene.inbr[,-7])
hom8 <- 1 - inbreedR::MLH(gene.inbr[,-8])
hom9 <- 1 - inbreedR::MLH(gene.inbr[,-9])
hom10 <- 1 - inbreedR::MLH(gene.inbr[,-10])
hom11 <- 1 - inbreedR::MLH(gene.inbr[,-11])
hom12 <- 1 - inbreedR::MLH(gene.inbr[,-12])
hom13 <- 1 - inbreedR::MLH(gene.inbr[,-13])
hom14 <- 1 - inbreedR::MLH(gene.inbr[,-14])
hom15 <- 1 - inbreedR::MLH(gene.inbr[,-15])
hom16 <- 1 - inbreedR::MLH(gene.inbr[,-16])
hom17 <- 1 - inbreedR::MLH(gene.inbr[,-17])
hom18 <- 1 - inbreedR::MLH(gene.inbr[,-18])
hom19 <- 1 - inbreedR::MLH(gene.inbr[,-19])
hom20 <- 1 - inbreedR::MLH(gene.inbr[,-20])
hom21 <- 1 - inbreedR::MLH(gene.inbr[,-21])
hom22 <- 1 - inbreedR::MLH(gene.inbr[,-22])

inbr_markerdropped <- data.frame(hom1, hom2, hom3, hom4, hom5, hom6, hom7, hom8, hom9, hom10,
                                 hom11, hom12, hom13, hom14, hom15, hom16, hom17, hom18, hom19, hom20, hom21, hom22)

saveRDS(inbr_markerdropped, "inbr_markerdropped.rds")
