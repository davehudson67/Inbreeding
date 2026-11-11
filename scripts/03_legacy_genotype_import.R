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


## Sort inbreeding information ##



## Prepare data for Microdrop
## import data
NewGeneData <- read.table("NewGeneData.txt", header = TRUE)

## change mising information to -9
NewGeneData[is.na(NewGeneData)] <- -9

## format data
colnames(NewGeneData) <- c("ID", as.character(seq(1:40)))

gene1 <- NewGeneData %>%
  select("ID", 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40)

gene2 <- NewGeneData %>%
  select("ID", 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41 )

colnames(gene1) <- c("ID", "Mel108", "Mel107", "Mel105", "Mel109", "m1", "m14", "m12", "Mel110", "Mel114", "Mel104", "Mel103",
                           "m10", "Mel101", "Mel102", "Mel112", "Mel111", "Mel115", "Mel113", "Mel117", "Mel116")

colnames(gene2) <- c("ID", "Mel108", "Mel107", "Mel105", "Mel109", "m1", "m14", "m12", "Mel110", "Mel114", "Mel104", "Mel103",
                     "m10", "Mel101", "Mel102", "Mel112", "Mel111", "Mel115", "Mel113", "Mel117", "Mel116")

## recombine data
GeneData <- bind_rows(gene1, gene2)

## arrange by ID
GeneData <- arrange(GeneData, ID)

## write df to file
write.table(GeneData, "NewGeneData.txt", sep = " ", row.names = FALSE, quote = FALSE)

########################################################################################################################################################################

## Load imputed data
ImpGeneData <- read.table("ImptD_1.txt", row.names = NULL)

## split data into loci
gene1 <- distinct(ImpGeneData, row.names, .keep_all = TRUE )
gene2 <- anti_join(ImpGeneData, gene1)

## rename loci to correct names
colnames(gene1) <- c("ID", "X1gxl.a", "X1gl.a", "X1gs.a", "X1ys.a", "m1.a", "m14.a", "m12.a", "X1yxl.a", "X2gl.a", "X1bxl.a", "X1bl.a",
                     "m10.a", "X1bs.a", "X1bm.a", "X2bxl.a", "X2bs.a", "X2gxl.a", "X2gs.a", "X2yl.a", "X2ys.a")

colnames(gene2) <- c("ID", "X1gxl.b", "X1gl.b", "X1gs.b", "X1ys.b", "m1.b", "m14.b", "m12.b", "X1yxl.b", "X2gl.b", "X1bxl.b", "X1bl.b",
                     "m10.b", "X1bs.b", "X1bm.b", "X2bxl.b", "X2bs.b", "X2gxl.b", "X2gs.b", "X2yl.b", "X2ys.b")

## combine datasets
GeneDataNew <- join(gene1, gene2, by = "ID")
rm(gene1)
rm(gene2)

## Load original gene data
gene.data<- read.table("Data/Gene_data.txt", header = TRUE)
gene.data<- gene.data[,c(1, 8:51)]

## combine
gene <- bind_rows(gene.data, GeneDataNew)

## remove columns that dont exist for all badgers
names(gene)
gene <- gene[, -c(12, 13, 44, 45)]

#########################################################################################################################################################################

## Load badger data
CH <- readRDS("Data/BadgersSGREADY220222.rds")

## combine data
CH <- inner_join(gene, CH)

## calculate inbreeding
gene <- CH[, c(2:45)]

## create matrix of indiciators for homozygous = 0 or heterozygous = 1 alleles 
inbr <- convert_raw(gene)

## Tests #######################################################################
## create allelic data frame

y <- alleles2loci(gene)
colnames(y)

# rename columns
colnames(y) <- c(paste("l", seq(1:20), sep = ""))

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
#y_red <- y[, -c(1, 3, 14)]
#t_red <- df2genind(y_red, sep = "/")
#badger.hwt_red<- hw.test(t_red)
#badger.hwt_red

## calculate measures of inbreeding...
#gene.inbr <- gene.inbr[, -c(1, 3, 14)]

## calculate MLH using inbreedR ##
het <- inbreedR::MLH(inbr)
hom <- 1 - het

hom <- data.frame(hom)

ggplot(data = hom, aes(x = hom)) +
  geom_histogram(bins = 30) +
  labs(x = "Multi-locus homozygosity", y = "Count") +
  theme(axis.text=element_text(size=12)) +
  theme(axis.title=element_text(size=15)) +
  theme(legend.text=element_text(size=10)) +
  theme(legend.title=element_text(size=15))

## calculate g2 statistic
?g2_microsats

g2 <- g2_microsats(inbr, 500, 500, CI = 0.95)
print(g2)
plot(g2)

## ad inbreeding measure to CH
CH$hom <- 1 - het
saveRDS(CH, "Data/NewGeneCH.rds")



