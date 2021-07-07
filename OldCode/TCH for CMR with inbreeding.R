library(tidyverse)
library(data.table)
rm(list=ls())

## read in data
CH<- read.table("Data/all.diag.results.txt", header = TRUE, row.names = NULL)
gene<-as.data.frame(readRDS("Data/Gene.info.rds"))

#hist(gene$f_inbreed)
## add categorical inbreeding variable
gene$Inb <- NA
for (i in 1:nrow(gene)){
  if(gene$f_inbreed[i] <= 0.2) {
    gene$Inb[i] <- 1
    #} else if (gene$f_inbreed[i] > 0.1 & gene$f_inbreed[i] < 0.2) {
    #  gene$Inb[i] <- 2
  } else {
    gene$Inb[i] <- 2
  }
}
gene$Inb <- as.factor(gene$Inb)
summary(gene$Inb)

## set date variable correctly
CH$date <- as.character(CH$date)
CH$date<- as.Date(CH$date, format="%d/%m/%Y")

## merge data frames
colnames(CH)[which(names(CH) == "tattoo")] <- "ID"
gene.data<-as.data.table(gene[,c(1,8:9)])
CH<- as.data.table(CH)
setindexv(CH, 'ID')
setindex(gene.data, 'ID')
CH <- merge(CH, gene.data, by='ID')
CH<- CH[,c(1:6,13:19,31:33)]
summary(CH)
CH<-arrange(CH, Inb)
CH<-droplevels(CH)
CH<-as.data.table(CH)

## keep only badgers with known age/sex
CH<-CH[complete.cases(CH[,8])] #sex
CH<-CH[complete.cases(CH[,10])] #age
CH<-arrange(CH, ID)

## remove badgers that are only recovered dead
counts<-CH %>% 
  group_by(ID) %>%
  summarise(no_rows = length(ID)) %>%
  tbl_df %>% print(n=nrow(.))

CH <- data.table(CH)
counts <- data.table(counts)
setindexv(CH, c('ID'))
setindexv(counts, c('ID'))
CH <- merge(CH, counts, by="ID")
#See if any counts are <=1 and it is pm
CH$toofew <- as.factor(ifelse(CH$no_rows <= 1 & CH$pm == 'Yes', 1, 0))
summary(CH$toofew)
CH <- CH[CH$toofew==0]
rm(counts)
CH <- droplevels(CH)

## create a trap month variable
CH$trap_month <- as.numeric(format(CH$date, "%m"))
CH$trap_month <- as.factor(CH$trap_month)
summary(CH$trap_month)

## put in date order
CH<- CH[order(CH$date),]

## change NAs
CH$statpak[is.na(CH$statpak)]<-3
CH$brock[is.na(CH$brock)]<-3
CH$IFNgamma[is.na(CH$IFNgamma)]<-3
CH$cult_SUM[is.na(CH$cult_SUM)]<--1

## remove occassions when no tests
#CH <- as.data.table(CH)
#CH$no.tests<- if_else(CH$statpak==3 & CH$brock==3 & CH$cult_SUM==-1 & CH$IFNgamma==3, 1, 0)
#CH<- CH[CH$no.tests==0]
#CH$no.tests<-NULL
#summary(CH)

## define what is infected
CH$infected<- as.factor(if_else((CH$statpak==1 | CH$brock==1 | CH$cult_SUM>0 | CH$IFNgamma==1),1,0))
CH$infected.as.cub<- as.factor(if_else((CH$statpak==1 | CH$brock==1 | CH$cult_SUM>0 | CH$IFNgamma==1) & CH$age=="CUB", 1, 0))
summary(CH$infected)
summary(CH$infected.as.cub)

## adjust the variable so now its ever infected or never infected
CH<- arrange(CH,ID, desc(infected))
for (i in 1:nrow(CH)){
  ifelse(CH$ID[i]==(CH$ID[i-1]),CH$infected[i]<-(CH$infected[i-1]),CH$infected[i])
}
CH$infected<- as.factor(CH$infected)

## copy infected as Cub across all capture history of each individual
CH<- arrange(CH, ID, desc(infected.as.cub))
for (i in 1:nrow(CH)){
  ifelse(CH$ID[i]==(CH$ID[i-1]),CH$infected.as.cub[i]<-(CH$infected.as.cub[i-1]),CH$infected.as.cub[i])
}
CH$infected.as.cub<- as.factor(CH$infected.as.cub)
summary(CH)

## split trap month into 4 different 'seasons', year split into 4 quarters
CH$trap_season <- NA
for (i in 1:nrow(CH)){
  if(CH$trap_month[i] == 1 | CH$trap_month[i] == 2 | CH$trap_month[i] == 3){
    CH$trap_season[i] <- 1
  } else if (CH$trap_month[i] == 4 | CH$trap_month[i] == 5 | CH$trap_month[i] == 6){
    CH$trap_season[i] <- 2
  } else if (CH$trap_month[i] == 7 | CH$trap_month[i] == 8 | CH$trap_month[i] == 9){
    CH$trap_season[i] <- 3
  } else if (CH$trap_month[i] == 10 | CH$trap_month[i] == 11 | CH$trap_month[i] == 12){
    CH$trap_season[i] <- 4
  }
}

## create occasion variable by combining capture year and trap season
CH <- unite(CH, "occasion", captureyear, trap_season, sep = ".", remove = FALSE)
CH$occasion <- as.factor(CH$occasion)
CH$occasion <- as.numeric(CH$occasion)

## adjust the occasion for 2 badger captures that are recorded as after their pm occasion
CH[8814,7] <- 81
CH[9670,7] <- 91

## create death occasion
CH$death <- ifelse(CH$pm=="Yes", CH$occasion, NA)
CH <- arrange(CH, ID, desc(pm))
for (i in 1:nrow(CH)){
  ifelse(CH$ID[i]==(CH$ID[i-1]),CH$death[i]<-CH$death[i-1], CH$death[i])
}

## create captured/recovered variable (recovered=2, captured=1) 
CH$capRec <- ifelse(CH$pm=='Yes',2,1)

## adjust duplicated badgers which were captured and died in the same season
CH <- as.data.table(CH)
CH <- arrange(CH, ID, desc(date))
CH$duplicate <- as.factor(ifelse(duplicated(CH[,c(1,7)]), 1, 0))
for (i in 1:nrow(CH)){
  ifelse((CH$ID[i]==(CH$ID[i-1]) & CH$duplicate==1), CH$occasion[i]<-CH$occasion[i]-1, CH$occasion[i])
}

## create birth occasion
CH$birth <- NA
CH <- arrange(CH, ID, date)
for (i in 1:nrow(CH)){
  if(CH$age[i] == "CUB" & CH$trap_season[i] == 1){
    CH$birth[i] <- CH$occasion[i] - 1
  } else if (CH$age[i] == "CUB" & CH$trap_season[i] == 2){
    CH$birth[i] <- CH$occasion[i] - 1
  } else if (CH$age[i] == "CUB" & CH$trap_season[i] == 3){
    CH$birth[i] <- CH$occasion[i] - 2
  } else if (CH$age[i] == "CUB" & CH$trap_season[i] == 4){
    CH$birth[i] <- CH$occasion[i] - 3
  }
}
for (i in 1:nrow(CH)){
  ifelse(CH$ID[i]==(CH$ID[i-1]),CH$birth[i]<-CH$birth[i-1], CH$birth[i])
}

CH <- droplevels(CH)

#Create status variable (infected/uninfected throughout life)
CH$status.cub <- NA
for(i in 1:nrow(CH)){
  if(CH$sex[i]=="MALE" & CH$infected.as.cub[i]==1){
    CH$status.cub[i] <- "CPM"
  } else if (CH$sex[i]=="FEMALE" & CH$infected.as.cub[i]==1){
    CH$status.cub[i] <- "CPF"
  } else if (CH$sex[i]=="MALE" & CH$infected[i]==0){
    CH$status.cub[i] <- "NPM"
  } else if (CH$sex[i]=="FEMALE" & CH$infected[i]==0){
    CH$status.cub[i] <- "NPF"
  } else {
    CH$status.cub[i] <- 0
  }
}

CH$status <- NA
for(i in 1:nrow(CH)){
  if(CH$sex[i]=="MALE" & CH$infected[i]==1){
    CH$status[i] <- "PM"
  } else if (CH$sex[i]=="FEMALE" & CH$infected[i]==1){
    CH$status[i] <- "PF"
  } else if (CH$sex[i]=="MALE" & CH$infected[i]==0){
    CH$status[i] <- "NPM"
  } else {
    CH$status[i] <- "NPF"
  }
}

CH$status <- as.factor(CH$status)
CH$status.cub <- as.factor(CH$status.cub)
summary(CH)

## create capture history
CapHist<- array(0, dim=c(length(levels(CH$ID)), length(unique(CH$occasion))))
rownames(CapHist)<- levels(CH$ID)
colnames(CapHist)<- 1:length(unique(CH$occasion))

## total number of observations
n.obs <- length(levels(CH$ID))*length(unique(CH$occasion))

## for loop to create the capture history array
CH <- as.data.table(CH)
captures <- CH[CH$capRec==1]
for(i in 1:n.obs){
  CapHist[captures$ID[i], captures$occasion[i]] <- captures$capRec[i]
}

## covariate information
CVinfo <- CH
CVinfo<- distinct(CVinfo, CVinfo$ID, .keep_all = TRUE)

## combine dataframes
CapHist <- as.data.frame(CapHist)
CapHist$ID <- rownames(CapHist)
CHCV <- merge(CVinfo, CapHist, by = 'ID')

saveRDS(CHCV, file = "Data/CHCVall.rds")
