library(tidyverse)
library(data.table)
rm(list=ls())

CH<- read.table("Data/all.diag.results.txt", header = TRUE, row.names = NULL)
gene<-as.data.frame(readRDS("Data/Gene.info.rds"))

#hist(gene$f_inbreed)

#Categorise inbreed coefficient
gene$Inb <- NA
for (i in 1:nrow(gene)){
if(gene$f_inbreed[i] <= 0.15) {
  gene$Inb[i] <- 1
#} else if (gene$f_inbreed[i] > 0.1 & gene$f_inbreed[i] < 0.2) {
#  gene$Inb[i] <- 2
} else {
  gene$Inb[i] <- 2
}
}
gene$Inb <- as.factor(gene$Inb)
summary(gene$Inb)

#Set date variable correctly
CH$date <- as.character(CH$date)
CH$date<- as.Date(CH$date, format="%d/%m/%Y")

#Merge data frames
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

#Keep only badgers with known age/sex
CH<-CH[complete.cases(CH[,8])] #sex
CH<-CH[complete.cases(CH[,10])] #age

#Create a trap month variable
CH$trap_month <- as.numeric(format(CH$date, "%m"))
CH$trap_month <- as.factor(CH$trap_month)
summary(CH$trap_month)

#Put in date order
CH<- CH[order(CH$date),]
CH <- CH[CH$captureyear>=1996]

#Change NAs
CH$statpak[is.na(CH$statpak)]<-3
CH$brock[is.na(CH$brock)]<-3
CH$IFNgamma[is.na(CH$IFNgamma)]<-3
CH$cult_SUM[is.na(CH$cult_SUM)]<--1

#Remove occassions when no tests
#CH$no.tests<- if_else(CH$statpak==3 & CH$brock==3 & CH$cult_SUM==-1 & CH$IFNgamma==3 & CH$pm!="Yes",1,0 )
#CH<- CH[CH$no.tests==0]
#CH$no.tests<-NULL
#summary(CH)

#Define what is infected
CH$infected<- as.factor(if_else((CH$statpak==1 | CH$brock==1 | CH$cult_SUM>0 | CH$IFNgamma==1),1,0))
CH$infected.as.cub<- as.factor(if_else((CH$statpak==1 | CH$brock==1 | CH$cult_SUM>0 | CH$IFNgamma==1) & CH$age=="CUB", 1, 0))

summary(CH$infected)
summary(CH$infected.as.cub)

#Adjust the variable so now its ever infected or never infected
CH<- arrange(CH,ID, desc(infected))
for (i in 1:nrow(CH)){
  ifelse(CH$ID[i]==(CH$ID[i-1]),CH$infected[i]<-(CH$infected[i-1]),CH$infected[i])
}
CH$infected<- as.factor(CH$infected)

#Copy infected as Cub across all capture history of each individual
CH<- arrange(CH, ID, desc(infected.as.cub))
for (i in 1:nrow(CH)){
  ifelse(CH$ID[i]==(CH$ID[i-1]),CH$infected.as.cub[i]<-(CH$infected.as.cub[i-1]),CH$infected.as.cub[i])
}
CH$infected.as.cub<- as.factor(CH$infected.as.cub)
summary(CH)

#Add if known death
CH<- arrange(CH, ID, desc(date))
CH$knowndeath <- ifelse(CH$pm=="Yes",1,0)
for (i in 1:nrow(CH)){
  ifelse(CH$ID[i]==(CH$ID[i-1]),CH$knowndeath[i]<-(CH$knowndeath[i-1]),CH$infected.as.cub[i])
}

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

#Split data
CH <- as.data.table(CH)
CH.CP<-as.data.frame(CH[CH$status.cub=="CPM" | CH$status.cub=="CPF"])
CH.P<-as.data.frame(CH[CH$status=="PM" | CH$status=="PF"])
CH.NP<-as.data.frame(CH[CH$status=="NPM" | CH$status=="NPF"])

## Postiive Badgers
CH.P <- droplevels(CH.P)
summary(CH.P)

#Split trap month into 4 different 'seasons', year split into 4 quarters
CH.P$trap_season <- NA
for (i in 1:nrow(CH.P)){
  if(CH.P$trap_month[i] == 1 | CH.P$trap_month[i] == 2 | CH.P$trap_month[i] == 3){
    CH.P$trap_season[i] <- 1
  } else if (CH.P$trap_month[i] == 4 | CH.P$trap_month[i] == 5 | CH.P$trap_month[i] == 6){
    CH.P$trap_season[i] <- 2
  } else if (CH.P$trap_month[i] == 7 | CH.P$trap_month[i] == 8 | CH.P$trap_month[i] == 9){
    CH.P$trap_season[i] <- 3
  } else if (CH.P$trap_month[i] == 10 | CH.P$trap_month[i] == 11 | CH.P$trap_month[i] == 12){
    CH.P$trap_season[i] <- 4
  }
}

#Create occasion variable by combining capture year and trap season
CH.P <- unite(CH.P, "occasion", captureyear, trap_season, sep = ".", remove = FALSE)
CH.P$occasion <- as.factor(CH.P$occasion)
CH.P$occasion <- as.numeric(CH.P$occasion)

CH.P$death <- ifelse(CH.P$pm=="Yes", CH.P$occasion, NA)
CH.P<- arrange(CH.P, ID, desc(date))
for (i in 1:nrow(CH.P)){
  ifelse(CH.P$ID[i]==(CH.P$ID[i-1]),CH.P$death[i]<-CH.P$death[i-1], CH.P$death[i])
}

CH.P <- arrange(CH.P, ID, date)
for (i in 1:nrow(CH.P)){
  if(CH.P$age[i] == "CUB" & CH.P$trap_season[i] == 1){
    CH.P$birth[i] <- CH.P$occasion[i]
  } else if (CH.P$age[i] == "CUB" & CH.P$trap_season[i] == 2){
    CH.P$birth[i] <- CH.P$occasion[i] - 1
  } else if (CH.P$age[i] == "CUB" & CH.P$trap_season[i] == 3){
    CH.P$birth[i] <- CH.P$occasion[i] - 2
  } else if (CH.P$age[i] == "CUB" & CH.P$trap_season[i] == 4){
    CH.P$birth[i] <- CH.P$occasion[i] - 3
  }
}

for (i in 1:nrow(CH.P)){
  ifelse(CH.P$ID[i]==(CH.P$ID[i-1]),CH.P$birth[i]<-CH.P$birth[i-1], CH.P$birth[i])
}

#Get rid of redundant levels and variables
CH.P <- droplevels(CH.P)

#Create age in quarter years variable
#Create first line for badger 1
#CH.P<- arrange(CH.P, ID, occasion)
#CH.P$age_quart.yr<- 1
#Continue for loop for all other badgers
#for (i in 1:nrow(CH.P)){
#  ifelse(CH.P$ID[i]!=(CH.P$ID[i-1]), CH.P$age_quart.yr[i]<-1,
#         CH.P$age_quart.yr[i]<-((CH.P$age_quart.yr[i-1])+(CH.P$occasion[i]-CH.P$occasion[i-1])))
#}

#Remove duplicate badgers (Same Badger caught in the same season)
CH.P<-as.data.table(CH.P)
CH.P$duplicate<- ifelse(duplicated(CH.P[,c(1,7)]), CH.P$duplicate<-1, CH.P$duplicate<- 0)
CH.P<-CH.P[CH.P$duplicate==0]
CH.P$duplicate<- NULL
summary(CH.P)
CH.P <- droplevels(CH.P)

#See how many occasions each badger is caught
#counts<-CH.P %>% 
#  group_by(ID) %>%
#  summarise(no_rows = length(ID)) %>%
#  tbl_df %>% print(n=nrow(.))

#CH.P<- data.table(CH.P)
#counts<- data.table(counts)
#setindexv(CH.P, c('ID'))
#setindexv(counts, c('ID'))
#CH.P<- merge(CH.P, counts, by="ID")

#See if any counts are <=2 and one of them is pm
#CH.P$toofew<- as.factor(ifelse(CH.P$no_rows<=1 & CH.P$knowndeath==1,1,0))
#summary(CH.P)
#CH.P<-CH.P[CH.P$toofew==0]
#rm(counts)
#CH.P <- droplevels(CH.P)

#CH.P2<- length(CH.P[,if(.N>2) .SD, by=ID])

#Create an age array of remaining badgers
#CH.P$occasion<-as.factor(CH.P$occasion)
#P.age.array<- array(NA, dim=c(length(levels(CH.P$ID)), length(levels(CH.P$occasion))))
#rownames(P.age.array)<- levels(CH.P$ID)
#colnames(P.age.array)<- levels(CH.P$occasion)

#Total number of observations
n.obs<- length(levels(CH.P$ID))*length(unique(CH.P$occasion))
#CH.P$trap_season<-as.numeric(CH.P$trap_season)

#for loop to create the age array
#for(i in 1:n.obs){
#  P.age.array[CH.P$ID[i], CH.P$occasion[i]]<- CH.P$age_quart.yr[i]+CH.P$trap_season[i]-1
#}

#Create a vector with occasion of first capture
#P.f <- apply(P.age.array, 1, function(x){which(!is.na(x))[1]})

#Fill age array with ages of badgers when not observed
#P.age.array <- t(apply(P.age.array, 1, function(x){
#  n <- which(!is.na(x))[1]
#  x[n:length(x)] <- x[n]:(x[n] + length(x) - n)
#  x
#}))

#CH.P$birth_yr<-as.numeric(CH.P$occasion)-CH.P$age_quart.yr

#Change NA in age array to 0 (Bugs doesnt like NAs)
#P.age.array[is.na(P.age.array)] <- 0

#State definition
#1 = test negative on all tests
#2 = test positive on one of any tests
#3 = positive culture from one site
#4 = positive culture from multiple sites
#5 = dead

CH.P$disease_state <- NA
for (i in 1:nrow(CH.P)){
  if(CH.P$statpak[i] == 1 | CH.P$brock[i] == 1 | CH.P$IFNgamma[i] == 1){
    CH.P$disease_state[i] <- 2
  } else if (CH.P$cult_SUM[i] == 1){
    CH.P$disease_state[i] <- 3
  } else if (CH.P$cult_SUM[i] > 1){
    CH.P$disease_state[i] <- 4
  } else if (CH.P$pm[i] == "Yes"){
    CH.P$disease_state[i] <- 5
  } else {
    CH.P$disease_state[i] <- 1
  }
}

CH.P$disease_state <- as.factor(CH.P$disease_state)

#Check disease state variable
summary(CH.P$disease_state)

#Order badgers into levels of mlh
#CH.P<-arrange(CH.P, Inb, ID)

#Create a state array of remaining badgers
P.state.array<- array(0, dim=c(length(levels(CH.P$ID)), length(unique(CH.P$occasion))))
rownames(P.state.array)<- levels(CH.P$ID)
colnames(P.state.array)<- 1:length(unique(CH.P$occasion))

#for loop to create the state array
for(i in 1:n.obs){
  P.state.array[CH.P$ID[i], CH.P$occasion[i]] <- CH.P$disease_state[i]
}

#Recode 'not seen' - 0's to 5
#P.rstatearray<- P.state.array
#P.rstatearray[P.rstatearray==0] <- 5

#Setup CJS capture history, 1 for seen, 0 for unseen.
P.CJS.array<- P.state.array
P.CJS.array[P.CJS.array>=5] <- 0
P.CJS.array[P.CJS.array>=1] <- 1
#order<-rownames(P.CJS.array)
#P.CJS.array<-as.data.frame(P.CJS.array)
#P.CJS.array$ID <- order
#P.CJS.array$ID<-order

#Covariate data
P.CVdata<- data.table(CH.P)
#P.CVdata<- arrange(P.CVdata, desc(Inb), ID)
P.CVdata<- distinct(P.CVdata, P.CVdata$ID, .keep_all = TRUE)

## combine dataframes
P.CJS.array <- as.data.frame(P.CJS.array)
P.CJS.array$ID <- rownames(P.CJS.array)
CH.P <- merge(P.CJS.array, P.CVdata, by='ID')

#Ensure all data is in the same order
#library(wrapr)
#p <- match_order(P.CJS.array$ID, P.CVdata$ID)
#P.CJS.array <- P.CJS.array[p,,drop=FALSE]
#P.age.array <- P.age.array[p,,drop=FALSE]
#P.dead <- P.CVdata$death

#saveRDS(P.age.array, file="P.age_array.rds")
#saveRDS(P.CVdata, file="P.CVdata.rds")
#saveRDS(P.CJS.array, file="P.CJS_array.rds")
#saveRDS(P.f, file = "P.f.rds")
#saveRDS(P.dead, file="P.dead.rds")

#Split into inbreeding categories
#summary(P.CVdata)
#low<-P.CJS.array[1:550,]
#mid<-P.CJS.array[192:576,]
#high<-P.CJS.array[551:617,]

#low.dead<-P.dead[1:550]
#mid.dead<-P.dead[192:576]
#high.dead<-P.dead[551:617]

#low.birth<-P.CVdata$birth_yr[1:550]
#mid.birth<-P.CVdata$birth_yr[192:576]
#high.birth<-P.CVdata$birth_yr[551:617]

saveRDS(CH.P, file="Data/CH.P.rds")
#saveRDS(P.CVdata, file="P.CVdata.rds")
#saveRDS(P.CJS.array, file="P.CJS_array.rds")
#saveRDS(P.f, file = "P.f.rds")
#saveRDS(P.dead, file="P.dead.rds")

#######################################################################################################
## Never Postiive Badgers
CH.NP <- droplevels(CH.NP)
summary(CH.NP)

#Split trap month into 4 different 'seasons', year split into 4 quarters
CH.NP$trap_season <- NA
for (i in 1:nrow(CH.NP)){
  if(CH.NP$trap_month[i] == 1 | CH.NP$trap_month[i] == 2 | CH.NP$trap_month[i] == 3){
    CH.NP$trap_season[i] <- 1
  } else if (CH.NP$trap_month[i] == 4 | CH.NP$trap_month[i] == 5 | CH.NP$trap_month[i] == 6){
    CH.NP$trap_season[i] <- 2
  } else if (CH.NP$trap_month[i] == 7 | CH.NP$trap_month[i] == 8 | CH.NP$trap_month[i] == 9){
    CH.NP$trap_season[i] <- 3
  } else if (CH.NP$trap_month[i] == 10 | CH.NP$trap_month[i] == 11 | CH.NP$trap_month[i] == 12){
    CH.NP$trap_season[i] <- 4
  }
}

#Create occasion variable by combining capture year and trap season
CH.NP <- unite(CH.NP, "occasion", captureyear, trap_season, sep = ".", remove = FALSE)
CH.NP$occasion <- as.factor(CH.NP$occasion)
CH.NP$occasion <- as.numeric(CH.NP$occasion)

CH.NP$death <- ifelse(CH.NP$pm=="Yes", CH.NP$occasion, NA)
CH.NP<- arrange(CH.NP, ID, desc(date))
for (i in 1:nrow(CH.NP)){
  ifelse(CH.NP$ID[i]==(CH.NP$ID[i-1]),CH.NP$death[i]<-CH.NP$death[i-1], CH.NP$death[i])
}

CH.NP <- arrange(CH.NP, ID, date)
for (i in 1:nrow(CH.NP)){
  if(CH.NP$age[i] == "CUB" & CH.NP$trap_season[i] == 1){
    CH.NP$birth[i] <- CH.NP$occasion[i]
  } else if (CH.NP$age[i] == "CUB" & CH.NP$trap_season[i] == 2){
    CH.NP$birth[i] <- CH.NP$occasion[i] - 1
  } else if (CH.NP$age[i] == "CUB" & CH.NP$trap_season[i] == 3){
    CH.NP$birth[i] <- CH.NP$occasion[i] - 2
  } else if (CH.NP$age[i] == "CUB" & CH.NP$trap_season[i] == 4){
    CH.NP$birth[i] <- CH.NP$occasion[i] - 3
  }
}

for (i in 1:nrow(CH.NP)){
  ifelse(CH.NP$ID[i]==(CH.NP$ID[i-1]),CH.NP$birth[i]<-CH.NP$birth[i-1], CH.NP$birth[i])
}


#Get rid of redundant levels and variables
CH.NP <- droplevels(CH.NP)

#Remove duplicate badgers (Same Badger caught in the same season)
CH.NP<-as.data.table(CH.NP)
CH.NP$duplicate<- ifelse(duplicated(CH.NP[,c(1,7)]), CH.NP$duplicate<-1, CH.NP$duplicate<- 0)
CH.NP<-CH.NP[CH.NP$duplicate==0]
CH.NP$duplicate<- NULL
summary(CH.NP)
CH.NP <- droplevels(CH.NP)

#Total number of observations
n.obs<- length(levels(CH.NP$ID))*length(unique(CH.NP$occasion))

CH.NP$disease_state <- NA
for (i in 1:nrow(CH.NP)){
  if(CH.NP$statpak[i] == 1 | CH.NP$brock[i] == 1 | CH.NP$IFNgamma[i] == 1){
    CH.NP$disease_state[i] <- 2
  } else if (CH.NP$cult_SUM[i] == 1){
    CH.NP$disease_state[i] <- 2
  } else if (CH.NP$cult_SUM[i] > 1){
    CH.NP$disease_state[i] <- 3
  } else if (CH.NP$pm[i] == "Yes"){
    CH.NP$disease_state[i] <- 5
  } else {
    CH.NP$disease_state[i] <- 1
  }
}

CH.NP$disease_state <- as.factor(CH.NP$disease_state)

#Check disease state variable
summary(CH.NP$disease_state)

#Create a state array of remaining badgers
P.state.array<- array(0, dim=c(length(levels(CH.NP$ID)), length(unique(CH.NP$occasion))))
rownames(P.state.array)<- levels(CH.NP$ID)
colnames(P.state.array)<- 1:length(unique(CH.NP$occasion))

#for loop to create the state array
for(i in 1:n.obs){
  P.state.array[CH.NP$ID[i], CH.NP$occasion[i]] <- CH.NP$disease_state[i]
}

#Setup CJS capture history, 1 for seen, 0 for unseen.
P.CJS.array<- P.state.array
P.CJS.array[P.CJS.array>=5] <- 0
P.CJS.array[P.CJS.array>=1] <- 1

#Covariate data
P.CVdata<- data.table(CH.NP)
P.CVdata<- distinct(P.CVdata, P.CVdata$ID, .keep_all = TRUE)

## combine dataframes
P.CJS.array <- as.data.frame(P.CJS.array)
P.CJS.array$ID <- rownames(P.CJS.array)

CH.NP <- merge(P.CJS.array, P.CVdata, by='ID')
saveRDS(CH.NP, file="Data/CH.NP.rds")

##########################################################################################################
## Cub Postiive Badgers
CH.CP <- droplevels(CH.CP)
summary(CH.CP)

#Split trap month into 4 different 'seasons', year split into 4 quarters
CH.CP$trap_season <- NA
for (i in 1:nrow(CH.CP)){
  if(CH.CP$trap_month[i] == 1 | CH.CP$trap_month[i] == 2 | CH.CP$trap_month[i] == 3){
    CH.CP$trap_season[i] <- 1
  } else if (CH.CP$trap_month[i] == 4 | CH.CP$trap_month[i] == 5 | CH.CP$trap_month[i] == 6){
    CH.CP$trap_season[i] <- 2
  } else if (CH.CP$trap_month[i] == 7 | CH.CP$trap_month[i] == 8 | CH.CP$trap_month[i] == 9){
    CH.CP$trap_season[i] <- 3
  } else if (CH.CP$trap_month[i] == 10 | CH.CP$trap_month[i] == 11 | CH.CP$trap_month[i] == 12){
    CH.CP$trap_season[i] <- 4
  }
}

#Create occasion variable by combining capture year and trap season
CH.CP <- unite(CH.CP, "occasion", captureyear, trap_season, sep = ".", remove = FALSE)
CH.CP$occasion <- as.factor(CH.CP$occasion)
CH.CP$occasion <- as.numeric(CH.CP$occasion)

CH.CP$death <- ifelse(CH.CP$pm=="Yes", CH.CP$occasion, NA)
CH.CP<- arrange(CH.CP, ID, desc(date))
for (i in 1:nrow(CH.CP)){
  ifelse(CH.CP$ID[i]==(CH.CP$ID[i-1]),CH.CP$death[i]<-CH.CP$death[i-1], CH.CP$death[i])
}

CH.CP <- arrange(CH.CP, ID, date)
for (i in 1:nrow(CH.CP)){
  if(CH.CP$age[i] == "CUB" & CH.CP$trap_season[i] == 1){
    CH.CP$birth[i] <- CH.CP$occasion[i]
  } else if (CH.CP$age[i] == "CUB" & CH.CP$trap_season[i] == 2){
    CH.CP$birth[i] <- CH.CP$occasion[i] - 1
  } else if (CH.CP$age[i] == "CUB" & CH.CP$trap_season[i] == 3){
    CH.CP$birth[i] <- CH.CP$occasion[i] - 2
  } else if (CH.CP$age[i] == "CUB" & CH.CP$trap_season[i] == 4){
    CH.CP$birth[i] <- CH.CP$occasion[i] - 3
  }
  
}

for (i in 1:nrow(CH.CP)){
  ifelse(CH.CP$ID[i]==(CH.CP$ID[i-1]),CH.CP$birth[i]<-CH.CP$birth[i-1], CH.CP$birth[i])
}

#Get rid of redundant levels and variables
CH.CP <- droplevels(CH.CP)

#Remove duplicate badgers (Same Badger caught in the same season)
CH.CP<-as.data.table(CH.CP)
CH.CP$duplicate<- ifelse(duplicated(CH.CP[,c(1,7)]), CH.CP$duplicate<-1, CH.CP$duplicate<- 0)
CH.CP<-CH.CP[CH.CP$duplicate==0]
CH.CP$duplicate<- NULL
summary(CH.CP)
CH.CP <- droplevels(CH.CP)

#Total number of observations
n.obs<- length(levels(CH.CP$ID))*length(unique(CH.CP$occasion))

CH.CP$disease_state <- NA
for (i in 1:nrow(CH.CP)){
  if(CH.CP$statpak[i] == 1 | CH.CP$brock[i] == 1 | CH.CP$IFNgamma[i] == 1){
    CH.CP$disease_state[i] <- 2
  } else if (CH.CP$cult_SUM[i] == 1){
    CH.CP$disease_state[i] <- 2
  } else if (CH.CP$cult_SUM[i] > 1){
    CH.CP$disease_state[i] <- 3
  } else if (CH.CP$pm[i] == "Yes"){
    CH.CP$disease_state[i] <- 5
  } else {
    CH.CP$disease_state[i] <- 1
  }
}

CH.CP$disease_state <- as.factor(CH.CP$disease_state)

#Check disease state variable
summary(CH.CP$disease_state)

#Create a state array of remaining badgers
P.state.array<- array(0, dim=c(length(levels(CH.CP$ID)), length(unique(CH.CP$occasion))))
rownames(P.state.array)<- levels(CH.CP$ID)
colnames(P.state.array)<- 1:length(unique(CH.CP$occasion))

#for loop to create the state array
for(i in 1:n.obs){
  P.state.array[CH.CP$ID[i], CH.CP$occasion[i]] <- CH.CP$disease_state[i]
}

#Setup CJS capture history, 1 for seen, 0 for unseen.
P.CJS.array<- P.state.array
P.CJS.array[P.CJS.array>=5] <- 0
P.CJS.array[P.CJS.array>=1] <- 1

#Covariate data
P.CVdata<- data.table(CH.CP)
P.CVdata<- distinct(P.CVdata, P.CVdata$ID, .keep_all = TRUE)

## combine dataframes
P.CJS.array <- as.data.frame(P.CJS.array)
P.CJS.array$ID <- rownames(P.CJS.array)

CH.CP <- merge(P.CJS.array, P.CVdata, by='ID')
saveRDS(CH.CP, file="Data/CH.CP.rds")

