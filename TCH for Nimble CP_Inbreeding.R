library(tidyverse)
library(data.table)
rm(list=ls())

CH<- read.table("all.diag.results.txt", header = TRUE, row.names = NULL)
gene<-as.data.table(readRDS("Gene.info.rds"))

#Categorise inbreed coefficient
gene<-arrange(gene, f_inbreed)
gene$Inb.cat<-as.factor(rep(1:3, times=c(644,644,645)))

CH$date <- as.character(CH$date)
CH$date<- as.Date(CH$date, format="%d/%m/%Y")

#Merge data frames
colnames(CH)[which(names(CH) == "tattoo")] <- "ID"
gene.data<-as.data.table(gene[,c(1,9)])
CH<- as.data.table(CH)
setindexv(CH, 'ID')
setindex(gene.data, 'ID')
CH <- merge(CH, gene.data, by='ID')

CH<- CH[,c(1:6,13:19,31:32)]
summary(CH)
CH<-droplevels(CH)

#Keep only badgers with known age/sex
CH<-CH[complete.cases(CH[,8])] #sex
CH<-CH[complete.cases(CH[,10])] #age

#Create a trap month variable
CH$trap_month <- as.numeric(format(CH$date, "%m"))
CH$trap_month <- as.factor(CH$trap_month)
summary(CH$trap_month)

#Put in date order
CH<- CH[order(CH$date),]

#Change NAs
CH$statpak[is.na(CH$statpak)]<-3
CH$brock[is.na(CH$brock)]<-3
CH$IFNgamma[is.na(CH$IFNgamma)]<-3
CH$cult_SUM[is.na(CH$cult_SUM)]<--1
#Remove occassions when no tests
CH$no.tests<- if_else(CH$statpak==3 & CH$brock==3 & CH$cult_SUM==-1 & CH$IFNgamma==3 & CH$pm!="Yes",1,0 )
CH<- CH[CH$no.tests==0]
CH$no.tests<-NULL
summary(CH)

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
CH$status<- ifelse(CH$sex=="MALE" & CH$infected==1,"Pm",0)
CH$status<- ifelse(CH$sex=="FEMALE" & CH$infected==1,"Pf",CH$status)
CH$status<- ifelse(CH$sex=="MALE" & CH$infected==0,"NPm",CH$status)
CH$status<- ifelse(CH$sex=="FEMALE" & CH$infected==0,"NPf",CH$status)
CH$status<- as.factor(CH$status)

#Create status variable (infected.as.cub/uninfected throughout life)
CH$status.cub<- ifelse(CH$sex=="MALE" & CH$infected.as.cub==1,"CPm",0)
CH$status.cub<- ifelse(CH$sex=="FEMALE" & CH$infected.as.cub==1,"CPf",CH$status.cub)
CH$status.cub<-as.factor(CH$status.cub)
summary(CH)

#See how many occasions each badger is caught
#counts<-CH %>% 
#  group_by(ID) %>%
#  summarise(no_rows = length(ID)) %>%
#  tbl_df %>% print(n=nrow(.))

#Merge count info with CH
#CH<- data.table(CH)
#counts<- data.table(counts)
#setindexv(CH, c('ID'))
#setindexv(counts, c('ID'))
#CH<- merge(CH, counts, by="ID")

#See if any counts are <=2
#CH$toofew<- as.factor(ifelse(CH$no_rows<=1 & CH$occasion==0,1,0))
#summary(CH)
#CH<-CH[CH$toofew==0]
#rm(counts)
#CH <- droplevels(CH)
#Keep only NP and CP badgers
#CH<-as.data.table(CH)
#CH<- CH[CH$status=="NPf"|CH$status=="NPm"|CH$status.cub=="CPf"|CH$status.cub=="CPm"]
summary(CH)

CH <- as.data.table(CH)
CH <- droplevels(CH)
summary(CH)

#Split trap month into 4 different 'seasons', year split into 4 quarters
CH$trap_season<- ifelse(CH$trap_month==1| CH$trap_month==2|CH$trap_month==3, "1", 0)
CH$trap_season<- ifelse(CH$trap_month==6|CH$trap_month==4|CH$trap_month==5, "2", CH$trap_season)
CH$trap_season<- ifelse(CH$trap_month==9|CH$trap_month==7|CH$trap_month==8, "3", CH$trap_season)
CH$trap_season<- ifelse(CH$trap_month==12|CH$trap_month==10|CH$trap_month==11, "4", CH$trap_season)

#Create occasion variable by combining capture year and trap season
CH <- unite(CH, "occasion", captureyear, trap_season, sep = ".", remove = FALSE)
CH$occasion <- as.factor(CH$occasion)
CH$occasion <- as.numeric(CH$occasion)

CH$death <- ifelse(CH$pm=="Yes", CH$occasion, NA)
CH<- as.data.table(CH)
CH<- arrange(CH, ID, desc(date))
for (i in 1:nrow(CH)){
  ifelse(CH$ID[i]==(CH$ID[i-1]),CH$death[i]<-CH$death[i-1], CH$death[i])
}

summary(CH)
#Get rid of redundant levels and variables
CH <- droplevels(CH)

#See how many occasions each badger is caught
#counts<-CH.CP %>% 
#  group_by(ID) %>%
#  summarise(no_rows = length(ID)) %>%
#  tbl_df %>% print(n=nrow(.))

#CH.CP<- data.table(CH.CP)
#counts<- data.table(counts)
#setindexv(CH.CP, c('ID'))
#setindexv(counts, c('ID'))
#CH.CP<- merge(CH.CP, counts, by="ID")

#See if any counts are <=2 and one of them is pm
#CH.CP$toofew<- as.factor(ifelse(CH.CP$no_rows.y<=1,1,0))
#summary(CH.CP)
#CH.CP<-CH.CP[CH.CP$toofew==0]
#rm(counts)
#CH.CP <- droplevels(CH.CP)

#Create age in quarter years variable
#Create first line for badger 1
CH<- arrange(CH, ID, occasion)
CH$age_quart.yr<- 1
#Continue for loop for all other badgers
for (i in 1:nrow(CH)){
  ifelse(CH$ID[i]!=(CH$ID[i-1]), CH$age_quart.yr[i]<-1,
         CH$age_quart.yr[i]<-((CH$age_quart.yr[i-1])+(CH$occasion[i]-CH$occasion[i-1])))
}

#Remove duplicate badgers (Same Badger caught in the same season)
CH<-as.data.table(CH)
CH$duplicate<- ifelse(duplicated(CH[,c(1,7)]), CH$duplicate<-1, CH$duplicate<- 0)
CH<-CH[CH$duplicate==0]
CH$duplicate<- NULL
summary(CH)
CH <- droplevels(CH)

#See how many occasions each badger is caught
counts<-CH %>% 
  group_by(ID) %>%
  summarise(no_rows = length(ID)) %>%
  tbl_df %>% print(n=nrow(.))

CH<- data.table(CH)
counts<- data.table(counts)
setindexv(CH, c('ID'))
setindexv(counts, c('ID'))
CH<- merge(CH, counts, by="ID")

#See if any counts are <=2 and one of them is pm
CH$toofew<- as.factor(ifelse(CH$no_rows<=1 & CH$knowndeath==1,1,0))
summary(CH)
CH<-CH[CH$toofew==0]
rm(counts)
CH <- droplevels(CH)

#CH2<- length(CH[,if(.N>2) .SD, by=ID])

#Create an age array of remaining badgers
CH$occasion<-as.factor(CH$occasion)
CP.age.array<- array(NA, dim=c(length(levels(CH$ID)), length(levels(CH$occasion))))
rownames(CP.age.array)<- levels(CH$ID)
colnames(CP.age.array)<- levels(CH$occasion)

#Total number of observations
n.obs<- length(levels(CH$ID))*length(levels(CH$occasion))
CH$trap_season<-as.numeric(CH$trap_season)

#for loop to create the age array
for(i in 1:n.obs){
  CP.age.array[CH$ID[i], CH$occasion[i]]<- CH$age_quart.yr[i]+CH$trap_season[i]-1
}

#Create a vector with occasion of first capture
CP.f <- apply(CP.age.array, 1, function(x){which(!is.na(x))[1]})

#Fill age array with ages of badgers when not observed
CP.age.array <- t(apply(CP.age.array, 1, function(x){
  n <- which(!is.na(x))[1]
  x[n:length(x)] <- x[n]:(x[n] + length(x) - n)
  x
}))

CH$birth_yr<-as.numeric(CH$occasion)-CH$age_quart.yr

#Change NA in age array to 0 (Bugs doesnt like NAs)
#CP.age.array[is.na(CP.age.array)] <- 0

#State definition
#1 = test negative on all tests
#2 = test positive on one of any tests
#3 = positive culture from one site
#4 = positive culture from multiple sites
#5 = dead
CH$disease_state <- 1
CH$disease_state <- ifelse(CH$statpak==1 | CH$brock==1,2,CH$disease_state)
CH$disease_state <- ifelse(CH$cult_SUM==1,3,CH$disease_state)
CH$disease_state <- ifelse(CH$cult_SUM>1,4,CH$disease_state)
CH$disease_state <- ifelse(CH$pm=="Yes",5, CH$disease_state)
CH$disease_state <- as.factor(CH$disease_state)
#CHeck disease state variable
summary(CH$disease_state)

#Create a state array of remaining badgers
CP.state.array<- array(0, dim=c(length(levels(CH$ID)), length(levels(CH$occasion))))
rownames(CP.state.array)<- levels(CH$ID)
colnames(CP.state.array)<- levels(CH$occasion)

#for loop to create the state array
for(i in 1:n.obs){
  CP.state.array[CH$ID[i], CH$occasion[i]]<- CH$disease_state[i]
}

#Recode 'not seen' - 0's to 5
#CP.rstatearray<- CP.state.array
#CP.rstatearray[CP.rstatearray==0] <- 5


#Setup CJS capture history, 1 for seen, 0 for unseen.
CP.CJS.array<- CP.state.array
CP.CJS.array[CP.CJS.array>=5] <- 0
CP.CJS.array[CP.CJS.array>=1] <- 1

#Covariate data
CP.CVdata<- data.table(CH)
#CP.CVdata<- arrange(CP.CVdata, ID, -as.numeric(occasion))
CP.CVdata<- distinct(CP.CVdata, CP.CVdata$ID, .keep_all = TRUE)

CP.dead <- CP.CVdata$death

#Split into inbreeding categories
summary(CP.CVdata)
low<-CP.CJS.array[1:585,]
mid<-CP.CJS.array[586:1160,]
high<-CP.CJS.array[1161:1748,]

low.dead<-CP.dead[1:585]
mid.dead<-CP.dead[586:1160]
high.dead<-CP.dead[1161:1748]

low.birth<-CP.CVdata$birth_yr[1:585]
mid.birth<-CP.CVdata$birth_yr[586:1160]
high.birth<-CP.CVdata$birth_yr[1161:1748]


saveRDS(CP.age.array, file="CP.age_array.rds")
saveRDS(CP.CVdata, file="CP.CVdata.rds")
saveRDS(CP.CJS.array, file="CP.CJS_array.rds")
saveRDS(CP.f, file = "CP.f.rds")
saveRDS(CP.dead, file="CP.dead.rds")
