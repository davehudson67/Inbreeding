library(tidyverse)
library(data.table)
rm(list=ls())

CH<- read.table("all.diag.results.txt", header = TRUE, row.names = NULL)
gene<-as.data.frame(readRDS("Gene.info.rds"))

hist(gene$f_inbreed)
#Categorise inbreed coefficient
gene$Inb <- NA
for (i in 1:nrow(gene)){
if(gene$f_inbreed[i] <= 0.18) {
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
gene.data<-as.data.table(gene[,c(1,9)])
CH<- as.data.table(CH)
setindexv(CH, 'ID')
setindex(gene.data, 'ID')
CH <- merge(CH, gene.data, by='ID')
CH<- CH[,c(1:6,13:19,31:32)]
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
#summary(CH)

#Split data
CH <- as.data.table(CH)
CH.CP<-as.data.frame(CH[CH$status.cub=="CPm" | CH$status.cub=="CPf"])
CH.P<-as.data.frame(CH[CH$status=="Pm" | CH$status=="Pf"])

CH.P <- droplevels(CH.P)
summary(CH.P)

#Split trap month into 4 different 'seasons', year split into 4 quarters
CH.P$trap_season<- ifelse(CH.P$trap_month==1| CH.P$trap_month==2|CH.P$trap_month==3, "1", 0)
CH.P$trap_season<- ifelse(CH.P$trap_month==6|CH.P$trap_month==4|CH.P$trap_month==5, "2", CH.P$trap_season)
CH.P$trap_season<- ifelse(CH.P$trap_month==9|CH.P$trap_month==7|CH.P$trap_month==8, "3", CH.P$trap_season)
CH.P$trap_season<- ifelse(CH.P$trap_month==12|CH.P$trap_month==10|CH.P$trap_month==11, "4", CH.P$trap_season)

#Create occasion variable by combining capture year and trap season
CH.P <- unite(CH.P, "occasion", captureyear, trap_season, sep = ".", remove = FALSE)
CH.P$occasion <- as.factor(CH.P$occasion)
CH.P$occasion <- as.numeric(CH.P$occasion)

CH.P$death <- ifelse(CH.P$pm=="Yes", CH.P$occasion, NA)
CH.P<- as.data.table(CH.P)
CH.P<- arrange(CH.P, ID, desc(date))
for (i in 1:nrow(CH.P)){
  ifelse(CH.P$ID[i]==(CH.P$ID[i-1]),CH.P$death[i]<-CH.P$death[i-1], CH.P$death[i])
}

summary(CH.P)
#Get rid of redundant levels and variables
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
#CH.P$toofew<- as.factor(ifelse(CH.P$no_rows.y<=1,1,0))
#summary(CH.P)
#CH.P<-CH.P[CH.P$toofew==0]
#rm(counts)
#CH.P <- droplevels(CH.P)

#Create age in quarter years variable
#Create first line for badger 1
CH.P<- arrange(CH.P, ID, occasion)
CH.P$age_quart.yr<- 1
#Continue for loop for all other badgers
for (i in 1:nrow(CH.P)){
  ifelse(CH.P$ID[i]!=(CH.P$ID[i-1]), CH.P$age_quart.yr[i]<-1,
         CH.P$age_quart.yr[i]<-((CH.P$age_quart.yr[i-1])+(CH.P$occasion[i]-CH.P$occasion[i-1])))
}

#Remove duplicate badgers (Same Badger caught in the same season)
CH.P<-as.data.table(CH.P)
CH.P$duplicate<- ifelse(duplicated(CH.P[,c(1,7)]), CH.P$duplicate<-1, CH.P$duplicate<- 0)
CH.P<-CH.P[CH.P$duplicate==0]
CH.P$duplicate<- NULL
summary(CH.P)
CH.P <- droplevels(CH.P)

#See how many occasions each badger is caught
counts<-CH.P %>% 
  group_by(ID) %>%
  summarise(no_rows = length(ID)) %>%
  tbl_df %>% print(n=nrow(.))

CH.P<- data.table(CH.P)
counts<- data.table(counts)
setindexv(CH.P, c('ID'))
setindexv(counts, c('ID'))
CH.P<- merge(CH.P, counts, by="ID")

#See if any counts are <=2 and one of them is pm
CH.P$toofew<- as.factor(ifelse(CH.P$no_rows<=1 & CH.P$knowndeath==1,1,0))
summary(CH.P)
CH.P<-CH.P[CH.P$toofew==0]
rm(counts)
CH.P <- droplevels(CH.P)

#CH.P2<- length(CH.P[,if(.N>2) .SD, by=ID])

#Create an age array of remaining badgers
CH.P$occasion<-as.factor(CH.P$occasion)
P.age.array<- array(NA, dim=c(length(levels(CH.P$ID)), length(levels(CH.P$occasion))))
rownames(P.age.array)<- levels(CH.P$ID)
colnames(P.age.array)<- levels(CH.P$occasion)

#Total number of observations
n.obs<- length(levels(CH.P$ID))*length(levels(CH.P$occasion))
CH.P$trap_season<-as.numeric(CH.P$trap_season)

#for loop to create the age array
for(i in 1:n.obs){
  P.age.array[CH.P$ID[i], CH.P$occasion[i]]<- CH.P$age_quart.yr[i]+CH.P$trap_season[i]-1
}

#Create a vector with occasion of first capture
P.f <- apply(P.age.array, 1, function(x){which(!is.na(x))[1]})

#Fill age array with ages of badgers when not observed
P.age.array <- t(apply(P.age.array, 1, function(x){
  n <- which(!is.na(x))[1]
  x[n:length(x)] <- x[n]:(x[n] + length(x) - n)
  x
}))

CH.P$birth_yr<-as.numeric(CH.P$occasion)-CH.P$age_quart.yr

#Change NA in age array to 0 (Bugs doesnt like NAs)
#P.age.array[is.na(P.age.array)] <- 0

#State definition
#1 = test negative on all tests
#2 = test positive on one of any tests
#3 = positive culture from one site
#4 = positive culture from multiple sites
#5 = dead
CH.P$disease_state <- 1
CH.P$disease_state <- ifelse(CH.P$statpak==1 | CH.P$brock==1,2,CH.P$disease_state)
CH.P$disease_state <- ifelse(CH.P$cult_SUM==1,3,CH.P$disease_state)
CH.P$disease_state <- ifelse(CH.P$cult_SUM>1,4,CH.P$disease_state)
CH.P$disease_state <- ifelse(CH.P$pm=="Yes",5, CH.P$disease_state)
CH.P$disease_state <- as.factor(CH.P$disease_state)
#CH.Peck disease state variable
summary(CH.P$disease_state)

#Order badgers into levels of mlh
CH.P<-arrange(CH.P, Inb, ID)

#Create a state array of remaining badgers
P.state.array<- array(0, dim=c(length(levels(CH.P$ID)), length(levels(CH.P$occasion))))
rownames(P.state.array)<- levels(CH.P$ID)
colnames(P.state.array)<- levels(CH.P$occasion)

#for loop to create the state array
for(i in 1:n.obs){
  P.state.array[CH.P$ID[i], CH.P$occasion[i]]<- CH.P$disease_state[i]
}

#Recode 'not seen' - 0's to 5
#P.rstatearray<- P.state.array
#P.rstatearray[P.rstatearray==0] <- 5


#Setup CJS capture history, 1 for seen, 0 for unseen.
P.CJS.array<- P.state.array
P.CJS.array[P.CJS.array>=5] <- 0
P.CJS.array[P.CJS.array>=1] <- 1
order<-rownames(P.CJS.array)
P.CJS.array<-as.data.frame(P.CJS.array)
P.CJS.array$ID <- order

#P.CJS.array$ID<-order

#Covariate data
P.CVdata<- data.table(CH.P)
P.CVdata<- arrange(P.CVdata, Inb, ID)
P.CVdata<- distinct(P.CVdata, P.CVdata$ID, .keep_all = TRUE)

#Ensure all data is in the same order
library(wrapr)
p <- match_order(P.CJS.array$ID, P.CVdata$ID)
P.CJS.array <- P.CJS.array[p,,drop=FALSE]
P.age.array <- P.age.array[p,,drop=FALSE]

P.dead <- P.CVdata$death

#saveRDS(P.age.array, file="P.age_array.rds")
#saveRDS(P.CVdata, file="P.CVdata.rds")
#saveRDS(P.CJS.array, file="P.CJS_array.rds")
#saveRDS(P.f, file = "P.f.rds")
#saveRDS(P.dead, file="P.dead.rds")

#Split into inbreeding categories
summary(P.CVdata)
low<-P.CJS.array[1:550,]
#mid<-P.CJS.array[192:576,]
high<-P.CJS.array[551:617,]

low.dead<-P.dead[1:550]
#mid.dead<-P.dead[192:576]
high.dead<-P.dead[551:617]

low.birth<-P.CVdata$birth_yr[1:550]
#mid.birth<-P.CVdata$birth_yr[192:576]
high.birth<-P.CVdata$birth_yr[551:617]


saveRDS(P.age.array, file="P.age_array.rds")
saveRDS(P.CVdata, file="P.CVdata.rds")
saveRDS(P.CJS.array, file="P.CJS_array.rds")
saveRDS(P.f, file = "P.f.rds")
saveRDS(P.dead, file="P.dead.rds")
