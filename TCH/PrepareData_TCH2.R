library(tidyverse)
library(lubridate)

rm(list=ls())

## read in data
CH <- readRDS("Data/BadgersCombinedData300621.rds")

CH <- arrange(CH, tattoo)

## keep only required variables
CH <- select(CH, tattoo, date, pm, occ, capture_yr, sex, age, year_fc, age_fc_yrs, trap_season, statpak, GAMMA, Cult_Sum, brock, v1s, v2s)

## adjust columns
colnames(CH)[which(names(CH) == "tattoo")] <- "ID"

## create infected variables
## change NAs
#CH$statpak[is.na(CH$statpak)] <- 3
#CH$brock[is.na(CH$brock)] <- 3
#CH$IFNgamma[is.na(CH$IFNgamma)] <- 3
#CH$cult_SUM[is.na(CH$cult_SUM)] <- -1

#CH <- CH %>%
#  mutate(infected = if_else((CH$statpak == 1 | CH$brock == 1 | CH$cult_SUM > 0 | CH$IFNgamma ==1 ), 1, 0)) %>%
#  mutate(infected_as_cub = if_else((CH$statpak == 1 | CH$brock == 1 | CH$cult_SUM > 0 | CH$IFNgamma == 1) & CH$age == "CUB", 1, 0))
#  mutate(no_test = if_else((CH$statpak == 3 & CH$brock == 3 & CH$cult_SUM == -1 & CH$IFNgamma == 3), 1, 0))

#CH <- CH %>%
#  group_by(ID) %>%
#  mutate(infected = max(infected)) %>%
#  mutate(infected_as_cub = max(infected_as_cub)) %>%
#  ungroup()

## adjust occasion to allow for the start of the study/badgers can be born before the first capture occasion
#CH$occ <- CH$occ + 5

## create birth occasion (need to decide on the +1 to determine when most births happen)
CH$occ <- as.numeric(CH$occ)

CH <- CH %>%
  group_by(ID) %>%
  mutate(birth = ifelse(age_fc_yrs == 1, occ - 4, occ)) %>%
  mutate(birth = min(birth)) %>%
  mutate(bts = min(trap_season[age == 0 | age == 1])) %>%
  mutate(birth = birth - bts + 1) %>%
  ungroup()

## check error message
#filter(CH, is.infinite(birth))

## create death occasion
CH <- CH %>%
  group_by(ID) %>%
  mutate(death = if_else(pm == "Yes", occ, 0)) %>%
  mutate(death = max(death)) %>%
  ungroup()

CH$death[CH$death == 0] <- 1000

## check each individual only has 1 pm date
CH %>%
  group_by(ID) %>%
  filter(pm == "Yes") %>%
  count(sort = TRUE)

## remove duplicate captures... individuals captured twice in same season.
#CH <- CH %>%
#  arrange(ID, pm) %>%
#  distinct(ID, occ, .keep_all = TRUE)

##-------------------------------------------- SPLIT DATA -----------------------------------------------##

## select individuals whose only entry is a recovered dead pm
CHdr <- CH %>%
  group_by(ID) %>%
  add_count(name = "captures") %>%
  filter(captures == 1 & pm == "Yes") %>%
  mutate(captures = captures - 1)

## remove all CHdr entries from CH
CH <- anti_join(CH, CHdr)

## create last_seen variable for this group
CHdr$last_seen <- CHdr$birth

## create max possible captures for this group
CHdr <- CHdr %>%
  group_by(ID) %>%
  mutate(max_captures = death - birth - 1) %>%
  mutate(max_captures = if_else(max_captures < 0, 0, max_captures)) %>%
  ungroup()

##------------------------------------------- DEAL WITH CH --------------------------------------------##

## check any individuals where birth >= death
CH %>%
  filter(birth >= death)

## add last seen alive
CH <- CH %>%
  group_by(ID) %>%
  mutate(last_seen = min(max(occ), death - 1))

## add number of captures 
CH <- CH %>%
  group_by(ID) %>%
  mutate(captures = n_distinct(occ)) %>%
  ungroup()

## create max possible captures (-1 removes death recovery occasion)
eos <- max(CH$occ)
CH <- CH %>%
  group_by(ID) %>%
  mutate(max_captures = if_else(death == 1000, eos - birth, death - birth - 1)) %>%
  mutate(max_captures = if_else(max_captures < 0, 0, max_captures)) %>%
  ungroup()

##------------------------------------------- REJOIN DATA -----------------------------------------------##

## combine dfs
CH <- bind_rows(CH, CHdr)
rm(CHdr)

##-------------------------------------------INFECTION STATUS-------------------------------------------##

## sort infection variables
CH$statpak <- as.numeric(CH$statpak)
CH$IFNgamma <- as.numeric(as.factor(CH$GAMMA)) - 1

CH$Number_of_positive_tests_today <- NULL

for(i in 1:nrow(CH)){
  CH$Number_of_positive_tests_today[i] <- sum(CH$statpak[i], CH$IFNgamma[i], CH$brock[i], CH$v1s[i], CH$v2s[i], CH$Cult_Sum[i], na.rm = TRUE)
}


CH$IFNgamma[is.na(CH$IFNgamma)] <- 0
#CH$IFNgamma <- CH$IFNgamma
CH$statpak[is.na(CH$statpak)] <- -1
CH$brock[is.na(CH$brock)] <- -1
CH$v1s[is.na(CH$v1s)] <- -1
CH$v2s[is.na(CH$v2s)] <- -1
CH$Cult_Sum[is.na(CH$Cult_Sum)] <- -1

for(i in 1:nrow(CH)){
  CH$Number_of_tests_today[i] <- sum(ifelse(CH$statpak[i] >= 0, 1, 0), ifelse(CH$IFNgamma[i] >= 0, 1, 0), 
                                     ifelse(CH$brock[i] >= 0, 1, 0), ifelse(CH$v1s[i] >= 0, 1, 0), ifelse(CH$v2s[i] >= 0, 1, 0), 
                                     ifelse(CH$Cult_Sum[i] == 0, 1, 0), ifelse(CH$Cult_Sum[i] > 0, CH$Cult_Sum[i], 0))
}



CH <- CH %>%
  group_by(ID) %>%
  mutate(infected_now = ifelse(statpak == 1 | IFNgamma == 1 | Cult_Sum > 0 | brock == 1 | v1s == 1 | v2s == 1, 1, 0)) %>%
  mutate(infected2 = ifelse(sum(statpak[statpak > 0], IFNgamma[IFNgamma > 0], Cult_Sum[Cult_Sum > 0], brock[brock > 0], v1s[v1s > 0], v2s[v2s > 0]) > 1, 1, 0)) %>%
 # mutate(infected2 = max(infected2)) %>%
 # mutate(infected_lifetime = max(infected)) %>%
  mutate(infected_as_cub = ifelse(infected_now == 1 & age <= 1, 1, 0)) %>%
 # mutate(infected_as_cub = max(infected_as_cub)) %>%
  ungroup()

CH <- CH %>%
  group_by(ID) %>%
  mutate(infected2 = max(infected2)) %>%
  mutate(infected_lifetime = max(infected_now)) %>%
  mutate(infected_as_cub = max(infected_as_cub)) %>%
  ungroup()


##----------------------------------------- run some checks ---------------------------------------------##
## check for any entries made after death
CH %>%
  group_by(ID) %>%
  filter(last_seen > death)

## any birth = death
CH %>%
  group_by(ID) %>%
  filter(birth == death)

## adjust these individuals
CH$death[CH$ID == "UN088"] <- 39
#CH$max_captures[CH$ID == "UN088"] <- 0
#CH$captures[CH$ID == "UN088"] <- 0
CH$death[CH$ID == "UN089"] <- 39
#CH$max_captures[CH$ID == "UN089"] <- 0
#CH$captures[CH$ID == "UN089"] <- 0
CH$death[CH$ID == "UN135"] <- 67
#CH$max_captures[CH$ID == "UN135"] <- 0
#CH$captures[CH$ID == "UN135"] <- 0
CH$death[CH$ID == "UN136"] <- 67
#CH$max_captures[CH$ID == "UN136"] <- 0
#CH$captures[CH$ID == "UN136"] <- 0
CH$death[CH$ID == "UN574"] <- 159


#CH$death[CH$ID == "UN129"] <- 56
#CH$max_captures[CH$ID == "UN129"] <- 0
#CH$captures[CH$ID == "UN136"] <- 0
#CH$death[CH$ID == "L081"] <- 40
#CH$max_captures[CH$ID == "L081"] <- 0
#CH$captures[CH$ID == "L081"] <- 0
#CH$last_seen[CH$ID == "L081"] <- 39

CH %>%
  group_by(ID) %>%
  filter(last_seen < birth)

CH %>%
  group_by(ID) %>%
  filter(captures > max_captures)

CH <- CH %>%
  group_by(ID) %>%
  mutate(captures = if_else(captures > max_captures, max_captures, captures))

## any last_seen = death
CH %>%
  group_by(ID) %>%
  filter(last_seen == death)

## change unknown death to NA
CH$death[CH$death == 1000] <- NA

## save out for capture history
saveRDS(CH, "Data/BadgersforCH300621.rds")

## keep unique entries
CH <- distinct(CH, ID, .keep_all = TRUE)

## check number of captures, these individuals were captured in their birth season so need to remove one season
CH %>%
  group_by(ID) %>%
  filter(last_seen - birth < captures)

CH <- CH %>%
  group_by(ID) %>%
  mutate(captures = if_else(last_seen - birth < captures, captures - 1, captures))

CH <- select(CH, ID, pm, sex, year_fc, birth, death, last_seen, captures, max_captures, infected2, infected_as_cub, infected_lifetime)

## save output
saveRDS(CH, "Data/BadgersNewREADY300621.rds")


summary(CH)
CH$infected_lifetime <- as.factor(CH$infected_lifetime)
CH$infected_as_cub <- as.factor(CH$infected_as_cub)
CH$infected <- as.factor(CH$infected)
CH$infected2 <- as.factor(CH$infected2)


