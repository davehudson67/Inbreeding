library(tidyverse)
library(lubridate)

rm(list=ls())

## read in data
CH <- read_csv("all.diag.results.csv")

## keep only required variables
CH <- select(CH, tattoo, date, pm, captureyear, sex, age, age_years)#, statpak, IFNgamma, cult_SUM, brock)

## set date variable correctly
CH$date <- dmy(CH$date)

## adjust columns
colnames(CH)[which(names(CH) == "tattoo")] <- "ID"

## keep only badgers with known age
CH <- CH %>%
  drop_na(age_years)

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

## create a trap month variable
CH <- CH %>%
  mutate(trap_month = as.numeric(month(date)))

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

CH$occ <- CH %>%
  arrange(date) %>%
  group_by(occasion) %>%
  group_indices()

## adjust occasion to allow for the start of the study/badgers can be born before the first capture occasion
CH$occ <- CH$occ + 5

## create birth occasion (need to decide on the +1 to determine when most births happen)
CH <- CH %>%
  group_by(ID) %>%
  mutate(birth = if_else(age == "CUB", occ, 1000)) %>%
  mutate(birth = min(birth)) %>%
  mutate(birth = birth - (min(trap_season[age == "CUB"])) + 1) %>%
  ungroup()

## check error message
filter(CH, is.infinite(birth))

## adjust birth for this individual
CH$birth[CH$ID=="UN574"] <- 154 - 4

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
CH$death[CH$ID == "UN088"] <- 35
#CH$max_captures[CH$ID == "UN088"] <- 0
#CH$captures[CH$ID == "UN088"] <- 0
CH$death[CH$ID == "UN089"] <- 35
#CH$max_captures[CH$ID == "UN089"] <- 0
#CH$captures[CH$ID == "UN089"] <- 0
CH$death[CH$ID == "UN135"] <- 63
#CH$max_captures[CH$ID == "UN135"] <- 0
#CH$captures[CH$ID == "UN135"] <- 0
CH$death[CH$ID == "UN136"] <- 63
#CH$max_captures[CH$ID == "UN136"] <- 0
#CH$captures[CH$ID == "UN136"] <- 0


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

## keep unique entries
CH <- distinct(CH, ID, .keep_all = TRUE)

## save output
saveRDS(CH, "CH.rds")


CH[52,]
122-119
CH %>%
  filter(captures > last_seen - birth)


CH <- CH %>%
  group_by(ID) %>%
  mutate(captures = if_else(captures > last_seen - birth, captures - 1, captures))
