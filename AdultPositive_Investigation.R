## load data
CH <- readRDS("Data/BadgersforCH300621.rds")
load("Data/inbreed_data.RData")

## combine data
CH <- inner_join(CH, mlh)

## select desired data
Adult_pos <- CH
  filter(infected_as_cub != infected_lifetime)

Adult_pos <- Adult_pos %>%
  group_by(ID) %>%
  mutate(Total_positive_tests = sum(ifelse(Cult_Sum[Cult_Sum > 0], 1, 0), brock[brock > 0], v1s[v1s > 0], v2s[v2s > 0], statpak[statpak > 0], IFNgamma[IFNgamma > 0])) %>%
  mutate(Total_positive_captures = sum(infected_now)) %>%
  mutate(First_positive_test = min(occ[infected_now == 1])) %>%
  select(ID, date, pm, occ, age, statpak, Cult_Sum, brock, v1s, v2s, IFNgamma, captures, infected_now, Total_positive_captures, Total_positive_tests, First_positive_test, 
         Number_of_positive_tests_today, Number_of_tests_today) %>%
  ungroup()
  
Adult_pos <- Adult_pos %>%
  group_by(ID) %>%
  arrange(occ) %>%
  mutate(Number_of_captures_after_first_positive = ifelse(occ > First_positive_test, 1, 0)) %>%
  mutate(Number_of_captures_after_first_positive = sum(Number_of_captures_after_first_positive)) %>%
  mutate(After_first_positive = ifelse(occ > First_positive_test, 1, 0))


Adult_pos <- Adult_pos %>%
  arrange(ID) %>%
  mutate(Number_of_positive_tests_after_first_positive_occ = sum(Number_of_positive_tests_today[After_first_positive == 1])) %>%
  mutate(Number_of_tests_after_first_positive_occ = sum(Number_of_tests_today[After_first_positive == 1])) %>%
  mutate(Proportion_of_positive_tests_after_first_positive = Number_of_positive_tests_after_first_positive_occ/Number_of_tests_after_first_positive_occ)

Adult_pos <- Adult_pos %>%
  group_by(ID) %>%
  mutate(Age_at_first_positive = ifelse(First_positive_test == occ, age, 0)) %>%
  mutate(Age_at_first_positive = max(Age_at_first_positive))


summary(Adult_pos$Age_at_first_positive)

ggplot(Adult_pos) +
  geom_bar(aes(x = as.factor(Age_at_first_positive)))


Adult_pos <- distinct(Adult_pos, ID, .keep_all = TRUE)


## plots
ggplot(Adult_pos) +
  geom_bar(aes(x = as.factor(captures), fill = "firebrick")) +
  labs(title = "Captures") +
  theme(legend.position = "none")

Adult_pos %>%
  filter(Number_of_captures_after_first_positive > 0) %>%
  ggplot() +
  geom_histogram(aes(x = Proportion_of_positive_tests_after_first_positive), binwidth = 0.05, fill = "firebrick", position = "dodge") +
  labs(title = "Proportion of positive tests after first testing positive", subtitle = "Badgers who were recaptured once or more after first test postive") +
  theme(legend.position = "none")

Adult_pos %>%
  filter(captures > 0) %>%
  ggplot() +
  geom_bar(aes(x = as.factor(captures), fill = "firebrick")) +
  labs(title = "Captures") +
  theme(legend.position = "none")





  
Adult_pos$Number_of_captures_after_first_positive <- as.factor(Adult_pos$Number_of_captures_after_first_positive)  
summary(Adult_pos)  
  
  

check$Always_positive_after_first_positive <- as.factor(check$Always_positive_after_first_positive)
summary(check)  




group_by(ID) %>%
  arrange(occ) %>%
  


ggplot(check) +
  geom_bar( aes(x = check$Always_positive_after_first_positive))

min(hisTest$captures_after_positive)

check$Always_positive_after_first_positive <- as.factor(check$Always_positive_after_first_positive)
summary(check)


check <- Test %>%
  distinct(ID, .keep_all = TRUE) %>%
  group_by(ID) %>%
  filter(captures_after_positive == 1) %>%
  mutate(Always_positive_after_first_positive = ifelse(Number_of_Positives  < captures_after_positive, 0, 1))


hisTest$Always_positive_after_first_positive <- as.factor(hisTest$Always_positive_after_first_positive)
summary(hisTest)


filter(Number_of_Positives < captures)

dev.off()



