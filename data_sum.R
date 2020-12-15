##### Initialization #####

library(tidyverse)

##### Read in Data #####

data.0 <- read_csv("data.csv",col_types="cccccdcdd") 

##### Summarize Data #####

# Number of dates
data.0.sum <- data.0 %>%
  select(owner,location,variety,date) %>%
  distinct() %>%
  group_by(owner,location,variety) %>%
  count() %>%
  rename(dates_0=n)

# Number of dates meeting screening criteria #1
# W >= 1.0 Mg/ha for >= 3 points

data.1.list <- data.0 %>%
  filter(W>=1) %>%
  drop_na() %>%
  group_by(owner,location,variety,date) %>%
  count() %>%
  filter(n>=3) %>%
  select(-n)

data.1 <- left_join(data.1.list,
                    data.0,
                    by = c("owner","location","variety","date")) %>%
  select(owner,study,year,location,variety,rate_n_kgha,date,W,N)

data.1.sum <- data.1 %>%
  select(owner,location,variety,date) %>%
  distinct() %>%
  group_by(owner,location,variety) %>%
  count() %>%
  rename(dates_1=n)

# Number of dates meeting screening criteria #2
# sd(W) >= 1.0 Mg/ha

data.2.list <- data.1 %>% 
  group_by(owner,location,variety,date) %>%
  summarize_at(vars(W),lst(mean,sd)) %>%
  rowwise() %>%
  mutate(cv=sd/mean) %>%
  ungroup() %>%
  # filter(sd>=1.0) %>%
  # select(-c(mean,sd))
  filter(cv>=0.1) %>%
  select(-c(mean,sd,cv))

data.2 <- left_join(data.2.list,
                    data.1,
                    by = c("owner","location","variety","date")) %>%
  select(owner,study,year,location,variety,rate_n_kgha,date,W,N)

data.2.sum <- data.2 %>%
  select(owner,location,variety,date) %>%
  distinct() %>%
  group_by(owner,location,variety) %>%
  count() %>%
  rename(dates_2=n)

# Combined summary

data.sum <- data.0.sum %>%
  left_join(data.1.sum, c("owner","location","variety")) %>%
  left_join(data.2.sum, c("owner","location","variety"))

##### END #####