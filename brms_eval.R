library(tidyverse)
library(brms)
library(tidybayes)

m0 <- readRDS("Models/model-1_Bohman_Minnesota_Clearwater-DakotaRusset-Easton-RussetBurbank-Umatilla.rds")
m1 <- readRDS("Models/model-1_Bohman_Minnesota_RussetBurbank.rds")
m2 <- readRDS("Models/model-1_Bohman_Minnesota_Clearwater.rds")
m3 <- readRDS("Models/model-1_Bohman_Minnesota_Umatilla.rds")
m4 <- readRDS("Models/model-1_Bohman_Minnesota_DakotaRusset.rds")
m5 <- readRDS("Models/model-1_Bohman_Minnesota_Easton.rds")

get_variables(m1)

m1$prior

m1 %>% 
  spread_draws(b_Bmax_Intercept, r_date__Bmax[date,]) %>% 
  mutate(date_Bmax = b_Bmax_Intercept + r_date__Bmax) %>% 
  ggplot(aes(x = date_Bmax, y = date)) +
  geom_halfeyeh()

m1 %>% 
  spread_draws(b_Si_Intercept, r_date__Si[date,]) %>% 
  mutate(date_Si = b_Si_Intercept + r_date__Si) %>% 
  ggplot(aes(x = date_Si, y = date)) +
  geom_halfeyeh()


