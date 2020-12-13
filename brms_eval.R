library(tidyverse)
library(brms)
library(tidybayes)

m1.0 <- readRDS("Models/model-1_Bohman_Minnesota_Clearwater-DakotaRusset-Easton-RussetBurbank-Umatilla.rds")
m1.1 <- readRDS("Models/model-1_Bohman_Minnesota_RussetBurbank.rds")
m1.2 <- readRDS("Models/model-1_Bohman_Minnesota_Clearwater.rds")
m1.3 <- readRDS("Models/model-1_Bohman_Minnesota_Umatilla.rds")
m1.4 <- readRDS("Models/model-1_Bohman_Minnesota_DakotaRusset.rds")
m1.5 <- readRDS("Models/model-1_Bohman_Minnesota_Easton.rds")

m2.0 <- readRDS("Models/model-2_Bohman_Minnesota_Clearwater-DakotaRusset-Easton-RussetBurbank-Umatilla.rds")

m3.0 <- readRDS("Models/model-3_Bohman_Minnesota_Clearwater-DakotaRusset-Easton-RussetBurbank-Umatilla.rds")

m2.0$fit
m2.0

# the fmin() function used in Stan isn't defined in R, so we need to create it so that when we try to use brms to make predictions, it knows what to do with the fmin()
fmin <- function(x,y){
  pmin(x,y)
}

# now we can do cool stuff like take the orignal data, generate predicted outcomes, and plot those against the real values
m2.0$data %>% 
  add_predicted_draws(m2.0) %>% 
  median_hdci() %>% 
  ggplot(aes(x = W, y = N)) +
  geom_pointinterval(aes(xmin = .lower, xmax = .upper, x = .prediction), alpha = 0.5) +
  geom_point(color = "forestgreen")

# or we can do a posterior check

pp_check(m2.0, nsamples = 30)

# looks pretty good!!


# looking at group level draws --------------------------------------------

# Model 2.0

get_variables(m2.0)

m2.0$prior

m2.0 %>% 
  spread_draws(b_alpha1_Intercept, r_variety__alpha1[variety,]) %>% 
  mutate(variety_alpha1 = b_alpha1_Intercept + r_variety__alpha1) %>% 
  ggplot(aes(x = variety_alpha1, y = variety)) +
  geom_halfeyeh()

m2.0 %>% 
  spread_draws(b_alpha2_Intercept, r_variety__alpha2[variety,]) %>% 
  mutate(variety_alpha2 = b_alpha2_Intercept + r_variety__alpha2) %>% 
  ggplot(aes(x = variety_alpha2, y = variety)) +
  geom_halfeyeh()

left_join(
  m2.0 %>% 
    spread_draws(b_alpha1_Intercept, r_variety__alpha1[variety,]) %>% 
    mutate(variety_alpha1 = b_alpha1_Intercept + r_variety__alpha1),
  m2.0 %>% 
    spread_draws(b_alpha2_Intercept, r_variety__alpha2[variety,]) %>% 
    mutate(variety_alpha2 = b_alpha2_Intercept + r_variety__alpha2),
  by = c(".chain", ".iteration", ".draw", "variety")) %>% 
  ggplot(aes(x = variety_alpha1, y = variety_alpha2, color=variety)) +
  geom_point(alpha=0.1) +
  scale_color_brewer(palette = "Set1")

# Model 1.0

get_variables(m1.0)

m1.0$prior

m1.0 %>% 
  spread_draws(b_Bmax_Intercept, r_date__Bmax[date,]) %>% 
  mutate(date_Bmax = b_Bmax_Intercept + r_date__Bmax) %>% 
  ggplot(aes(x = date_Bmax, y = date)) +
  geom_halfeyeh()

m1.0 %>% 
  spread_draws(b_Si_Intercept, r_date__Si[date,]) %>% 
  mutate(date_Si = b_Si_Intercept + r_date__Si) %>% 
  ggplot(aes(x = date_Si, y = date)) +
  geom_halfeyeh()



