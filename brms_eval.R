library(tidyverse)
library(brms)
library(tidybayes)
library(shinystan)

m1.0 <- readRDS("Models/model-1_Bohman_Minnesota_Clearwater-DakotaRusset-Easton-RussetBurbank-Umatilla.rds")
m1.1 <- readRDS("Models/model-1_Bohman_Minnesota_RussetBurbank.rds")
m1.2 <- readRDS("Models/model-1_Bohman_Minnesota_Clearwater.rds")
m1.3 <- readRDS("Models/model-1_Bohman_Minnesota_Umatilla.rds")
m1.4 <- readRDS("Models/model-1_Bohman_Minnesota_DakotaRusset.rds")
m1.5 <- readRDS("Models/model-1_Bohman_Minnesota_Easton.rds")

m2.0 <- readRDS("Models/model-2_Bohman_Minnesota_Clearwater-DakotaRusset-Easton-RussetBurbank-Umatilla.rds")

# m2.0.shiny <- as.shinystan(m2.0$fit)
# launch_shinystan(m2.0.shiny)

m3.0 <- readRDS("Models/model-3_Bohman_Minnesota_Clearwater-DakotaRusset-Easton-RussetBurbank-Umatilla.rds")

m3.0

# m3.0$fit
# m3.0$prior

# pairs(m3.0)
# m3.0.shiny <- as.shinystan(m3.0$fit)
# launch_shinystan(m3.0.shiny)

# the fmin() function used in Stan isn't defined in R, so we need to create it so that when we try to use brms to make predictions, it knows what to do with the fmin()
fmin <- function(x,y){
  pmin(x,y)
}

# now we can do cool stuff like take the orignal data, generate predicted outcomes, and plot those against the real values

# m3.0$data %>% 
#   add_predicted_draws(m3.0) %>% 
#   median_hdci() %>% 
#   rowid_to_column() %>% 
#   ggplot(aes(x = rowid, y = W)) +
#   geom_pointinterval(aes(ymin = .lower, ymax = .upper, y = .prediction), alpha = 0.5) +
#   geom_point(color = "forestgreen")
# 
# m3.0$data %>% 
#   add_predicted_draws(m3.0) %>% 
#   median_hdci() %>% 
#   ggplot(aes(x = W, y = N)) +
#   geom_pointintervalh(aes(xmin = .lower, xmax = .upper, x = .prediction), alpha = 0.5) +
#   geom_point(color = "forestgreen") +
#   scale_y_log10()

# or we can do a posterior check

# pp_check(m3.0, nsamples = 30)

# looks pretty good!!


# looking at group level draws --------------------------------------------

# Model 3.0

#get_variables(m3.0)

m3.0 %>% 
  spread_draws(b_alpha1_Intercept, r_variety__alpha1[variety,]) %>% 
  mutate(variety_alpha1 = b_alpha1_Intercept + r_variety__alpha1) %>% 
  ggplot(aes(x = variety_alpha1, y = variety)) +
  geom_halfeyeh()

m3.0 %>% 
  spread_draws(b_alpha2_Intercept, r_variety__alpha2[variety,]) %>% 
  mutate(variety_alpha2 = b_alpha2_Intercept + r_variety__alpha2) %>% 
  ggplot(aes(x = variety_alpha2, y = variety)) +
  geom_halfeyeh()

left_join(
  m3.0 %>% 
    spread_draws(b_alpha1_Intercept, r_variety__alpha1[variety,]) %>% 
    mutate(variety_alpha1 = b_alpha1_Intercept + r_variety__alpha1),
  m3.0 %>% 
    spread_draws(b_alpha2_Intercept, r_variety__alpha2[variety,]) %>% 
    mutate(variety_alpha2 = b_alpha2_Intercept + r_variety__alpha2),
  by = c(".chain", ".iteration", ".draw", "variety")) %>% 
  ggplot(aes(x = variety_alpha1, y = variety_alpha2, color=variety)) +
  geom_point(alpha=0.005) +
  geom_smooth(formula="y~x",method="lm") +
  theme_classic() +
  scale_color_brewer(palette = "Set1")

# this is how you would go about calculating the difference between alpha values by variety. instead of Bmax, it would be one of the alpha values, and instead of date, it would be variety. the key here is spread_draws to get the global mean and offset for each group, then mutate to make it into the actual group-level estimate, then compare_levels to get every pairwise difference.
m3.0 %>% 
  spread_draws(b_alpha1_Intercept, r_variety__alpha1[variety,]) %>% 
  mutate(variety_alpha1 = b_alpha1_Intercept + r_variety__alpha1) %>% 
  # filter(str_detect(date, "^19")) %>% # this is just for the example since there are too many dates
  compare_levels(variety_alpha1, by = variety) %>% 
  ggplot(aes(x = variety_alpha1, y = variety)) +
  geom_halfeyeh()

m3.0 %>% 
  spread_draws(b_alpha2_Intercept, r_variety__alpha2[variety,]) %>% 
  mutate(variety_alpha2 = b_alpha2_Intercept + r_variety__alpha2) %>% 
  # filter(str_detect(date, "^19")) %>% # this is just for the example since there are too many dates
  compare_levels(variety_alpha2, by = variety) %>% 
  ggplot(aes(x = variety_alpha2, y = variety)) +
  geom_halfeyeh()

# Model 1.0

get_variables(m1.0)

m1.0$prior

m1.0 %>% 
  spread_draws(b_Bmax_Intercept, r_date__Bmax[date,]) %>% 
  mutate(date_Bmax = b_Bmax_Intercept + r_date__Bmax) %>% 
  ggplot(aes(x = date_Bmax, y = date)) +
  geom_halfeyeh()

# this is how you would go about calculating the difference between alpha values by variety. instead of Bmax, it would be one of the alpha values, and instead of date, it would be variety. the key here is spread_draws to get the global mean and offset for each group, then mutate to make it into the actual group-level estimate, then compare_levels to get every pairwise difference.
m1.0 %>% 
  spread_draws(b_Bmax_Intercept, r_date__Bmax[date,]) %>% 
  mutate(date_Bmax = b_Bmax_Intercept + r_date__Bmax) %>% 
  filter(str_detect(date, "^19")) %>% # this is just for the example since there are too many dates
  compare_levels(date_Bmax, by = date) %>% 
  ggplot(aes(x = date_Bmax, y = date)) +
  geom_halfeyeh()



# actually comparing alpha1 and alpha2 across varieties -------------------


m3.0 %>% 
  spread_draws(b_alpha1_Intercept, r_variety__alpha1[variety,]) %>% 
  mutate(variety_alpha1 = b_alpha1_Intercept + r_variety__alpha1) %>% 
  compare_levels(variety_alpha1, by = variety) %>% 
  ggplot(aes(x = variety_alpha1, y = variety)) +
  geom_halfeyeh()

m3.0 %>% 
  spread_draws(b_alpha2_Intercept, r_variety__alpha2[variety,]) %>% 
  mutate(variety_alpha2 = b_alpha2_Intercept + r_variety__alpha2) %>% 
  compare_levels(variety_alpha2, by = variety) %>% 
  ggplot(aes(x = variety_alpha2, y = variety)) +
  geom_halfeyeh()

# model 2
m2.0 %>% 
  spread_draws(b_alpha1_Intercept, r_variety__alpha1[variety,]) %>% 
  mutate(variety_alpha1 = b_alpha1_Intercept + r_variety__alpha1) %>% 
  compare_levels(variety_alpha1, by = variety) %>% 
  ggplot(aes(x = variety_alpha1, y = variety)) +
  geom_halfeyeh()

m2.0 %>% 
  spread_draws(b_alpha2_Intercept, r_variety__alpha2[variety,]) %>% 
  mutate(variety_alpha2 = b_alpha2_Intercept + r_variety__alpha2) %>% 
  compare_levels(variety_alpha2, by = variety) %>% 
  ggplot(aes(x = variety_alpha2, y = variety)) +
  geom_halfeyeh()



