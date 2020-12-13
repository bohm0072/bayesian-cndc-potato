library(tidyverse)
library(brms)
library(tidybayes)

m0 <- readRDS("Models/model-1_Bohman_Minnesota_Clearwater-DakotaRusset-Easton-RussetBurbank-Umatilla.rds")
m1 <- readRDS("Models/model-1_Bohman_Minnesota_RussetBurbank.rds")
m2 <- readRDS("Models/model-1_Bohman_Minnesota_Clearwater.rds")
m3 <- readRDS("Models/model-1_Bohman_Minnesota_Umatilla.rds")
m4 <- readRDS("Models/model-1_Bohman_Minnesota_DakotaRusset.rds")
m5 <- readRDS("Models/model-1_Bohman_Minnesota_Easton.rds")

m1$fit


m1

# the fmin() function used in Stan isn't defined in R, so we need to create it so that when we try to use brms to make predictions, it knows what to do with the fmin()
fmin <- function(x,y){
  pmin(x,y)
}

# now we can do cool stuff like take the orignal data, generate predicted outcomes, and plot those against the real values
m1$data %>% 
  add_predicted_draws(m1) %>% 
  median_hdci() %>% 
  ggplot(aes(x = N, y = W)) +
  geom_pointinterval(aes(ymin = .lower, ymax = .upper, y = .prediction), alpha = 0.5) +
  geom_point(color = "forestgreen")

# or we can do a posterior check

pp_check(m1, nsamples = 30)

# looks pretty good!!
