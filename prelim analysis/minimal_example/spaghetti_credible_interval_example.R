library(brms)
library(tidyverse)
library(tidybayes)

model <- read_rds("brms/models/model_060221.rds")


# nonlinear formula functions ----------------------------------------------------------

fmin <- function(x,y){
  pmin(x,y)
}

calc_nc <- function(W = seq(from = 1, to = 10, length.out = 100), alpha1, alpha2){
  Nc <- alpha1*(W^(-alpha2))
  d <- tibble(Nc = Nc, W = W)
  return(d)
}


# spread draws ------------------------------------------------------------

# NOTE!!!!!! Here we take 500 draws because we want a smoother interval
alphas <- model %>% 
  spread_draws(b_alpha1_Intercept, 
               `r_location:variety__alpha1`[`location:variety`,],
               r_location__alpha1[location,], n = 500, seed = 1) %>% 
  mutate(lv_alpha1 = b_alpha1_Intercept + `r_location:variety__alpha1` + r_location__alpha1) %>% 
  left_join(
    model %>% 
      spread_draws(b_alpha2_Intercept, 
                   `r_location:variety__alpha2`[`location:variety`,],
                   r_location__alpha2[location,], n = 500, seed = 1) %>% 
      mutate(lv_alpha2 = b_alpha2_Intercept + `r_location:variety__alpha2` + r_location__alpha2) 
  ) %>% 
  select(.draw, lv = `location:variety`, lv_alpha1, lv_alpha2) %>% 
  ungroup() %>% 
  filter(is.na(str_match(lv,location))==F)


# estimate Nc across W ----------------------------------------------------

alphas <- alphas %>% 
  rowwise() %>% 
  mutate(values = list(calc_nc(W = seq(from = 1, to = 30, length.out = 3000), alpha1 = lv_alpha1, alpha2 = lv_alpha2))) %>% 
  unnest(values) %>% 
  ungroup()


# pick one variety x location --------------------------------------------

reference <- alphas %>% 
  filter(lv == "Minnesota_Russet.Burbank")

# I'm centering and scaling Nc here, which should be similar to the % difference used in the manuscript figures. Not exactly sure how that's being calculated, but I think it's basically the same thing, saying how different a Nc value is from the median Nc value at that level of W
reference <- reference %>% 
  group_by(W) %>% 
  mutate(Nc_scale = scale(Nc) %>% as.vector()) 

ref_med_qi <- reference %>%
  group_by(W) %>% 
  median_qi(Nc_scale, .width = 0.9)


# plot --------------------------------------------------------------------

# first we select 100 draws, then we plot them as a spaghetti plot, but we also include the median + 90% interval for the full 500 draws

reference %>% 
  group_by(.draw) %>% 
  nest() %>% 
  ungroup() %>% 
  slice_sample(n = 100) %>% 
  unnest() %>% 
  ggplot(aes(x = W, y = Nc_scale)) +
  geom_line(alpha = 0.2, aes(group = .draw)) +
  geom_line(data = ref_med_qi, color = "forestgreen") +
  geom_ribbon(data = ref_med_qi, aes(ymin = .lower, ymax = .upper), alpha = 0.1, color = NA, fill = "forestgreen") +
  theme_minimal()

