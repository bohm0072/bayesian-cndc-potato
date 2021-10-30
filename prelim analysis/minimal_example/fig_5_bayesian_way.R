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

alphas <- model %>% 
  spread_draws(b_alpha1_Intercept, 
               `r_location:variety__alpha1`[`location:variety`,],
               r_location__alpha1[location,], n = 100, seed = 1) %>% 
  mutate(lv_alpha1 = b_alpha1_Intercept + `r_location:variety__alpha1` + r_location__alpha1) %>% 
  left_join(
    model %>% 
      spread_draws(b_alpha2_Intercept, 
                   `r_location:variety__alpha2`[`location:variety`,],
                   r_location__alpha2[location,], n = 100, seed = 1) %>% 
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


# define reference and compare --------------------------------------------

reference <- alphas %>% 
  filter(lv == "Minnesota_Russet.Burbank") %>% 
  rename(ref_Nc = Nc)

ref_summary <- reference %>% 
  select(.draw, ref_Nc, W) %>% 
  right_join(alphas, by = c(".draw", "W")) %>% 
  mutate(Nc_diff = Nc - ref_Nc) %>% 
  group_by(location, lv, W) %>% 
  median_qi(Nc_diff, .width = 0.9) %>% 
  mutate(crosses_0 = .lower < 0 & .upper > 0)


# plot --------------------------------------------------------------------

ref_summary %>% 
  ggplot(aes(x = W, y = Nc_diff, color = crosses_0)) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.2, color = NA) +
  geom_line(aes(group = lv), size= 1.2) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.2) +
  facet_wrap(vars(lv)) +
  theme_bw() + 
  scale_color_manual(values = c("red", "blue"))
  
