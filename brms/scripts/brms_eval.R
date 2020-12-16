library(tidyverse)
library(brms)
library(tidybayes)
# library(shinystan)

# read in data -------------------
data_cndc <- read_csv("data/analysis/data_cndc.csv",col_types="cccccccdcdd"); data = data_cndc
data_cndc_index <- read_csv("data/analysis/data_cndc_index.csv",col_types="ccccccc"); data_index = data_cndc_index

# read in model fit results ------------------

# readRDS("brms/models/m0000_Bohman_Minnesota_RussetBurbank.rds")
# readRDS("brms/models/m0001_Bohman_Minnesota_RussetBurbank.rds")
# readRDS("brms/models/m0002_Bohman_Minnesota_All.rds")
# readRDS("brms/models/m0003_Giletto_Canada_All.rds")
# readRDS("brms/models/m0004_Giletto_Argentina_All.rds")
# readRDS("brms/models/m0005_All.rds")
m0006 <- readRDS("brms/models/m0006_All.rds"); m0006; model = m0006

# pairs(m0006)

# looking at group level draws --------------------------------------------

# the fmin() function used in Stan isn't defined in R, so we need to create it so that when we try to use brms to make predictions, it knows what to do with the fmin()
fmin <- function(x,y){
  pmin(x,y)
}

f.eval <- function(model,data){
  
  p1 <- model %>%
    spread_draws(b_alpha1_Intercept, `r_group__alpha1`[`group`,]) %>%
    mutate(`group_alpha1` = b_alpha1_Intercept + `r_group__alpha1`) %>%
    mutate_at(vars(group),as.character) %>% 
    ggplot(aes(x = `group_alpha1`, y = `group`)) +
    stat_halfeye()

  p2 <- model %>%
    spread_draws(b_alpha2_Intercept, `r_group__alpha2`[`group`,]) %>%
    mutate(`group_alpha2` = b_alpha2_Intercept + `r_group__alpha2`) %>%
    mutate_at(vars(group),as.character) %>% 
    ggplot(aes(x = `group_alpha2`, y = `group`)) +
    stat_halfeye()

  p3 <- left_join(
    model %>%
      spread_draws(b_alpha1_Intercept, `r_group__alpha1`[`group`,]) %>%
      mutate(`group_alpha1` = b_alpha1_Intercept + `r_group__alpha1`),
    model %>%
      spread_draws(b_alpha2_Intercept, `r_group__alpha2`[`group`,]) %>%
      mutate(`group_alpha2` = b_alpha2_Intercept + `r_group__alpha2`),
    by = c(".chain", ".iteration", ".draw", "group")) %>%
    mutate_at(vars(group),as.character) %>% 
    ggplot(aes(x = `group_alpha1`, y = `group_alpha2`, color=`group`)) +
    geom_point(alpha=0.005) +
    geom_smooth(formula="y~x",method="lm") +
    theme_classic() +
    scale_color_brewer(palette = "Set1")

  # this is how you would go about calculating the difference between alpha values by variety.
  p4 <- model %>%
    spread_draws(b_alpha1_Intercept, `r_group__alpha1`[`group`,]) %>%
    mutate(`group_alpha1` = b_alpha1_Intercept + `r_group__alpha1`) %>%
    compare_levels(`group_alpha1`, by = `group`) %>%
    ggplot(aes(x = `group_alpha1`, y = `group`)) +
    stat_halfeye()

  p5 <- model %>%
    spread_draws(b_alpha2_Intercept, `r_group__alpha2`[`group`,]) %>%
    mutate(`group_alpha2` = b_alpha2_Intercept + `r_group__alpha2`) %>%
    compare_levels(`group_alpha2`, by = `group`) %>%
    ggplot(aes(x = `group_alpha2`, y = `group`)) +
    stat_halfeye()
  
  out <- list(p1,p2,p3,p4,p5)
  return(out)
  
  eval1 <- data %>%
    left_join(
      left_join(
        model %>%
          spread_draws(b_Bmax_Intercept, r_index__Bmax[index,]) %>%
          mutate(index_Bmax = b_Bmax_Intercept + r_index__Bmax),
        model %>%
          spread_draws(b_Si_Intercept, r_index__Si[index,]) %>%
          mutate(index_Si = b_Si_Intercept + r_index__Si),
        by = c(".chain", ".iteration", ".draw", "index")) %>%
        select(.chain,.iteration,.draw,index,index_Bmax,index_Si) %>%
        rename(Bmax=index_Bmax,Si=index_Si) %>%
        mutate_at(vars(index),as.character),
      by="index") %>% 
    left_join(
      left_join(
        model %>%
          spread_draws(b_alpha1_Intercept, r_group__alpha1[group,]) %>%
          mutate(group_alpha1 = b_alpha1_Intercept + r_group__alpha1),
        model %>%
          spread_draws(b_alpha2_Intercept, r_group__alpha2[group,]) %>%
          mutate(group_alpha2 = b_alpha2_Intercept + r_group__alpha2),
        by = c(".chain", ".iteration", ".draw", "group")) %>%
        select(.chain,.iteration,.draw,group,group_alpha1,group_alpha2) %>%
        rename(alpha1=group_alpha1,alpha2=group_alpha2) %>%
        mutate_at(vars(group),as.character),
      by=c(".chain",".iteration",".draw","group")) %>%
    mutate(Nc=alpha1*(Bmax^(-alpha2)))
  
  eval2 <- eval1 %>%
    group_by(index) %>%
    summarise_at(vars(Bmax,Si,Nc),mean) %>%
    arrange(as.numeric(index))
  
  eval3 <- expand.grid(
    index=eval2 %>% pull(index),
    N=seq(0,7,0.1),stringsAsFactors=F
  ) %>%
    as_tibble() %>%
    arrange(as.numeric(index)) %>%
    left_join(
      eval2,
      by="index"
    ) %>%
    mutate(W=fmin(Bmax + Si * (N - Nc), Bmax)) 
  
  p6 <- eval3 %>%
    ggplot() +
    geom_line(aes(x=W,y=N,group=index),alpha=0.25) +
    theme_classic() +
    # scale_x_continuous(limits=c(0,50)) +
    # scale_y_continuous(limits=c(0,6)) +
    facet_wrap(vars(as.numeric(index))); p6
  
  p6 +
    geom_point(data=data_cndc,aes(x=W,y=N,group=index))
    
  
}

f.eval(m0006)

# older stuff -------------------------------------

# m2.0.shiny <- as.shinystan(m2.0$fit)
# launch_shinystan(m2.0.shiny)

# m3.0$fit
# m3.0$prior

# pairs(m3.0)
# m3.0.shiny <- as.shinystan(m3.0$fit)
# launch_shinystan(m3.0.shiny)

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



# trying to compare curves ------------------------------------------------

#Nc=alpha1*(W^(-alpha2))

# function to calculate curves for a range of W values
calc_nc <- function(W = seq(from = 1, to = 10, length.out = 100), alpha1, alpha2){
  Nc <- alpha1*(W^(-alpha2))
  d <- tibble(Nc = Nc, W = W)
  return(d)
}

# get alpha values for all varieties
a <- m3.0 %>% 
  spread_draws(b_alpha1_Intercept, r_variety__alpha1[variety,], n = 100, seed = 1) %>% 
  mutate(variety_alpha1 = b_alpha1_Intercept + r_variety__alpha1) %>% 
  left_join(
    m3.0 %>% 
      spread_draws(b_alpha2_Intercept, r_variety__alpha2[variety,], n = 100, seed = 1) %>% 
      mutate(variety_alpha2 = b_alpha2_Intercept + r_variety__alpha2) 
  ) %>% 
  select(.draw, variety, variety_alpha1, variety_alpha2)

# generate Nc estimates for every variety/alpha combination across all W values
a <- a %>% 
  rowwise() %>% 
  mutate(values = list(calc_nc(W = seq(from = 1, to = 10, length.out = 1000), alpha1 = variety_alpha1, alpha2 = variety_alpha2))) %>% 
  unnest(values)

# plot one curve for each draw+variety combo
a %>% 
  ggplot(aes(x = W, y = Nc, color = variety, group = interaction(.draw, variety))) +
  geom_line(alpha = 0.2) +
  scale_color_viridis_d() +
  theme_minimal()

# same as above, but median and 95% credible intervals
a %>% 
  group_by(variety, W) %>% 
  median_hdci(Nc) %>% 
  ggplot(aes(x = W, y = Nc, color = variety, fill = variety)) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.2, color = NA) +
  geom_line() +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  theme_minimal()

# actually calculating the difference between two varieties and plotting that
a %>% 
  select(.draw, variety, Nc, W) %>% 
  filter(variety %in% c("Russet.Burbank", "Easton")) %>% 
  pivot_wider(names_from = variety, values_from = Nc) %>% 
  mutate(diff = Russet.Burbank - Easton) %>% 
  select(.draw, W, diff) %>% 
  group_by(W) %>% 
  median_hdci(diff) %>% 
  ggplot(aes(x = W, y = diff)) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.2, color = NA) +
  geom_line() +
  theme_minimal() +
  labs(y = "Diff in Nc between Russet Burbank and Easton")


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



