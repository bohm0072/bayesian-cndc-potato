# initialization -----------------

library(tidyverse)
library(brms)
library(tidybayes)

# read in data -------------------
data <- read_csv("data/analysis/data_cndc.csv",col_types="cccccccdcdd"); data_cndc = data
# data_cndc_index <- read_csv("data/analysis/data_cndc_index.csv",col_types="ccccccc"); #data_index = data_cndc_index

# read in model fit results ------------------

# model7 <- readRDS("brms/models/m0007.rds"); model7
# model8 <- readRDS("brms/models/m0008.rds"); model8
model9 <- readRDS("brms/models/m0009.rds"); model9
model10 <- readRDS("brms/models/m0010.rds"); model10

# the fmin() function used in Stan isn't defined in R, so we need to create it so that when we try to use brms to make predictions, it knows what to do with the fmin()
fmin <- function(x,y){
  pmin(x,y)
}

# evaluating model fit  -------------------------------------------

# model=model9
# model=model10

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
  
  p3 <- model %>%
    spread_draws(b_Bmax_Intercept, `r_index__Bmax`[`index`,]) %>%
    mutate(`index_Bmax` = b_Bmax_Intercept + `r_index__Bmax`) %>%
    mutate_at(vars(index),as.character) %>% 
    ggplot(aes(x = `index_Bmax`, y = `index`)) +
    stat_halfeye()
  
  p3.1 <- model %>%
    spread_draws(b_Bmax_Intercept, `r_index__Bmax`[`index`,]) %>%
    mutate(`index_Bmax` = b_Bmax_Intercept + `r_index__Bmax`) %>%
    mutate_at(vars(index),as.character) %>% 
    ggplot(aes(x = `index_Bmax`)) +
    stat_halfeye()
    
  p4 <- model %>%
    spread_draws(b_Si_Intercept, `r_index__Si`[`index`,]) %>%
    mutate(`index_Si` = b_Si_Intercept + `r_index__Si`) %>%
    mutate_at(vars(index),as.character) %>% 
    ggplot(aes(x = `index_Si`, y = `index`)) +
    stat_halfeye()
  
  p4.1 <- model %>%
    spread_draws(b_Si_Intercept, `r_index__Si`[`index`,]) %>%
    mutate(`index_Si` = b_Si_Intercept + `r_index__Si`) %>%
    mutate_at(vars(index),as.character) %>% 
    ggplot(aes(x = `index_Si`)) +
    stat_halfeye()
  
  p5 <- left_join(
    model %>%
      spread_draws(b_alpha1_Intercept, `r_group__alpha1`[`group`,]) %>%
      mutate(`group_alpha1` = b_alpha1_Intercept + `r_group__alpha1`),
    model %>%
      spread_draws(b_alpha2_Intercept, `r_group__alpha2`[`group`,]) %>%
      mutate(`group_alpha2` = b_alpha2_Intercept + `r_group__alpha2`),
    by = c(".chain", ".iteration", ".draw", "group")) %>%
    mutate_at(vars(group),as.character) %>% 
    ggplot(aes(x = `group_alpha1`, y = `group_alpha2`, color=`group`)) +
    geom_point(alpha=0.01) +
    geom_smooth(formula="y~x",method="lm") +
    theme_classic() +
    scale_color_brewer(palette = "Set3")

  # this is how you would go about calculating the difference between alpha values by variety.
  p6 <- model %>%
    spread_draws(b_alpha1_Intercept, `r_group__alpha1`[`group`,]) %>%
    mutate(`group_alpha1` = b_alpha1_Intercept + `r_group__alpha1`) %>%
    compare_levels(`group_alpha1`, by = `group`) %>%
    ggplot(aes(x = `group_alpha1`, y = `group`)) +
    stat_halfeye()

  p7 <- model %>%
    spread_draws(b_alpha2_Intercept, `r_group__alpha2`[`group`,]) %>%
    mutate(`group_alpha2` = b_alpha2_Intercept + `r_group__alpha2`) %>%
    compare_levels(`group_alpha2`, by = `group`) %>%
    ggplot(aes(x = `group_alpha2`, y = `group`)) +
    stat_halfeye()
  
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
    mutate(W=fmin(Bmax + Si * (N - Nc), Bmax)) %>%
    left_join(
      data %>%
        select(index,group) %>%
        distinct,
      by="index"
      ) %>%
    relocate(group,.before=index)
  
  p8 <- eval3 %>%
    ggplot() +
    geom_line(aes(x=W,y=N,group=index),alpha=1.0) +
    theme_classic() +
    facet_wrap(vars(as.numeric(index))) +
    geom_point(data=data_cndc,aes(x=W,y=N,group=index),alpha=0.25)
  
  p9 <- eval3 %>%
    ggplot() +
    geom_line(aes(x=W,y=N,group=index),alpha=0.5) +
    theme_classic() +
    facet_wrap(vars(as.numeric(group))) +
    scale_x_continuous(limits=c(0,NA))
  
  out <- list(p1,
              p2,
              p3,
              p3.1,
              p4,
              p4.1,
              p5,
              p6,
              p7,
              p8,
              p9)
  
  return(out)
    
  
}

# eval8 <- f.eval(model8,data_cndc)
# eval9 <- f.eval(model9,data_cndc)

f.eval2 <- function(model,data){
  
  p1.1 <- model %>%
    spread_draws(b_alpha1_Intercept, `r_location__alpha1`[`location`,], `r_location:variety__alpha1`[`location:variety`,]) %>%
    rowwise() %>%
    filter(is.na(str_match(`location:variety`,location))==F) %>%
    ungroup() %>%
    mutate(`group_alpha1` = b_alpha1_Intercept + `r_location__alpha1` + `r_location:variety__alpha1`) %>%
    mutate_at(vars(`location:variety`),as.character) %>%
    ggplot(aes(x = `group_alpha1`, y = `location:variety`)) +
    stat_halfeye()
  
  p1.2 <- model %>%
    spread_draws(b_alpha1_Intercept, `r_location__alpha1`[`location`,], `r_location:variety__alpha1`[`location:variety`,]) %>%
    rowwise() %>%
    filter(is.na(str_match(`location:variety`,location))==F) %>%
    ungroup() %>%
    mutate(`group_alpha1` = b_alpha1_Intercept + `r_location__alpha1`) %>%
    mutate_at(vars(location),as.character) %>%
    ggplot(aes(x = `group_alpha1`, y = `location`)) +
    stat_halfeye()
  
  p2.1 <- model %>%
    spread_draws(b_alpha2_Intercept, `r_location__alpha2`[`location`,], `r_location:variety__alpha2`[`location:variety`,]) %>%
    rowwise() %>%
    filter(is.na(str_match(`location:variety`,location))==F) %>%
    ungroup() %>%
    mutate(`group_alpha2` = b_alpha2_Intercept + `r_location__alpha2` + `r_location:variety__alpha2`) %>%
    mutate_at(vars(`location:variety`),as.character) %>%
    ggplot(aes(x = `group_alpha2`, y = `location:variety`)) +
    stat_halfeye()
  
  p2.2 <- model %>%
    spread_draws(b_alpha2_Intercept, `r_location__alpha2`[`location`,], `r_location:variety__alpha2`[`location:variety`,]) %>%
    rowwise() %>%
    filter(is.na(str_match(`location:variety`,location))==F) %>%
    ungroup() %>%
    mutate(`group_alpha2` = b_alpha2_Intercept + `r_location__alpha2`) %>%
    mutate_at(vars(location),as.character) %>%
    ggplot(aes(x = `group_alpha2`, y = `location`)) +
    stat_halfeye()
  
  p3.1 <- model %>%
    spread_draws(b_Bmax_Intercept, `r_index__Bmax`[`index`,]) %>%
    mutate(`index_Bmax` = b_Bmax_Intercept + `r_index__Bmax`) %>%
    mutate_at(vars(index),as.character) %>% 
    ggplot(aes(x = `index_Bmax`, y = `index`)) +
    stat_halfeye()
  
  p3.2 <- model %>%
    spread_draws(b_Bmax_Intercept, `r_index__Bmax`[`index`,]) %>%
    mutate(`index_Bmax` = b_Bmax_Intercept + `r_index__Bmax`) %>%
    mutate_at(vars(index),as.character) %>% 
    ggplot(aes(x = `index_Bmax`)) +
    stat_halfeye()
  
  p4.1 <- model %>%
    spread_draws(b_Si_Intercept, `r_index__Si`[`index`,]) %>%
    mutate(`index_Si` = b_Si_Intercept + `r_index__Si`) %>%
    mutate_at(vars(index),as.character) %>% 
    ggplot(aes(x = `index_Si`, y = `index`)) +
    stat_halfeye()
  
  p4.2 <- model %>%
    spread_draws(b_Si_Intercept, `r_index__Si`[`index`,]) %>%
    mutate(`index_Si` = b_Si_Intercept + `r_index__Si`) %>%
    mutate_at(vars(index),as.character) %>% 
    ggplot(aes(x = `index_Si`)) +
    stat_halfeye()
  
  p5 <- left_join(
    model %>%
      spread_draws(b_alpha1_Intercept, `r_location__alpha1`[`location`,], `r_location:variety__alpha1`[`location:variety`,]) %>%
      rowwise() %>%
      filter(is.na(str_match(`location:variety`,location))==F) %>%
      ungroup() %>%
      mutate(`group_alpha1` = b_alpha1_Intercept + `r_location__alpha1`),
    model %>%
      spread_draws(b_alpha2_Intercept, `r_location__alpha2`[`location`,], `r_location:variety__alpha2`[`location:variety`,]) %>%
      rowwise() %>%
      filter(is.na(str_match(`location:variety`,location))==F) %>%
      ungroup() %>%
      mutate(`group_alpha2` = b_alpha2_Intercept + `r_location__alpha2`),
    by = c(".chain", ".iteration", ".draw", "location", "location:variety")) %>%
    mutate_at(vars(location,`location:variety`),as.character) %>% 
    ggplot(aes(x = `group_alpha1`, y = `group_alpha2`, color=`group`)) +
    geom_point(alpha=0.01) +
    geom_smooth(formula="y~x",method="lm") +
    theme_classic() +
    scale_color_brewer(palette = "Set3")
  
  # this is how you would go about calculating the difference between alpha values by variety.
  p6 <- model %>%
    spread_draws(b_alpha1_Intercept, `r_group__alpha1`[`group`,]) %>%
    mutate(`group_alpha1` = b_alpha1_Intercept + `r_group__alpha1`) %>%
    compare_levels(`group_alpha1`, by = `group`) %>%
    ggplot(aes(x = `group_alpha1`, y = `group`)) +
    stat_halfeye()
  
  p7 <- model %>%
    spread_draws(b_alpha2_Intercept, `r_group__alpha2`[`group`,]) %>%
    mutate(`group_alpha2` = b_alpha2_Intercept + `r_group__alpha2`) %>%
    compare_levels(`group_alpha2`, by = `group`) %>%
    ggplot(aes(x = `group_alpha2`, y = `group`)) +
    stat_halfeye()
  
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
    mutate(W=fmin(Bmax + Si * (N - Nc), Bmax)) %>%
    left_join(
      data %>%
        select(index,group) %>%
        distinct,
      by="index"
    ) %>%
    relocate(group,.before=index)
  
  p8 <- eval3 %>%
    ggplot() +
    geom_line(aes(x=W,y=N,group=index),alpha=1.0) +
    theme_classic() +
    facet_wrap(vars(as.numeric(index))) +
    geom_point(data=data_cndc,aes(x=W,y=N,group=index),alpha=0.25)
  
  p9 <- eval3 %>%
    ggplot() +
    geom_line(aes(x=W,y=N,group=index),alpha=0.5) +
    theme_classic() +
    facet_wrap(vars(as.numeric(group))) +
    scale_x_continuous(limits=c(0,NA))
  
  out <- list(p1,
              p2,
              p3,
              p3.1,
              p4,
              p4.1,
              p5,
              p6,
              p7,
              p8,
              p9)
  
  return(out)
  
  
}

eval10 <- f.eval2(model9,data_cndc)

# trying to compare curves ------------------------------------------------

#Nc=alpha1*(W^(-alpha2))

# function to calculate curves for a range of W values
calc_nc <- function(W = seq(from = 1, to = 10, length.out = 100), alpha1, alpha2){
  Nc <- alpha1*(W^(-alpha2))
  d <- tibble(Nc = Nc, W = W)
  return(d)
}

# get alpha values for all group (variety x location)
a <- model %>% 
  spread_draws(b_alpha1_Intercept, r_group__alpha1[group,], n = 100, seed = 1) %>% 
  mutate(group_alpha1 = b_alpha1_Intercept + r_group__alpha1) %>% 
  left_join(
    model %>% 
      spread_draws(b_alpha2_Intercept, r_group__alpha2[group,], n = 100, seed = 1) %>% 
      mutate(group_alpha2 = b_alpha2_Intercept + r_group__alpha2) 
  ) %>% 
  select(.draw, group, group_alpha1, group_alpha2)

# generate Nc estimates for every variety/alpha combination across all W values
a <- a %>% 
  rowwise() %>% 
  mutate(values = list(calc_nc(W = seq(from = 1, to = 30, length.out = 3000), alpha1 = group_alpha1, alpha2 = group_alpha2))) %>% 
  unnest(values)

# plot one curve for each draw+variety combo
a %>% 
  ggplot(aes(x = W, y = Nc, color = as.character(group), group = interaction(.draw, group))) +
  geom_line(alpha = 0.2) +
  scale_color_viridis_d() +
  theme_minimal()

# same as above, but median and 95% credible intervals
a %>% 
  group_by(group, W) %>% 
  median_hdci(Nc) %>% 
  ggplot(aes(x = W, y = Nc, color = as.character(group), fill = as.character(group))) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.2, color = NA) +
  geom_line() +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  theme_minimal()

# actually calculating the difference between two varieties and plotting that
# Group 1: MN Russet Burbank
# Group 6: Canada Russet Burbank

a %>% 
  select(.draw, group, Nc, W) %>% 
  filter(group %in% c("1", "6")) %>% 
  pivot_wider(names_from = group, values_from = Nc) %>% 
  mutate(diff = `1` - `6`) %>% 
  select(.draw, W, diff) %>% 
  group_by(W) %>% 
  median_hdci(diff) %>% 
  ggplot(aes(x = W, y = diff)) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.2, color = NA) +
  geom_line() +
  theme_minimal() +
  labs(y = "Diff in Nc between MN Russet Burbank and Canada Russet Burbank")


# older stuff -------------------------------------

# actually comparing alpha1 and alpha2 across varieties

# m3.0 %>% 
#   spread_draws(b_alpha1_Intercept, r_variety__alpha1[variety,]) %>% 
#   mutate(variety_alpha1 = b_alpha1_Intercept + r_variety__alpha1) %>% 
#   compare_levels(variety_alpha1, by = variety) %>% 
#   ggplot(aes(x = variety_alpha1, y = variety)) +
#   geom_halfeyeh()
# 
# m3.0 %>% 
#   spread_draws(b_alpha2_Intercept, r_variety__alpha2[variety,]) %>% 
#   mutate(variety_alpha2 = b_alpha2_Intercept + r_variety__alpha2) %>% 
#   compare_levels(variety_alpha2, by = variety) %>% 
#   ggplot(aes(x = variety_alpha2, y = variety)) +
#   geom_halfeyeh()
# 
# # model 2
# m2.0 %>% 
#   spread_draws(b_alpha1_Intercept, r_variety__alpha1[variety,]) %>% 
#   mutate(variety_alpha1 = b_alpha1_Intercept + r_variety__alpha1) %>% 
#   compare_levels(variety_alpha1, by = variety) %>% 
#   ggplot(aes(x = variety_alpha1, y = variety)) +
#   geom_halfeyeh()
# 
# m2.0 %>% 
#   spread_draws(b_alpha2_Intercept, r_variety__alpha2[variety,]) %>% 
#   mutate(variety_alpha2 = b_alpha2_Intercept + r_variety__alpha2) %>% 
#   compare_levels(variety_alpha2, by = variety) %>% 
#   ggplot(aes(x = variety_alpha2, y = variety)) +
#   geom_halfeyeh()

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


# Model 1.0

# get_variables(m1.0)
# 
# m1.0$prior
# 
# m1.0 %>% 
#   spread_draws(b_Bmax_Intercept, r_date__Bmax[date,]) %>% 
#   mutate(date_Bmax = b_Bmax_Intercept + r_date__Bmax) %>% 
#   ggplot(aes(x = date_Bmax, y = date)) +
#   geom_halfeyeh()
# 
# # this is how you would go about calculating the difference between alpha values by variety. instead of Bmax, it would be one of the alpha values, and instead of date, it would be variety. the key here is spread_draws to get the global mean and offset for each group, then mutate to make it into the actual group-level estimate, then compare_levels to get every pairwise difference.
# m1.0 %>% 
#   spread_draws(b_Bmax_Intercept, r_date__Bmax[date,]) %>% 
#   mutate(date_Bmax = b_Bmax_Intercept + r_date__Bmax) %>% 
#   filter(str_detect(date, "^19")) %>% # this is just for the example since there are too many dates
#   compare_levels(date_Bmax, by = date) %>% 
#   ggplot(aes(x = date_Bmax, y = date)) +
#   geom_halfeyeh()



##### END #####