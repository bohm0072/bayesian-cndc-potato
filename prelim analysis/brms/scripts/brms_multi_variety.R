library(tidyverse)
library(brms)
library(tidybayes)


# data prep ---------------------------------------------------------------

data <- read_csv("data.csv",col_types="cccccdcdd") 

data

f.data.filter <- function(owner,location,variety,W_min,date_min,sd_min){
  
  d <- data %>% 
    filter(owner%in%owner) %>% #"Bohman") %>% #=="Bohman") %>% #
    filter(location%in%location) %>% #"Minnesota") %>% #=="Minnesota") %>% #
    filter(variety%in%variety) #"Russet Burbank") #=="Russet Burbank") #
  
  d <- d %>%
    filter(W>=1) %>%
    drop_na()
  
  filter_date <- d %>%
    group_by(date) %>%
    count() %>%
    filter(n>=3) %>%
    select(date) %>%
    pull()
  
  d <- d %>%
    filter(date%in%filter_date)
  
  d.sum <- d %>% 
    group_by(date) %>%
    summarize_at(vars(W),lst(mean,sd)) %>%
    rowwise() %>%
    mutate(cv=sd/mean) %>%
    ungroup()
  
  filter_sd <- d.sum %>%
    filter(sd>=1.0) %>%
    select(date) %>%
    pull()
  
  d <- d %>%
    filter(date%in%filter_sd)
  
  return(d)
  
}

d <- f.data.filter(owner="Bohman",
                   location="Minnesota",
                   W_min=1.0,  # Required to conform with convention that no critical points are defined for W <= 1.0
                   date_min=3, # Required to properly fit linear-plateu model
                   sd_min=1.0) # Required to properly fit linear-plateau model

d

d.plot <- d %>%
  ggplot(aes(x = W, y = N, color = variety)) +
  geom_point() +
  facet_wrap(vars(date))

d.plot

# model specification -----------------------------------------------------


formula_2 <- bf(W ~ fmin(Bmax + Si * (N - (alpha1*(Bmax^(-alpha2)))), Bmax),
                Bmax + Si ~ 1 + (1|variety/date),
                alpha1 + alpha2 ~ 1 + (1|variety),
                nl = T)

get_prior(formula_2, family = gaussian, data = d)

# prior class      coef        group resp dpar  nlpar bound       source
# student_t(3, 0, 5.3) sigma                                                    default
# (flat)     b                                  alpha1            default
# (flat)     b Intercept                        alpha1       (vectorized)
# student_t(3, 0, 5.3)    sd                                  alpha1            default
# student_t(3, 0, 5.3)    sd                variety           alpha1       (vectorized)
# student_t(3, 0, 5.3)    sd Intercept      variety           alpha1       (vectorized)
# (flat)     b                                  alpha2            default
# (flat)     b Intercept                        alpha2       (vectorized)
# student_t(3, 0, 5.3)    sd                                  alpha2            default
# student_t(3, 0, 5.3)    sd                variety           alpha2       (vectorized)
# student_t(3, 0, 5.3)    sd Intercept      variety           alpha2       (vectorized)
# (flat)     b                                    Bmax            default
# (flat)     b Intercept                          Bmax       (vectorized)
# student_t(3, 0, 5.3)    sd                                    Bmax            default
# student_t(3, 0, 5.3)    sd                variety             Bmax       (vectorized)
# student_t(3, 0, 5.3)    sd Intercept      variety             Bmax       (vectorized)
# student_t(3, 0, 5.3)    sd           variety:date             Bmax       (vectorized)
# student_t(3, 0, 5.3)    sd Intercept variety:date             Bmax       (vectorized)
# (flat)     b                                      Si            default
# (flat)     b Intercept                            Si       (vectorized)
# student_t(3, 0, 5.3)    sd                                      Si            default
# student_t(3, 0, 5.3)    sd                variety               Si       (vectorized)
# student_t(3, 0, 5.3)    sd Intercept      variety               Si       (vectorized)
# student_t(3, 0, 5.3)    sd           variety:date               Si       (vectorized)
# student_t(3, 0, 5.3)    sd Intercept variety:date               Si       (vectorized)

priors <- c(set_prior("normal(5,1)", nlpar = "alpha1", lb = 0, ub = 10),
            set_prior("normal(0.5,0.1)", nlpar = "alpha2", lb = 0, ub = 1),
            set_prior("normal(10,10)", nlpar = "Bmax", lb = 1, ub = 30),
            set_prior("normal(6,1)", nlpar = "Si", lb = 0),
            set_prior("cauchy(0,0.1)", class = "sd", nlpar = "Bmax"),
            set_prior("cauchy(0,0.1)", class = "sd", nlpar = "Si"),
            set_prior("cauchy(0,0.1)", class = "sd", nlpar = "alpha1"),
            set_prior("cauchy(0,0.1)", class = "sd", nlpar = "alpha2")
            )

get_prior(formula_2, family = gaussian, data = d, prior = priors)


# model fitting -----------------------------------------------------------

m2 <- brm(formula_2, data = d, family = gaussian, prior = priors,
          cores = 4, chains = 4, iter = 500, warmup = 200,
          control = list(adapt_delta = 0.99, max_treedepth = 15),
          sample_prior = "only")

m2

fmin <- function(x,y){
  pmin(x,y)
}

d %>% 
  add_predicted_draws(m2, n = 100) %>% 
  ggplot(aes(x = .prediction)) +
  geom_histogram() +
  scale_x_log10()

# Estimating group-specific effects ---------------------------------------

# this example is getting the Bmax estimate for each date. this isn't of interest, but the same pattern will apply when you start trying to do alpha values for each variety.

get_variables(m2)

m1 %>% 
  spread_draws(b_Bmax_Intercept, r_date__Bmax[date,]) %>% 
  mutate(date_Bmax = b_Bmax_Intercept + r_date__Bmax) %>% 
  ggplot(aes(x = date_Bmax, y = date)) +
  geom_halfeyeh()


# trying min instead of fmin ----------------------------------------------


formula_3 <- bf(W ~ min(Bmax + Si * (N - (alpha1*(Bmax^(-alpha2)))), Bmax),
                Bmax + Si ~ 1 + (1|variety/date),
                alpha1 + alpha2 ~ 1 + (1|variety),
                nl = T)

m3 <- brm(formula_3, data = d, family = gaussian, prior = priors,
          cores = 4, chains = 4, iter = 500, warmup = 200,
          control = list(adapt_delta = 0.99, max_treedepth = 15),
          sample_prior = "only")

