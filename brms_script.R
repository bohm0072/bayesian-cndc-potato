library(tidyverse)
library(brms)
library(pushoverr)

data <- read_csv("data.csv",col_types="cccccdcdd") 

data

f.data.filter <- function(Owner,Location,Variety,W_min,date_min,sd_min){
  
  d <- data %>% 
    filter(owner%in%Owner) %>% #"Bohman") %>% #=="Bohman") %>% #
    filter(location%in%Location) %>% #"Minnesota") %>% #=="Minnesota") %>% #
    filter(variety%in%Variety) #"Russet Burbank") #=="Russet Burbank") #
  
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

d <- f.data.filter(Owner="Bohman",
                   Location="Minnesota",
                   Variety=c("Russet Burbank","Clearwater","Umatilla","Dakota Russet","Easton"),
                   W_min=1.0,  # Required to conform with convention that no critical points are defined for W <= 1.0
                   date_min=3, # Required to properly fit linear-plateu model
                   sd_min=1.0) # Required to properly fit linear-plateau model

d

d.plot <- d %>%
  ggplot(aes(x = W, y = N)) +
  geom_point() +
  facet_wrap(vars(date))

formula_1 <- bf(W ~ fmin(Bmax + Si * (N - (alpha1*(Bmax^(-alpha2)))), Bmax),
                Bmax + Si ~ 1 + (1|date), #(1|Date/year)
                alpha1 + alpha2 ~ 1,
                nl = T)

get_prior(formula_1, family = gaussian, data = d)

priors_1 <- c(set_prior("normal(5,1)", nlpar = "alpha1", lb = 0, ub = 10),
            set_prior("normal(0.5,0.1)", nlpar = "alpha2", lb = 0, ub = 1),
            set_prior("normal(10,10)", nlpar = "Bmax", lb = 1, ub = 30),
            set_prior("normal(6,1)", nlpar = "Si", lb = 0))

m1 <- brm(formula_1, data = d, family = gaussian, prior = priors_1,
          cores = 4, chains = 4, iter = 5000, warmup = 2000,
          control = list(adapt_delta = 0.99, max_treedepth = 15),
          silent=F)

m1

pushover(message=as.character(Sys.time()),
         title="brms finished!",
         app="aywd2zcms7exp9b1bzuu4k2kucvgw7",
         user="us7sp3iysh1q843wedps638pcvepib")


formula_2 <- bf(W ~ fmin(Bmax + Si * (N - (alpha1*(Bmax^(-alpha2)))), Bmax),
                Bmax + Si ~ 1 + (1 | variety/date),
                alpha1 + alpha2 ~ (1 | variety),
                nl = T)

get_prior(formula_2, family = gaussian, data = d)

priors_2 <- c(set_prior("normal(5,1)", nlpar = "alpha1", lb = 0, ub = 10),
              set_prior("normal(0.5,0.1)", nlpar = "alpha2", lb = 0, ub = 1),
              set_prior("normal(10,10)", nlpar = "Bmax", lb = 1, ub = 30),
              set_prior("normal(6,1)", nlpar = "Si", lb = 0))

m2 <- brm(formula_2, data = d, family = gaussian, prior = priors_2,
          cores = 4, chains = 4, iter = 5000, warmup = 2000,
          control = list(adapt_delta = 0.99, max_treedepth = 15),
          silent=F)

m2

pushover(message=as.character(Sys.time()),
         title="brms finished!",
         app="aywd2zcms7exp9b1bzuu4k2kucvgw7",
         user="us7sp3iysh1q843wedps638pcvepib")

# plot(m1)
# summary(m1)

# m1.data <- as.data.frame(m1) %>% as_tibble()
# 
# m1.summary <- m1.data %>%
#   summarise_all(mean) %>%
#   select(1,8:57) %>%
#   pivot_longer(2:51) %>%
#   mutate(B_Max=b_Bmax_Intercept+value)
  
  
  
  
