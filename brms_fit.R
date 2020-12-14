library(tidyverse)
library(brms)
library(tidybayes)
library(shinystan)

library(future)
plan(multisession)
# options(future=TRUE)
# options(future.seed=52624)

library(pushoverr)
# set_pushover_app("aywd2zcms7exp9b1bzuu4k2kucvgw7")
# set_pushover_user("us7sp3iysh1q843wedps638pcvepib")

data <- read_csv("data.csv",col_types="cccccdcdd") 

brms_fit <- function(data,model,Owner,Location,Variety){
  
  Variety.str <- Variety %>% sort() %>% str_replace(" ","") %>% str_c(collapse="-")
  
  model.name <- paste(model,Owner,Location,Variety.str,sep="_")
  
  f.data.filter <- function(data,Owner,Location,Variety,W_min,date_min,sd_min){
    
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
  
  d <- f.data.filter(data=data,
                     Owner=Owner,          # Required to select data to fit model to
                     Location=Location,    # Required to select data to fit model to
                     Variety=Variety,      # Required to select data to fit model to
                     W_min=1.0,            # Required to conform with convention that critical points are defined only for W > 1.0
                     date_min=3,           # Required to properly fit linear-plateu model
                     sd_min=1.0)           # Required to properly fit linear-plateau model
  
  if (model=="model-1"){
    
    brm_iter = 5000
    brm_warmup = 2000
    brm_adapt_delta = 0.99
    
    formula <- bf(W ~ fmin(Bmax + Si * (N - (alpha1*(Bmax^(-alpha2)))), Bmax),
                  Bmax + Si ~ 1 + (1|date),
                  alpha1 + alpha2 ~ 1,
                  nl = T)
    
    #get_prior(formula, family = gaussian, data = d)
    priors <- c(set_prior("normal(4.5,0.5)", nlpar = "alpha1", lb = 0, ub = 10),
                set_prior("normal(0.5,0.05)", nlpar = "alpha2", lb = 0, ub = 1),
                set_prior("normal(12,0.1)", nlpar = "Bmax", lb = 1), #"normal(12,0.3)" #"normal(12,0.6)" #"normal(12,4)"
                set_prior("normal(4.5,0.2)", nlpar = "Si", lb = 0),  #"normal(4.5,0.5)" #"normal(6,2)"
                set_prior("normal(3,0.1)", class = "sd", nlpar = "Bmax"), #"normal(3,0.3)"
                set_prior("normal(3,0.1)", class = "sd", nlpar = "Si") #"normal(3,0.3)"
                )
    
  } else if (model=="model-2"){
    
    brm_iter = 5000
    brm_warmup = 2000
    brm_adapt_delta = 0.99
    
    formula <- bf(W ~ fmin(Bmax + Si * (N - (alpha1*(Bmax^(-alpha2)))), Bmax),
                  Bmax + Si ~ 1 + (1 | variety/date),
                  alpha1 + alpha2 ~ (1 | variety),
                  nl = T)
    
    # get_prior(formula, family = gaussian, data = d)
    priors <- c(set_prior("normal(4.5,0.2)", nlpar = "alpha1", lb = 0, ub = 10),
                set_prior("normal(0.5,0.02)", nlpar = "alpha2", lb = 0, ub = 1),
                set_prior("normal(12,0.05)", nlpar = "Bmax", lb = 1), #"normal(12,0.3)" #"normal(12,0.6)" #"normal(12,4)"
                set_prior("normal(4.5,0.1)", nlpar = "Si", lb = 0),  #"normal(4.5,0.5)" #"normal(6,2)"
                set_prior("normal(0.02,0.01)", class = "sd", nlpar = "Bmax"),
                set_prior("normal(0.001,0.001)", class = "sd", nlpar = "Si"),
                set_prior("normal(0.2,0.1)", class = "sd", nlpar = "alpha1"),
                set_prior("normal(0.02,0.01)", class = "sd", nlpar = "alpha2"),
                set_prior("normal(3,0.1)", class = "sd", nlpar = "Si", group = "variety:date"),
                set_prior("normal(3,0.1)", class = "sd", nlpar = "Bmax", group = "variety:date")
                )
    
  } else if (model=="model-3"){
    
    brm_iter = 5000
    brm_warmup = 2000
    brm_adapt_delta = 0.99
    
    formula <- bf(W ~ fmin(Bmax + Si * (N - (alpha1*(Bmax^(-alpha2)))), Bmax),
                  Bmax + Si ~ 1 + (1 | variety:date),
                  alpha1 + alpha2 ~ (1 | variety),
                  nl = T)
    
    # get_prior(formula, family = gaussian, data = d)
    priors <- c(set_prior("student_t(3,3,0.01)", class = "sd", nlpar = "Si", group = "variety:date"),
                set_prior("student_t(3,3,0.01)", class = "sd", nlpar = "Bmax", group = "variety:date"),
                set_prior("student_t(3,0.20,0.05)", nlpar = "alpha1", class = "sd"), #"normal(0.12,0.005)"
                set_prior("student_t(3,0.05,0.02)", nlpar = "alpha2", class = "sd"), #"normal(0.04,0.005)"
                set_prior("normal(12.00,0.05)", nlpar = "Bmax", lb = 1),
                set_prior("normal(4.50,0.10)", nlpar = "Si", lb = 0),
                set_prior("normal(4.50,0.20)", nlpar = "alpha1", lb = 0), #"normal(4.50,0.01)"
                set_prior("normal(0.50,0.02)", nlpar = "alpha2", lb = 0, ub = 1), #"normal(0.50,0.001)"
                set_prior("student_t(3,1.3,0.4)", class = "sigma"))
    
  }
  
  pushover(message=model.name,
           title=paste("brms started!",as.character(Sys.time())),
           app="aywd2zcms7exp9b1bzuu4k2kucvgw7",
           user="us7sp3iysh1q843wedps638pcvepib")
  
  fit <- brm(formula, data = d, family = gaussian, prior = priors,
           cores = 4, #future = T, 
           chains = 4, iter = brm_iter, warmup = brm_warmup, #iter = 5000, warmup = 2000,
           control = list(adapt_delta = brm_adapt_delta, max_treedepth = 15),
           # sample_prior = "only",
           silent=F, 
           seed=52624,
           file=paste("Models/",model.name,sep=""))
  
  pushover(message=model.name,
           title=paste("brms finished!",as.character(Sys.time())),
           app="aywd2zcms7exp9b1bzuu4k2kucvgw7",
           user="us7sp3iysh1q843wedps638pcvepib")
  
  return(fit)
  
}

run_fits <- function(){
  
  # m1.1 <- brms_fit(data,model="model-1",Owner=c("Bohman"),Location=c("Minnesota"),Variety=c("Russet Burbank"))
  # m1.2 <- brms_fit(data,model="model-1",Owner=c("Bohman"),Location=c("Minnesota"),Variety=c("Clearwater"))
  # m1.3 <- brms_fit(data,model="model-1",Owner=c("Bohman"),Location=c("Minnesota"),Variety=c("Umatilla"))
  # m1.4 <- brms_fit(data,model="model-1",Owner=c("Bohman"),Location=c("Minnesota"),Variety=c("Dakota Russet"))
  # m1.5 <- brms_fit(data,model="model-1",Owner=c("Bohman"),Location=c("Minnesota"),Variety=c("Easton"))
  # m1.0 <- brms_fit(data,model="model-1",Owner=c("Bohman"),Location=c("Minnesota"),Variety=c("Russet Burbank","Clearwater","Umatilla","Dakota Russet","Easton"))
  # 
  # m2.0 <- brms_fit(data,model="model-2",Owner=c("Bohman"),Location=c("Minnesota"),Variety=c("Russet Burbank","Clearwater","Umatilla","Dakota Russet","Easton"))
  # 
  m3.0 <- brms_fit(data,model="model-3",Owner=c("Bohman"),Location=c("Minnesota"),Variety=c("Russet Burbank","Clearwater","Umatilla","Dakota Russet","Easton"))
  
  out <- list(# m1.0=m1.0,
              # m1.1=m1.1,
              # m1.2=m1.2,
              # m1.3=m1.3,
              # m1.4=m1.4,
              # m1.5=m1.5,
              # m2.0=m2.0,
              m3.0=m3.0)
  
  return(out)
  
}

fits %<-% run_fits()


