library(tidyverse)
library(brms)

library(future)
plan(multisession)
# options(future=TRUE)
# options(future.seed=52624)

library(pushoverr)
# set_pushover_app("aywd2zcms7exp9b1bzuu4k2kucvgw7")
# set_pushover_user("us7sp3iysh1q843wedps638pcvepib")

data <- read_csv("data.csv",col_types="cccccdcdd") 

brms_fit <- function(data,model,Owner,Location,Variety){
  
  str_replace(Variety," ","")
  str_c(Variety)
  
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
    
    formula <- bf(W ~ fmin(Bmax + Si * (N - (alpha1*(Bmax^(-alpha2)))), Bmax),
                  Bmax + Si ~ 1 + (1|date),
                  alpha1 + alpha2 ~ 1,
                  nl = T)
    
    #get_prior(formula, family = gaussian, data = d)
    priors <- c(set_prior("normal(5,1)", nlpar = "alpha1", lb = 0, ub = 10),
                set_prior("normal(0.5,0.1)", nlpar = "alpha2", lb = 0, ub = 1),
                set_prior("normal(10,10)", nlpar = "Bmax", lb = 1, ub = 30),
                set_prior("normal(6,1)", nlpar = "Si", lb = 0))
    
  } else if (model=="model-2"){
    
    formula <- bf(W ~ fmin(Bmax + Si * (N - (alpha1*(Bmax^(-alpha2)))), Bmax),
                  Bmax + Si ~ 1 + (1 | variety/date),
                  alpha1 + alpha2 ~ (1 | variety),
                  nl = T)
    
    # get_prior(formula_2, family = gaussian, data = d)
    priors <- c(set_prior("normal(5,1)", nlpar = "alpha1", lb = 0, ub = 10),
                set_prior("normal(0.5,0.1)", nlpar = "alpha2", lb = 0, ub = 1),
                set_prior("normal(10,10)", nlpar = "Bmax", lb = 1, ub = 30),
                set_prior("normal(6,1)", nlpar = "Si", lb = 0))
    
  }
  
  pushover(message=model.name,
           title=paste("brms started!",as.character(Sys.time())),
           app="aywd2zcms7exp9b1bzuu4k2kucvgw7",
           user="us7sp3iysh1q843wedps638pcvepib")
  
  fit <- brm(formula, data = d, family = gaussian, prior = priors,
           cores = 4, #future = T, 
           chains = 4, iter = 5000, warmup = 2000,
           control = list(adapt_delta = 0.99, max_treedepth = 15),
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
  
  m0 <- brms_fit(data,model="model-1",Owner=c("Bohman"),Location=c("Minnesota"),Variety=c("Russet Burbank","Clearwater","Umatilla","Dakota Russet","Easton"))
  m1 <- brms_fit(data,model="model-1",Owner=c("Bohman"),Location=c("Minnesota"),Variety=c("Russet Burbank"))
  m2 <- brms_fit(data,model="model-1",Owner=c("Bohman"),Location=c("Minnesota"),Variety=c("Clearwater"))
  m3 <- brms_fit(data,model="model-1",Owner=c("Bohman"),Location=c("Minnesota"),Variety=c("Umatilla"))
  m4 <- brms_fit(data,model="model-1",Owner=c("Bohman"),Location=c("Minnesota"),Variety=c("Dakota Russet"))
  m5 <- brms_fit(data,model="model-1",Owner=c("Bohman"),Location=c("Minnesota"),Variety=c("Easton"))
  
  out <- list(m1=m1,
              m2=m2,
              m3=m3,
              m4=m4,
              m5=m5)
  
  return(out)
  
}

fits %<-% run_fits()




