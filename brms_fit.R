library(tidyverse)
library(brms)
library(tidybayes)
library(shinystan)

library(future)
plan(multisession)

library(pushoverr)

data_cndc <- read_csv("data_cndc.csv",col_types="ccccccdcdd") 
data_cndc_index <- read_csv("data_cndc_index.csv",col_types="cccccc") 

brms_fit <- function(data,model.name,index.list,formula,priors="default",iter,warmup,adapt_delta){
  
  d <- data %>%
    filter(index%in%index.list)
  
  if (class(priors)[[1]]!="brmsprior"){priors <- get_prior(formula, family = gaussian, data = d)}
  
  pushover(message=model.name,
           title=paste("brms started!",as.character(Sys.time())),
           app="aywd2zcms7exp9b1bzuu4k2kucvgw7",
           user="us7sp3iysh1q843wedps638pcvepib")
  
  fit <- brm(formula = formula, data = d, family = gaussian, prior = priors,
             cores = 4, chains = 4, iter = iter, warmup = warmup,
             control = list(adapt_delta = adapt_delta, max_treedepth = 15),
             seed=52624,
             file=paste("Models/",model.name,sep=""))
  
  pushover(message=model.name,
           title=paste("brms finished!",as.character(Sys.time())),
           app="aywd2zcms7exp9b1bzuu4k2kucvgw7",
           user="us7sp3iysh1q843wedps638pcvepib")
  
  return(fit)
  
}

run_fits <- function(){
  
  # m0001 <- brms_fit(data=data_cndc,
  #                   model.name="m0001_Bohman_Minnesota_RussetBurbank",
  #                   index=data_cndc_index %>%
  #                     filter(owner=="Bohman") %>%
  #                     filter(location=="Minnesota") %>%
  #                     filter(variety=="Russet Burbank") %>%
  #                     pull(index),
  #                   formula=bf(W ~ fmin(Bmax + Si * (N - (alpha1*(Bmax^(-alpha2)))), Bmax),
  #                              Bmax + Si ~ 1 + (1|index),
  #                              alpha1 + alpha2 ~ 1,
  #                              nl = T),
  #                   priors=c(set_prior("normal(5.0,0.5)", nlpar = "alpha1", lb = 0),
  #                            set_prior("normal(0.40,0.05)", nlpar = "alpha2", lb = 0, ub = 1),
  #                            set_prior("normal(12.0,0.1)", nlpar = "Bmax", lb = 1),
  #                            set_prior("normal(4.5,0.1)", nlpar = "Si", lb = 0),
  #                            set_prior("normal(3.0,0.1)", class = "sd", nlpar = "Bmax"),
  #                            set_prior("normal(3.0,0.1)", class = "sd", nlpar = "Si")),
  #                   iter = 5000,
  #                   warmup = 2000,
  #                   adapt_delta = 0.99)
  # 
  # m0002 <- brms_fit(data=data_cndc,
  #                   model.name="m0002_Bohman_Minnesota_All",
  #                   index=data_cndc_index %>%
  #                     filter(owner=="Bohman") %>%
  #                     filter(location=="Minnesota") %>%
  #                     filter(variety%in%c("Clearwater","Dakota Russet","Easton","Russet Burbank","Umatilla")) %>%
  #                     pull(index),
  #                   formula=bf(W ~ fmin(Bmax + Si * (N - (alpha1*(Bmax^(-alpha2)))), Bmax),
  #                              Bmax + Si ~ 1 + (1|index),
  #                              alpha1 + alpha2 ~ 1,
  #                              nl = T),
  #                   priors=c(set_prior("normal(5.0,0.5)", nlpar = "alpha1", lb = 0),
  #                            set_prior("normal(0.40,0.05)", nlpar = "alpha2", lb = 0, ub = 1),
  #                            set_prior("normal(12.0,0.1)", nlpar = "Bmax", lb = 1),
  #                            set_prior("normal(4.5,0.1)", nlpar = "Si", lb = 0),
  #                            set_prior("normal(3.0,0.1)", class = "sd", nlpar = "Bmax"),
  #                            set_prior("normal(3.0,0.1)", class = "sd", nlpar = "Si")),
  #                   iter = 5000,
  #                   warmup = 2000,
  #                   adapt_delta = 0.99)

  m0003 <- brms_fit(data=data_cndc,
                    model.name="m0003_Giletto_Canada_All",
                    index=data_cndc_index %>%
                      filter(owner=="Giletto") %>%
                      filter(location=="Canada") %>%
                      filter(variety%in%c("Russet Burbank","Shepody")) %>%
                      pull(index),
                    formula=bf(W ~ fmin(Bmax + Si * (N - (alpha1*(Bmax^(-alpha2)))), Bmax),
                               Bmax + Si ~ 1 + (1|index),
                               alpha1 + alpha2 ~ 1,
                               nl = T),
                    priors=c(set_prior("normal(5.0,0.5)", nlpar = "alpha1", lb = 0),
                             set_prior("normal(0.40,0.05)", nlpar = "alpha2", lb = 0, ub = 1),
                             set_prior("normal(12.0,0.1)", nlpar = "Bmax", lb = 1),
                             set_prior("normal(4.5,0.1)", nlpar = "Si", lb = 0),
                             set_prior("normal(3.0,0.1)", class = "sd", nlpar = "Bmax"),
                             set_prior("normal(3.0,0.1)", class = "sd", nlpar = "Si")),
                    iter = 5000,
                    warmup = 2000,
                    adapt_delta = 0.99)
  
  # m0004 <- brms_fit(data=data_cndc,
  #                   model.name="m0004_Giletto_Argentina_All",
  #                   index=data_cndc_index %>% 
  #                     filter(owner=="Giletto") %>%
  #                     filter(location=="Argentina") %>%
  #                     filter(variety%in%c("Bannock Russet","Gem Russet","Innovator","Markies Russet","Umatilla Russet")) %>%
  #                     pull(index),
  #                   formula=bf(W ~ fmin(Bmax + Si * (N - (alpha1*(Bmax^(-alpha2)))), Bmax),
  #                              Bmax + Si ~ 1 + (1|index),
  #                              alpha1 + alpha2 ~ 1,
  #                              nl = T),
  #                   priors=c(set_prior("normal(5.0,0.5)", nlpar = "alpha1", lb = 0),
  #                            set_prior("normal(0.40,0.05)", nlpar = "alpha2", lb = 0, ub = 1),
  #                            set_prior("normal(14.0,0.1)", nlpar = "Bmax", lb = 1),
  #                            set_prior("normal(4.5,0.1)", nlpar = "Si", lb = 0),
  #                            set_prior("normal(6.0,0.1)", class = "sd", nlpar = "Bmax"),
  #                            set_prior("normal(3.0,0.1)", class = "sd", nlpar = "Si")),
  #                   iter = 5000,
  #                   warmup = 2000,
  #                   adapt_delta = 0.99)
  
}

fits %<-% run_fits()

#test <- function(){
  brm_iter = 5000
  brm_warmup = 2000
  brm_adapt_delta = 0.99
  
  if(model=="model-base"){
    
    formula <- bf(W ~ fmin(Bmax + Si * (N - (alpha1*(Bmax^(-alpha2)))), Bmax),
                  Bmax + Si ~ 1 + (1|index),
                  alpha1 + alpha2 ~ 1,
                  nl = T)
    
    #get_prior(formula, family = gaussian, data = data_cndc)
    
  } else
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
      
    } else 
      if (model=="model-2"){
        
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
        
      } else 
        if (model=="model-3"){
          
          brm_iter = 5000
          brm_warmup = 2000
          brm_adapt_delta = 0.99
          
          formula <- bf(W ~ fmin(Bmax + Si * (N - (alpha1*(Bmax^(-alpha2)))), Bmax),
                        Bmax + Si ~ 1 + (1 | variety:date),
                        alpha1 + alpha2 ~ (1 | variety),
                        nl = T)
          
          # get_prior(formula, family = gaussian, data = d)
          priors <- c(set_prior("student_t(3,3,0.1)", class = "sd", nlpar = "Si", group = "variety:date"),
                      set_prior("student_t(3,3,0.1)", class = "sd", nlpar = "Bmax", group = "variety:date"),
                      set_prior("normal(0.10,0.05)", nlpar = "alpha1", class = "sd"), #"normal(0.12,0.005)"
                      set_prior("normal(0.05,0.02)", nlpar = "alpha2", class = "sd"), #"normal(0.04,0.005)"
                      set_prior("normal(12.00,0.05)", nlpar = "Bmax", lb = 1),
                      set_prior("normal(4.50,0.10)", nlpar = "Si", lb = 0),
                      set_prior("normal(4.50,0.10)", nlpar = "alpha1", lb = 0), #"normal(4.50,0.01)"
                      set_prior("normal(0.50,0.01)", nlpar = "alpha2", lb = 0, ub = 1), #"normal(0.50,0.001)"
                      set_prior("student_t(3,0,1)", class = "sigma"))
          
        } else 
          if (model=="model-4"){
            
            brm_iter = 5000
            brm_warmup = 2000
            brm_adapt_delta = 0.99
            
            formula <- bf(W ~ fmin(Bmax + Si * (N - (alpha1*(Bmax^(-alpha2)))), Bmax),
                          Bmax + Si ~ 1 + (1 | date),
                          alpha1 + alpha2 ~ 1 + (1 | variety),
                          nl = T)
            
            # get_prior(formula, family = gaussian, data = d)
            priors <- c(set_prior("normal(3,0.1)", class = "sd", nlpar = "Si", group = "date"),
                        set_prior("normal(3,0.1)", class = "sd", nlpar = "Bmax", group = "date"),
                        set_prior("normal(0.12,0.04)", class = "sd", nlpar = "alpha1", group = "variety"),
                        set_prior("normal(0.06,0.02)", class = "sd", nlpar = "alpha2", group = "variety"),
                        set_prior("normal(12.00,0.05)", nlpar = "Bmax", lb = 1),
                        set_prior("normal(4.50,0.10)", nlpar = "Si", lb = 0),
                        set_prior("normal(4.50,0.01)", nlpar = "alpha1", lb = 0), 
                        set_prior("normal(0.50,0.001)", nlpar = "alpha2", lb = 0, ub = 1),
                        set_prior("student_t(3,0,1)", class = "sigma"))
            
          } else 
            if (model=="model-5"){
              
              brm_iter = 500
              brm_warmup = 200
              brm_adapt_delta = 0.99
              
              formula <- bf(W ~ fmin(Bmax + Si * (N - (alpha1*(Bmax^(-alpha2)))), Bmax),
                            Bmax + Si ~ 1 + (1 | date/variety),
                            alpha1 + alpha2 ~ 1 + (1 | variety),
                            nl = T)
              
              # get_prior(formula, family = gaussian, data = d)
              priors <- c(set_prior("student_t(3,0,1)", class = "sd", nlpar = "Si", group = "date"),
                          set_prior("student_t(3,0,1)", class = "sd", nlpar = "Bmax", group = "date"),
                          set_prior("normal(3,0.1)", class = "sd", nlpar = "Si", group = "date:variety"),
                          set_prior("normal(3,0.1)", class = "sd", nlpar = "Bmax", group = "date:variety"),
                          set_prior("normal(0.12,0.04)", class = "sd", nlpar = "alpha1", group = "variety"),
                          set_prior("normal(0.06,0.02)", class = "sd", nlpar = "alpha2", group = "variety"),
                          set_prior("normal(12.00,0.05)", nlpar = "Bmax", lb = 1),
                          set_prior("normal(4.50,0.10)", nlpar = "Si", lb = 0),
                          set_prior("normal(4.50,0.01)", nlpar = "alpha1", lb = 0), 
                          set_prior("normal(0.50,0.001)", nlpar = "alpha2", lb = 0, ub = 1),
                          set_prior("student_t(3,0,1)", class = "sigma"))
              
            }
}


