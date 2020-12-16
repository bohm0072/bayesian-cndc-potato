library(tidyverse)
library(brms)
library(tidybayes)
# library(shinystan)

library(future)
plan(multisession)

library(pushoverr)

data_cndc <- read_csv("data/analysis/data_cndc.csv",col_types="cccccccdcdd") 
data_cndc_index <- read_csv("data/analysis/data_cndc_index.csv",col_types="ccccccc") 

run_fits <- function(){
  
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
               file=paste("brms/models/",model.name,sep=""))
    
    pushover(message=model.name,
             title=paste("brms finished!",as.character(Sys.time())),
             app="aywd2zcms7exp9b1bzuu4k2kucvgw7",
             user="us7sp3iysh1q843wedps638pcvepib")
    
    return(fit)
    
  }
  
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

  #get_prior(formula, family = gaussian, data = data_cndc)
  
  # m0002 <- brms_fit(data=data_cndc,
  #                   model.name="m0002_Bohman_Minnesota_All",
  #                   index=data_cndc_index %>%
  #                     filter(
  #                       !(owner=="Bohman"&location=="Minnesota"&variety%in%c("Clearwater","Dakota Russet","Umatilla","Easton"))
  #                     ) %>%
  #                     pull(index),
  #                   formula=bf(W ~ fmin(Bmax + Si * (N - (alpha1*(Bmax^(-alpha2)))), Bmax),
  #                              Bmax + Si ~ 1 + (1|index),
  #                              alpha1 + alpha2 ~ 1 + (1|group),
  #                              nl = T),
  #                   priors=c(set_prior("normal(3.0,0.1)", nlpar = "Bmax", class = "sd", group ="index"),
  #                            set_prior("normal(3.0,0.1)", nlpar = "Si", class = "sd", group ="index"),
  #                            set_prior("normal(0.20,0.01)", nlpar = "alpha1", class = "sd", group ="group"),
  #                            set_prior("normal(0.12,0.01)", nlpar = "alpha2", class = "sd", group ="group"),
  #                            set_prior("normal(12.0,0.1)", nlpar = "Bmax", lb = 1),
  #                            set_prior("normal(4.5,0.1)", nlpar = "Si", lb = 0),
  #                            set_prior("normal(5.0,0.5)", nlpar = "alpha1", lb = 0),
  #                            set_prior("normal(0.40,0.05)", nlpar = "alpha2", lb = 0, ub = 1),
  #                            set_prior("student_t(3,1.35,0.05)", class = "sigma")
  #                            ),
  #                   iter = 500,
  #                   warmup = 200,
  #                   adapt_delta = 0.99)

  # m0003 <- brms_fit(data=data_cndc,
  #                   model.name="m0003_Giletto_Canada_All",
  #                   index=data_cndc_index %>%
  #                     filter(owner=="Giletto") %>%
  #                     filter(location=="Canada") %>%
  #                     filter(variety%in%c("Russet Burbank","Shepody")) %>%
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
  
  # m0005 <- brms_fit(data=data_cndc,
  #                   model.name="m0005_All",
  #                   index=data_cndc_index %>%
  #                     pull(index),
  #                   formula=bf(W ~ fmin(Bmax + Si * (N - (alpha1*(Bmax^(-alpha2)))), Bmax),
  #                              Bmax + Si ~ 1 + (1|index),
  #                              alpha1 + alpha2 ~ 1,
  #                              nl = T),
  #                   priors=c(set_prior("normal(4.5,0.5)", nlpar = "alpha1", lb = 0),
  #                            set_prior("normal(0.40,0.05)", nlpar = "alpha2", lb = 0, ub = 1),
  #                            set_prior("normal(10.0,0.1)", nlpar = "Bmax", lb = 1),
  #                            set_prior("normal(4.5,0.1)", nlpar = "Si", lb = 0),
  #                            set_prior("normal(4.5,0.1)", class = "sd", nlpar = "Bmax"),
  #                            set_prior("normal(3.5,0.1)", class = "sd", nlpar = "Si")),
  #                   iter = 500,
  #                   warmup = 200,
  #                   adapt_delta = 0.99)
  
  m0006 <- brms_fit(data=data_cndc,
                    model.name="m0006_All",
                    index=data_cndc_index %>%
                      filter(
                        !(owner=="Bohman"&location=="Minnesota"&variety%in%c("Clearwater","Dakota Russet","Umatilla","Easton"))
                      ) %>%
                      pull(index),
                    formula=bf(W ~ fmin(Bmax + abs(Si) * (N - (alpha1*(Bmax^(-alpha2)))), Bmax),
                               Bmax + Si ~ 1 + (1|index),
                               alpha1 + alpha2 ~ 1 + (1|group),
                               nl = T),
                    priors=c(set_prior("normal(8.6,0.4)", class = "sd", nlpar = "Bmax"),
                             set_prior("normal(4.1,0.3)", class = "sd", nlpar = "Si"),
                             set_prior("normal(0.10,0.02)", nlpar = "alpha1", class = "sd", group ="group"),
                             set_prior("normal(0.12,0.02)", nlpar = "alpha2", class = "sd", group ="group"),
                             set_prior("normal(11.8,0.2)", nlpar = "Bmax", lb = 1),
                             set_prior("normal(5.0,0.05)", nlpar = "Si", lb = 0),
                             set_prior("normal(4.6,0.3)", nlpar = "alpha1", lb = 0),
                             set_prior("normal(0.30,0.03)", nlpar = "alpha2", lb = 0, ub = 1),
                             set_prior("student_t(3,1.1,0.4)", class = "sigma")),
                    iter = 500,
                    warmup = 200,
                    adapt_delta = 0.99)
  
}

fits %<-% run_fits()


