# pre-init ----------------------------

# Due to compatibility issues with `rstan`, `renv` must be deactivated before `brms` fit
renv::hydrate(library="/Library/Frameworks/R.framework/Versions/4.1-arm64/Resources/library",update=T,sources="/Users/mini-87621949/GitHub/cndc_bayesian_eval/renv/library/R-4.1/aarch64-apple-darwin20")
renv::deactivate()

# initialization ----------------------

library(tidyverse)
library(brms)
library(tidybayes)
library(future)

library(googledrive)

# plan(multisession)

# load data ---------------------------

data_cndc <- read_csv("data/analysis/data_cndc.csv",col_types="cccccccdcdd") 
index_cndc <- read_csv("data/analysis/index_cndc.csv",col_types="ccccccc") 

# fit model ---------------------------

model.name = "model_070621"

fit_model <- function(data=data_cndc,
                      data_index=index_cndc,
                      model.name=model.name){
  
  index.list <- data_index %>%
    pull(index)
  
  d <- data %>%
    filter(index%in%index.list)
  
  formula <- bf(W ~ fmin(Bmax + Si * (N - (alpha1*(Bmax^(-alpha2)))), Bmax),
                Bmax + Si ~ 1 + (1|index),
                alpha1 + alpha2 ~ 1 + (1|location/variety),
                nl = T)
  
  priors <- c(set_prior("normal(7.0,1.0)", class = "sd", nlpar = "Bmax", group = "index"),
              set_prior("normal(1.0,0.1)", class = "sd", nlpar = "Si", group = "index"),
              set_prior("normal(0.10,0.02)", nlpar = "alpha1", class = "sd", group ="location"),
              set_prior("normal(0.05,0.02)", nlpar = "alpha2", class = "sd", group ="location"),
              set_prior("normal(0.05,0.01)", nlpar = "alpha1", class = "sd", group ="location:variety"),
              set_prior("normal(0.02,0.01)", nlpar = "alpha2", class = "sd", group ="location:variety"),
              set_prior("normal(8.0,0.1)", nlpar = "Bmax", lb = 1),
              set_prior("normal(6.0,0.1)", nlpar = "Si", lb = 0),
              set_prior("normal(5.3,0.1)", nlpar = "alpha1", lb = 0),
              set_prior("normal(0.40,0.01)", nlpar = "alpha2", lb = 0, ub = 1),
              set_prior("student_t(3,1.0,0.1)", class = "sigma"))
  
  if (class(priors)[[1]]!="brmsprior"){priors <- get_prior(formula, family = gaussian, data = d)}
  
  model <- brm(formula = formula, 
             data = d, 
             family = gaussian, 
             prior = priors,
             cores = 4, 
             chains = 4, 
             iter = 10000, 
             warmup = 3000,
             control = list(adapt_delta = 0.99, 
                            max_treedepth = 15),
             seed=52624,
             file=paste("manuscript/models/",model.name,sep=""),
             silent=F)
  
  return(model)
  
}

# model %<-% fit_model()
model <- fit_model()

# post model to Google Drive ----------

f.google.drive.model.sync <- function(model.name){
  
  model.path = paste0("manuscript/models/",model.name)
  
  drive_auth(email = "bohm0072@umn.edu")
  
  drive_put(media=paste0(model.path,".rds"),
            path=as_id("https://drive.google.com/drive/folders/1sSHW3PZW1lknSbnwIMDYFTXjcPXwiMYE"),
            # path="/My Drive/Research/Publications/Chapter 4/Version 2/GitHub/cndc_bayesian_eval/manuscript/models/",
            name=paste0(model.name,".rds"))
  
}

f.google.drive.model.sync(model.name)


# postscript --------------------------
# Due to compatibility issues with `rstan`, `renv` must be re-activated after `brms` fit
renv::activate()

# end ---------------------------------
