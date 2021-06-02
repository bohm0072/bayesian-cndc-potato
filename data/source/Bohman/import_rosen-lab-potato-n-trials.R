`f.import_rosen-lab-potato-n-trials` <- function(){
  
  # Given the `rosen-lab-potato-n-trials` GitHub repository
  # (https://github.com/bohm0072/rosen-lab-potato-n-trials), 
  # import data relevant to the analysis in this manuscript
  
  # Init ------------------------------------------------
  library(tidyverse)
  library(lubridate)
  library(stringr)
  
  source("data/source/Bohman/github.import.R")
 
  # Import ----------------------------------------------
  
  # If necessary create GitHub personal access token [PAT]; Note: this only needs to be done once.
  # source("data/source/Bohman/github.token.R")
  
  # Token, once created, can be retrieved from cache with this function
  token <- gitcreds::gitcreds_get()$password
  
  # Commit or Branch reference for GitHub repo version used in the analysis
  ref = read_json("data/source/Bohman/github.ref.json")$ref[[1]]
  
  account = "bohm0072"
  repo = "rosen-lab-potato-n-trials"
  
  dates <- f.import.github(account,repo,"Tables/dates.csv",ref,token) %>% read_csv(col_types="cccDDDD")
  experiments <- f.import.github(account,repo,"Tables/experiments.csv",ref,token) %>% read_csv(col_types="ccccc")
  trt <- f.import.github(account,repo,"Tables/trt.csv",ref,token) %>% read_csv(col_types="cccccccc")
  trt_n_sum <- f.import.github(account,repo,"Tables/trt_n_sum.csv",ref,token) %>% read_csv(col_types="ccccd")
  obs_tissue <- f.import.github(account,repo,"Tables/obs_tissue.csv",ref,token) %>% read_csv(col_types="ccccDccd")
  
  # Join Data -------------------------------------------
  
  join <- obs_tissue %>%
    left_join(experiments, by=c("owner","study","year","plot_id")) %>%
    left_join(trt, by=c("owner","study","year","trt_id")) %>%
    left_join(trt_n_sum, by=c("owner","study","year","trt_n")) %>%
    left_join(dates, by=c("owner","study","year")) 
  
  # Calculate Data --------------------------------------
  
  calc <- join %>%
    mutate(dae=date-date_emerge)
  
  # Filter Data -----------------------------------------
  
  # Filter out data that is not relevant to the critical N curve derivation
  # i.e., only studies and treatments with in-season whole plant sampling protocols
  
  filter <- calc %>%
    filter(
      (owner=="rosen lab"&
         study=="blcmr"&
         year%in%c("1991","1992")&
         trt_var%in%c("russet burbank")&
         trt_n%in%c("n.1","n.2","n.3","n.4","n.5","n.6","n.7","n.8","n.9","n.10"))|
        (owner=="rosen lab"&
           study=="vxn"&
           year%in%c("2014","2015")&
           trt_var%in%c("russet burbank","easton","dakota russet")&
           trt_n%in%c("n.1","n.2","n.3","n.4","n.5","n.x"))|
        (owner=="rosen lab"&
           study=="chloro"&
           year%in%c("2016")&
           trt_var%in%c("russet burbank")&
           trt_n%in%c("n.1","n.2","n.3","n.4"))|
        (owner=="rosen lab"&
           study=="sanjay"&
           year%in%c("2018","2019")&
           trt_var%in%c("russet burbank","clearwater","umatilla")&
           trt_n%in%c("n.a","n.b","n.c"))|
        (owner=="rosen lab"&
           study=="nni"&
           year%in%c("2019")&
           trt_var%in%c("russet burbank")&
           trt_n%in%c("n.1","n.2","n.3","n.4","n.5","n.6","n.7","n.8"))|
        (owner=="rosen lab"&
           study=="esn-nni"&
           year%in%c("2020")&
           trt_var%in%c("russet burbank")&
           trt_n%in%c("n.1","n.6","n.11","n.12","n.13","n.14","n.15","n.16"))
    )
  
  # Format Data -----------------------------------------
  
  format <- filter %>%
    select(owner,study,year,date,dae,plot_id,trt_id,trt_n,sum_rate_n_kgha,trt_var,measure,tissue,value)
  
  # Export Data -----------------------------------------
  
  data <- format
  
  return(data)
  
  # END -------------------------------------------------
   
}