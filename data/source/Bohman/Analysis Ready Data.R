f.ard <- function(){
  
  ##### Initialization #####
  library(tidyverse)
  library(lubridate)
  library(stringr)
  library(readxl)
  
  ##### Import Analysis Ready Data #####
  
  metadata_dates <- read_csv("data/source/Bohman/metadata_dates.csv", col_types="cccDDDD")
  metadata_exp <- read_csv("data/source/Bohman/metadata_exp.csv", col_types="cccccc")
  metadata_trt <- read_csv("data/source/Bohman/metadata_trt.csv", col_types="cccccc")
  metadata_trt_n_cum <- read_csv("data/source/Bohman/metadata_trt_n_cum.csv", col_types="ccccDd")
  
  obs_tissue <- read_csv("data/source/Bohman/obs_tissue.csv", col_types="ccccDccd")
  
  ##### Join Data ####
  
  join <- obs_tissue %>%
    left_join(metadata_exp, by=c("owner","study","year","plot_id")) %>%
    left_join(metadata_trt, by=c("owner","study","year","trt_id")) %>%
    left_join(metadata_trt_n_cum, by=c("owner","study","year","trt_n","date")) %>%
    left_join(metadata_dates, by=c("owner","study","year")) 
  
  ##### Calculate Data #####
  
  calc <- join %>%
    mutate(dae=date-date_emerge)
  
  ##### Filter Data #####
  
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
  
  ##### Format Data #####
  
  format <- filter %>%
    select(owner,study,year,date,dae,plot_id,rep_id,trt_id,trt_n,cum_rate_n_kgha,trt_var,measure,tissue,value)
  
  ##### Export Data #####
  
  data <- format
  
  return(data)
  
  # write_csv(data,"Data/Analysis/data.csv")
  
  ##### END #####
  
}