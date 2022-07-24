# initialization -----------------
library(tidyverse)

# read in data -------------------
data_cndc <- read_csv("data/analysis/data_cndc.csv",col_types="cccccccdcdd")

# Table 1 ------------
# Summary of experimental data used in this study

f.t1 <- function(){
 
  t1_yr <- data_cndc %>%
    select(owner,study,location,variety,year) %>%
    distinct() %>%
    group_by(owner,study,location,variety) %>%
    summarize_at(vars(year),~n()) %>%
    ungroup() %>%
    select(-study) %>%
    group_by(owner,location,variety) %>%
    summarize_at(vars(year),~sum(.)) 
  
  # Adjust for growing season overlapping multiple calendar years for Argentina
  t1_yr <- t1_yr %>%
    mutate(year=case_when(
      location=="Argentina" ~ as.integer(year - 1),
      T ~ year
    ))
  
  t1_ix <- data_cndc %>%
    select(owner,location,variety,index) %>%
    distinct() %>%
    group_by(owner,location,variety) %>%
    summarize_at(vars(index),~n())
  
  t1_ob <- data_cndc %>%
    distinct() %>%
    drop_na() %>%
    group_by(owner,location,variety) %>%
    summarize_at(vars(W),~n())
  
  t1 <- t1_yr %>%
    left_join(t1_ix, by=c("owner","location","variety")) %>%
    left_join(t1_ob, by=c("owner","location","variety"))
  
  return(t1)
  
}
t1 <- f.t1()
write_csv(t1,"manuscript/tables/table1.csv")


# Appendix Table 1 ------------
# This is data_cndc, the set of experimental data used to populate the model
# Formatted for presentation...

f.appx_t1 <- function(){
  appx_t1 <- data_cndc %>%
    select(c(location,variety,index,date,study,year,rate_n_kgha,W,N)) %>%
    drop_na() %>%
    mutate(study=case_when(
      location=="Minnesota" & study=="blcmr" ~ "MN-1",
      location=="Minnesota" & study=="vxn" ~ "MN-2",
      location=="Minnesota" & study=="chloro" ~ "MN-3",
      location=="Minnesota" & study=="sanjay" ~ "MN-4",
      location=="Minnesota" & study=="nni" ~ "MN-5",
      location=="Minnesota" & study=="esn-nni" ~ "MN-6",
      location=="Argentina" & study=="Giletto" ~ "", #NA_character_,
      T ~ study
    )) %>%
    # mutate(rate_n_kgha=round(rate_n_kgha,0)) %>%
    mutate(across(.cols=c(rate_n_kgha),.fns=~format(round(., 0), nsmall = 0)))%>%
    # mutate(across(.cols=c(W,N),.fns=~round(.,3))) %>%
    mutate(across(.cols=c(W,N),.fns=~format(round(., 2), nsmall = 2)))%>%
    rename(
      biomass_Mgha = W,
      n_pct = N
    ) %>%
    mutate(n=row_number()) %>%
    relocate(n)
}
appx_t1 <- f.appx_t1()
write_csv(appx_t1,"manuscript/tables/tableS1.csv")

# end ----------------
