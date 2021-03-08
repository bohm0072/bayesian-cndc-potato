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
  
  t1_ix <- data_cndc %>%
    select(owner,location,variety,index) %>%
    distinct() %>%
    group_by(owner,location,variety) %>%
    summarize_at(vars(index),~n())
  
  t1_ob <- data_cndc %>%
    distinct() %>%
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

appx_t1 <- data_cndc
write_csv(appx_t1,"manuscript/tables/appx_table1.csv")


# end ----------------
