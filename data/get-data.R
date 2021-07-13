# Summarize and export data for subsequent analysis of Bayesian methods for CNDC
# Using framework from Makowski et al. (2020) - doi:10.1016/j.eja.2020.126076
# Data from:
#   Rosen Lab at University of Minnesota - unpublished in its entirety
#   Giletto et. al (2020), Appendix A - doi:10.1016/j.eja.2020.126114
#   Ben Abdallah et. al (2016), Table 2 - doi:10.1007/s11540-016-9331-y, and Unpublished

# Initialization -------------------------------------
library(tidyverse)
library(lubridate)
library(stringr)
library(readxl)

# Bohman Data ----------------------------------------
# Previously unpublished data in its entirety from Rosen Lab at University of Minnesota

f.bohman <- function(){
  
  source("data/source/Bohman/import_rosen-lab-potato-n-trials.R")
  `rosen-lab-potato-n-trials` <- `f.import_rosen-lab-potato-n-trials`()
  
  d <- `rosen-lab-potato-n-trials` %>%
    filter(tissue=="wholeplant") %>%
    pivot_wider(names_from="measure",
                values_from="value") %>%
    rename(rate_n_kgha=sum_rate_n_kgha) %>%
    rename(variety=trt_var) %>%
    group_by(owner,study,year,date,variety,trt_n,rate_n_kgha) %>%
    summarize_at(vars(biomdry_Mgha,n_pct),~mean(.,na.rm=T)) %>%
    ungroup() %>%
    mutate(owner="Bohman") %>%
    mutate(variety=str_to_title(variety)) %>%
    mutate(location="Minnesota") %>%
    select(owner,study,year,location,variety,rate_n_kgha,date,W=biomdry_Mgha,N=n_pct)
  
  # Update `variety` for uniform naming across studies evaluated
  d <- d %>%
    mutate(variety=case_when(
      variety=="Umatilla" ~ "Umatilla Russet",
      T ~ variety
    ))
  
  return(d)
  
}

bohman <- f.bohman()

# Giletto Data ---------------------------------------
# Data from Giletto et. al (2020), Appendix A - doi:10.1016/j.eja.2020.126114

f.giletto <- function(){
  
  appendix_a <- read_excel("data/source/Giletto/Giletto-2020.xlsx", sheet="Appendix A", col_names=F)
  appendix_b <- read_excel("data/source/Giletto/Giletto-2020.xlsx", sheet="Appendix B", col_names=F)
  
  # data = "appendix_a"
  # # #measure = "W_sh"
  # rows = c(4:8)
  # # #cols = c("D","E","F","G")
  # n_rate = c("0","80","150","250")
  # dap = c(45,60,75,90,110)
  # date_plant = mdy("2003-10-24")
  # year = "2003-04"
  # variety = "Innovator"
  # location = "Argentina"
  # study = "Giletto"
  
  f.get <- function(data,rows,n_rate,dap,date_plant,year,variety,location,study){
    
    data <- get(data)
    
    ExcelLetters <- c(LETTERS[1:26],paste("A",LETTERS[1:26],sep=""))
    
    f.get.W <- function(){
      
      f.get.W_sh <- function(){
        
        cols_W_sh <- c("D","E","F","G")
        
        d_W_sh <- data %>%
          select(match(cols_W_sh, ExcelLetters)) %>%
          slice(rows)
        
        colnames(d_W_sh) <- n_rate
        
        d_W_sh <- d_W_sh %>%
          mutate(dap=dap,
                 year=year,
                 variety=variety,
                 study=study,
                 location=location) %>%
          mutate_all(as.character) %>%
          pivot_longer(names_to="n_rate",
                       values_to="W_sh", #measure
                       cols=seq(1:length(cols_W_sh))) %>% 
          mutate(year=str_sub(year,1L,4L)) %>%
          # mutate(Date=ymd(paste(year,"01-01",sep="-"))) %>%
          # mutate(Date=Date+as.numeric(dap)) %>%
          mutate(Date=date_plant+as.numeric(dap)) %>%
          select(-c(dap,year)) %>%
          select(study,location,variety,n_rate,Date,W_sh)
        
        return(d_W_sh)
        
      }
      f.get.W_t <- function(){
        
        cols_W_t <- c("S","T","U","V")
        
        d_W_t <- data %>%
          select(match(cols_W_t, ExcelLetters)) %>%
          slice(rows)
        
        colnames(d_W_t) <- n_rate
        
        d_W_t <- d_W_t %>%
          mutate(dap=dap,
                 year=year,
                 variety=variety,
                 study=study,
                 location=location) %>%
          mutate_all(as.character) %>%
          pivot_longer(names_to="n_rate",
                       values_to="W_t", #measure
                       cols=seq(1:length(cols_W_t))) %>% 
          mutate(year=str_sub(year,1L,4L)) %>%
          # mutate(Date=ymd(paste(year,"01-01",sep="-"))) %>%
          # mutate(Date=Date+as.numeric(dap)) %>%
          mutate(Date=date_plant+as.numeric(dap)) %>%
          select(-c(dap,year)) %>%
          select(study,location,variety,n_rate,Date,W_t)
        
        return(d_W_t)
        
      }
      
      W_sh <- f.get.W_sh()
      W_t <- f.get.W_t()
      
      W <- bind_cols(W_sh,
                     W_t %>% select(W_t))
      
      return(W)
      
    }
    f.get.N <- function(){
      
      f.get.N_sh <- function(){
        
        cols_N_sh <- c("I","J","K","L")
        
        d_N_sh <- data %>%
          select(match(cols_N_sh, ExcelLetters)) %>%
          slice(rows)
        
        colnames(d_N_sh) <- n_rate
        
        d_N_sh <- d_N_sh %>%
          mutate(dap=dap,
                 year=year,
                 variety=variety,
                 study=study,
                 location=location) %>%
          mutate_all(as.character) %>%
          pivot_longer(names_to="n_rate",
                       values_to="N_sh", #measure
                       cols=seq(1:length(cols_N_sh))) %>% 
          mutate(year=str_sub(year,1L,4L)) %>%
          # mutate(Date=ymd(paste(year,"01-01",sep="-"))) %>%
          # mutate(Date=Date+as.numeric(dap)) %>%
          mutate(Date=date_plant+as.numeric(dap)) %>%
          select(-c(dap,year)) %>%
          select(study,location,variety,n_rate,Date,N_sh)
        
        return(d_N_sh)
        
      }
      f.get.N_t <- function(){
        
        cols_N_t <- c("X","Y","Z","AA")
        
        d_N_t <- data %>%
          select(match(cols_N_t, ExcelLetters)) %>%
          slice(rows)
        
        colnames(d_N_t) <- n_rate
        
        d_N_t <- d_N_t %>%
          mutate(dap=dap,
                 year=year,
                 variety=variety,
                 study=study,
                 location=location) %>%
          mutate_all(as.character) %>%
          pivot_longer(names_to="n_rate",
                       values_to="N_t", #measure
                       cols=seq(1:length(cols_N_t))) %>% 
          mutate(year=str_sub(year,1L,4L)) %>%
          # mutate(Date=ymd(paste(year,"01-01",sep="-"))) %>%
          # mutate(Date=Date+as.numeric(dap)) %>%
          mutate(Date=date_plant+as.numeric(dap)) %>%
          select(-c(dap,year)) %>%
          select(study,location,variety,n_rate,Date,N_t)
        
        return(d_N_t)
        
      }
      
      N_sh <- f.get.N_sh()
      N_t <- f.get.N_t()
      
      N <- bind_cols(N_sh,
                     N_t %>% select(N_t))
      
      return(N)
      
    }
    
    d <- bind_cols(f.get.W(),
                   f.get.N() %>% select(N_sh,N_t)) %>%
      mutate_at(vars(W_sh,W_t,N_sh,N_t),as.numeric) %>%
      rowwise() %>%
      mutate(W=sum(W_sh,W_t,na.rm=F)) %>%
      mutate(N=sum(W_sh*N_sh,W_t*N_t,na.rm=F)/W) %>%
      ungroup() %>%
      mutate(N=round(N,1)) %>%
      rename(rate_n_kgha=n_rate) %>%
      mutate_at(vars(rate_n_kgha),as.numeric) %>%
      select(study,location,variety,rate_n_kgha,Date,W,N)
    
    return(d)
    
  }
  
  f.get.list <- function(){
    
    bind_rows(
      tibble(
        data = list("appendix_a"),
        rows = list(c(4:8),c(9:12),c(13:17),c(18:22)),
        n_rate = list(c("0","80","150","250")),
        dap = list(c(45,60,75,90,110),c(55,61,83,110),c(47,62,75,91,105),c(45,58,74,95,116)),
        date_plant=list(ymd("2003-10-24"),ymd("2004-10-04"),ymd("2005-09-30"),ymd("2006-10-18")),
        year = list("2003-04","2004-05","2005-06","2006-07"),
        variety = list("Innovator"),
        location = list("Argentina"),
        study = list("Giletto")
      ),
      tibble(
        data = list("appendix_a"),
        rows = list(c(23:27),c(28:31),c(32:36),c(37:41)),
        n_rate = list(c("0","80","150","250")),
        dap = list(c(45,60,75,90,110),c(55,61,83,110),c(47,62,75,91,105),c(45,58,74,95,116)),
        date_plant=list(ymd("2003-10-24"),ymd("2004-10-04"),ymd("2005-09-30"),ymd("2006-10-18")),
        year = list("2003-04","2004-05","2005-06","2006-07"),
        variety = list("Gem Russet"),
        location = list("Argentina"),
        study = list("Giletto")
      ),
      tibble(
        data = list("appendix_a"),
        rows = list(c(42:45),c(46:50),c(51:55)),
        n_rate = list(c("0","80","150","250")),
        dap = list(c(55,61,83,110),c(47,62,75,91,105),c(45,58,74,95,116)),
        date_plant=list(ymd("2004-10-04"),ymd("2005-09-30"),ymd("2006-10-18")),
        year = list("2004-05","2005-06","2006-07"),
        variety = list("Umatilla Russet"),
        location = list("Argentina"),
        study = list("Giletto")
      ),
      tibble(
        data = list("appendix_a"),
        rows = list(c(56:60),c(61:64),c(65:69)),
        n_rate = list(c("0","80","150","250")),
        dap = list(c(45,60,75,90,110),c(55,61,83,110),c(47,62,75,91,105)),
        date_plant=list(ymd("2003-10-24"),ymd("2004-10-04"),ymd("2005-09-30")),
        year = list("2003-04","2004-05","2005-06"),
        variety = list("Bannock Russet"),
        location = list("Argentina"),
        study = list("Giletto")
      ),
      tibble(
        data = list("appendix_a"),
        rows = list(c(70:73),c(74:78)),
        n_rate = list(c("0","80","150","250")),
        dap = list(c(55,61,83,110),c(47,62,75,91,105)),
        date_plant=list(ymd("2004-10-04"),ymd("2005-09-30")),
        year = list("2004-05","2005-06"),
        variety = list("Markies Russet"),
        location = list("Argentina"),
        study = list("Giletto")
      ),
      tibble(
        data = list("appendix_b"),
        rows = list(c(4:13),c(14:23),c(24:32),c(33:42)),
        n_rate = list(c("0","50","100","250")),
        dap = list(c(44,50,57,64,74,81,87,94,107,114),c(41,49,56,59,68,73,80,97,104,110),c(40,49,54,61,69,74,81,90,98),c(32,39,46,53,63,69,76,83,88,95)),
        date_plant=list(ymd("1995-06-02"),ymd("1995-05-20"),ymd("1997-06-05"),ymd("1997-05-30")),
        year = list("1995","1995","1997","1997"),#list("1900"),
        variety = list("Shepody"),
        location = list("Canada"), #list("Canada-Jacksonville","Canada-London","Canada-Hartland","Canada-Drummond"),
        study = list("Jacksonville","London","Hartland","Drummond") #list("Giletto")
      ),
      tibble(
        data = list("appendix_b"),
        rows = list(c(43:52),c(53:62),c(63:71),c(72:80)),
        n_rate = list(c("0","50","100","250")),
        dap = list(c(44,50,57,64,74,81,87,94,107,114),c(41,49,56,59,68,73,80,81,97,104),c(40,49,54,61,69,74,81,90,98),c(39,46,53,63,69,76,83,88,95)),
        date_plant=list(ymd("1995-06-02"),ymd("1995-05-20"),ymd("1997-06-05"),ymd("1997-05-30")),
        year = list("1995","1995","1997","1997"),#list("1900"),
        variety = list("Russet Burbank"),
        location = list("Canada"), #list("Canada-Jacksonville","Canada-London","Canada-Hartland","Canada-Drummond"),
        study = list("Jacksonville","London","Hartland","Drummond") #list("Giletto")
      )
    )
    
  }
  
  get.list <- f.get.list()
  
  d  <- pmap(list(get.list$data,
                  get.list$rows,
                  get.list$n_rate,
                  get.list$dap,
                  get.list$date_plant,
                  get.list$year,
                  get.list$variety,
                  get.list$location,
                  get.list$study), 
             ~f.get(data = ..1,
                    rows = ..2,
                    n_rate = ..3,
                    dap = ..4,
                    date_plant = ..5,
                    year = ..6,
                    variety = ..7,
                    location = ..8,
                    study = ..9)) %>% bind_rows()
  
  d <- d %>%
    mutate(owner="Giletto") %>%
    # mutate(study=case_when(location == "Argentina" ~ "Giletto",
    #                        str_detect(location,"Canada")==T ~ "Belanger",
    #                        T ~ NA_character_)) %>%
    mutate(year=as.character(year(Date))) %>%
    rename(date=Date) %>%
    select(owner,study,year,location,variety,rate_n_kgha,date,W,N)
  
  return(d)
  
}

giletto <- f.giletto()

# Ben Abdallah Data ----------------------------------
# Data from Ben Abdallah et. al (2016), Table 2 - doi:10.1007/s11540-016-9331-y, and Unpublished

f.ben_abdallah <- function(){
  
  d <- read_xlsx("data/source/Ben Abdallah/Ben Abdallah-2016.xlsx", sheet=2, col_names=T) %>%
    mutate(owner="Ben Abdallah",
           location="Belgium") %>%
    rename(year=Year,
           date=Date,
           variety=Cultivar,
           study=Site,
           rate_n_kgha=`Applied N (kg/ha)`,
           W=`Total Biomass        (t DM/ha)`,
           N=`N%`) %>%
    mutate_at(vars(date), as_date) %>%
    mutate_at(vars(owner,study,year,location,variety), as.character) %>%
    mutate_at(vars(rate_n_kgha,W,N), as.numeric) %>%
    select(owner,study,year,location,variety,rate_n_kgha,date,W,N)
  
  return(d)
  
}

ben_abdallah <- f.ben_abdallah()

# Combine Data ---------------------------------------

data <- bind_rows(
  bohman,
  giletto,
  ben_abdallah
)

##### Filter Data for CNDC Fit #####

f.cndc <- function(data){
  
  index <- data %>%
    select(owner,study,location,variety,date) %>%
    distinct() %>%
    arrange(owner,study,location,variety,date) %>%
    mutate(index=row_number()) %>%
    relocate(index,.before="owner")
  
  index <- left_join(
    index,
    index %>%
      select(owner,location,variety) %>%
      distinct() %>%
      mutate(group=row_number()),
    by=c("owner","location","variety")
  )  %>%
    relocate(group,.after="index")
  
  data.0 <- data %>%
    left_join(index,by=c("owner","study","location","variety","date"))  %>%
    relocate(c(index,group),.before="owner")
  
  # Summarize number of dates for each group
  data.0.sum <- data.0 %>%
    select(group,date) %>%
    distinct() %>%
    group_by(group) %>%
    count() %>%
    ungroup() %>%
    rename(count_0=n)
  
  # Screening criteria
  #   W >= 1.0 Mg/ha for >= 3 points
  
  data.1.list <- data.0 %>%
    filter(W>=1) %>%
    filter(N>=0) %>%
    drop_na() %>%
    group_by(index) %>%
    count() %>%
    ungroup() %>%
    filter(n>=3) %>%
    select(-n)
  
  data.1 <- left_join(data.1.list,
                      data.0,
                      by = c("index"))
  
  data.1.sum <- data.1 %>%
    select(group,date) %>%
    distinct() %>%
    group_by(group) %>%
    count() %>%
    ungroup() %>%
    rename(count_1=n)
  
  # Combined summary
  
  index.sum <- index %>%
    select(group,owner,location,variety) %>%
    distinct() %>%
    left_join(data.0.sum, c("group")) %>%
    left_join(data.1.sum, c("group")) %>%
    rename("count.index.data"="count_0") %>%
    rename("count.index.data_cndc"="count_1")
  
  # data.cndc <- data.2
  data.cndc <- data.1
  
  index.cndc <- data.cndc %>%
    select(index,group,owner,study,location,variety,date) %>%
    distinct() %>%
    arrange(group,index)
  
  out <- list(data=data,
              data.cndc=data.cndc,
              index=index,
              index.cndc=index.cndc,
              index.sum=index.sum)
  
  return(out)
  
}

d <- f.cndc(data)

##### Export Data #####

write_csv(d$data, "data/analysis/data.csv")
write_csv(d$data.cndc, "data/analysis/data_cndc.csv")
write_csv(d$index, "data/analysis/index.csv")
write_csv(d$index.cndc, "data/analysis/index_cndc.csv")
write_csv(d$index.sum, "data/analysis/index_sum.csv")

##### END #####