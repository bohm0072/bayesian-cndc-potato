# Summarize and export data for preliminary analysis of Bayesian methods for CNDC
# Using framework from Makowski et al. (2020) - doi:10.1016/j.eja.2020.126076

##### Initialization #####
library(tidyverse)
library(lubridate)
library(stringr)

##### Rosen Data #####
# Unpublished data from Rosen lab

f.bohman <- function(){
  
  source("../Analysis/Scripts/Analysis/Analysis Ready Data.R")
  ard <- f.ard()
  
  d <- ard %>%
    filter(tissue=="wholeplant") %>%
    # filter(trt_var=="russet burbank") %>%
    pivot_wider(names_from="measure",
                values_from="value") %>%
    rename(rate_n_kgha=cum_rate_n_kgha) %>%
    rename(variety=trt_var) %>%
    group_by(owner,study,year,date,variety,rate_n_kgha) %>%
    summarize_at(vars(biomdry_Mgha,n_pct),~mean(.,na.rm=T)) %>%
    ungroup() %>%
    mutate(study="Bohman") %>%
    mutate(variety=str_to_title(variety)) %>%
    # mutate(variety="Russet Burbank") %>%
    mutate(location="Minnesota") %>%
    select(study,location,variety,rate_n_kgha,Date=date,W=biomdry_Mgha,N=n_pct)
  
  return(d)
  
}

bohman <- f.bohman()

##### Giletto Data #####
# Data from Giletto et. al (2020), Appendix A - doi:10.1016/j.eja.2020.126114

f.giletto <- function(){
  
  appendix_a <- read_excel("1-s2.0-S1161030120301210-mmc1.xlsx", sheet="Appendix A", col_names=F)
  appendix_b <- read_excel("1-s2.0-S1161030120301210-mmc1.xlsx", sheet="Appendix B", col_names=F)
  
  # data = "appendix_a"
  # # #measure = "W_sh"
  # rows = c(4:8)
  # # #cols = c("D","E","F","G")
  # n_rate = c("0","80","150","250")
  # dap = c(45,60,75,90,110)
  # year = "2003-04"
  # variety = "Innovator"
  # location = "Argentina"
  # study = "Giletto"
  
  f.get <- function(data,rows,n_rate,dap,year,variety,location,study){
    
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
          mutate(Date=ymd(paste(year,"01-01",sep="-"))) %>%
          mutate(Date=Date+as.numeric(dap)) %>%
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
          mutate(Date=ymd(paste(year,"01-01",sep="-"))) %>%
          mutate(Date=Date+as.numeric(dap)) %>%
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
          mutate(Date=ymd(paste(year,"01-01",sep="-"))) %>%
          mutate(Date=Date+as.numeric(dap)) %>%
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
          mutate(Date=ymd(paste(year,"01-01",sep="-"))) %>%
          mutate(Date=Date+as.numeric(dap)) %>%
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
      mutate(W=sum(W_sh,W_t,na.rm=T)) %>%
      mutate(N=sum(W_sh*N_sh,W_t*N_t,na.rm=T)/W) %>%
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
        year = list("2004-05","2005-06"),
        variety = list("Markies Russet"),
        location = list("Argentina"),
        study = list("Giletto")
      ),
      tibble(
        data = list("appendix_a"),
        rows = list(c(4:8),c(9:12),c(13:17),c(18:22)),
        n_rate = list(c("0","80","150","250")),
        dap = list(c(45,60,75,90,110),c(55,61,83,110),c(47,62,75,91,105),c(45,58,74,95,116)),
        year = list("2003-04","2004-05","2005-06","2006-07"),
        variety = list("Innovator"),
        location = list("Argentina"),
        study = list("Giletto")
      ),
      tibble(
        data = list("appendix_b"),
        rows = list(c(4:13),c(14:23),c(24:32),c(33:42)),
        n_rate = list(c("0","50","100","250")),
        dap = list(c(44,50,57,64,74,81,87,94,107,114),c(41,49,56,59,68,73,80,97,104,110),c(40,49,54,61,69,74,81,90,98),c(32,39,46,53,63,69,76,83,88,95)),
        year = list("1900"),
        variety = list("Shepody"),
        location = list("Canada-Jacksonville","Canada-London","Canada-Hartland","Canada-Drummond"),
        study = list("Giletto")
      ),
      tibble(
        data = list("appendix_b"),
        rows = list(c(43:52),c(53:62),c(63:71),c(72:80)),
        n_rate = list(c("0","50","100","250")),
        dap = list(c(44,50,57,64,74,81,87,94,107,114),c(41,49,56,59,68,73,80,81,97,104),c(40,49,54,61,69,74,81,90,98),c(39,46,53,63,69,76,83,88,95)),
        year = list("1900"),
        variety = list("Russet Burbank"),
        location = list("Canada-Jacksonville","Canada-London","Canada-Hartland","Canada-Drummond"),
        study = list("Giletto")
      )
    )
    
  }
  
  get.list <- f.get.list()
  
  d  <- pmap(list(get.list$data,
                  get.list$rows,
                  get.list$n_rate,
                  get.list$dap,
                  get.list$year,
                  get.list$variety,
                  get.list$location,
                  get.list$study), 
             ~f.get(data = ..1,
                    rows = ..2,
                    n_rate = ..3,
                    dap = ..4,
                    year = ..5,
                    variety = ..6,
                    location = ..7,
                    study = ..8)) %>% bind_rows()
  
}

giletto <- f.giletto()

##### Combine and Format Data #####

data <- bind_rows(
  bohman,
  giletto
)

data <- data %>%
  filter(W>=1.0)

##### Export Data #####

write_csv(data,"data.csv")

##### END #####