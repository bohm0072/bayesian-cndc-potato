library(tidyverse)
library(lubridate)

f.import <- function(var,loc){
  
  path = paste0("data/source/Weather Spark/Digitized/",var,"/",loc,".csv")
  
  import <- read_csv(path,col_names=F)
  
  data <- bind_rows(
    import,
    import %>% slice_head() %>% mutate(X1=X1+365),
    import %>% slice_tail() %>% mutate(X1=X1-365)
  ) %>%
    approx(xout=c(1:365)) %>%
    as_tibble() %>%
    rename("julian_day"=1,value=2) %>%
    mutate(date=as_date(julian_day-1,origin=mdy(01011900))) %>%
    mutate(month=month(date)) %>%
    mutate(location=loc) %>%
    mutate(measure=var) %>%
    select(date,month,location,measure,value)
  
  return(data)
  
}

var.list <- list("RAIN","SNOW","SRAD","LIGHT","TMAX","TMIN")

import <- bind_rows(
  map(.x=var.list,.f=~f.import(.,loc="Becker")) %>% bind_rows(),
  map(.x=var.list,.f=~f.import(.,loc="Balcarce")) %>% bind_rows(),
  map(.x=var.list,.f=~f.import(.,loc="Gembloux")) %>% bind_rows(),
  map(.x=var.list,.f=~f.import(.,loc="Saint-Leonard")) %>% bind_rows()
)

data <- import %>% 
  pivot_wider(names_from=measure,
              values_from=value) %>%
  mutate(across(.cols=c(RAIN,SNOW),.fns=~.*25.4/31)) %>% #Convert inch to mm
  mutate(PRECIP=RAIN+SNOW) %>% #Calc. total Precip as Rain + Snow
  select(-c(RAIN,SNOW)) %>%
  mutate(across(.cols=c(TMAX,TMIN),.fns=~(.-32)*5/9)) %>% #Convert deg F to deg C
  mutate(across(.cols=c(SRAD),.fns=~.*3.6)) %>% #Convert kWh/m2 to MJ/m2
  rowwise() %>%
  mutate(TMEAN=mean(c(TMAX,TMIN))) %>% #Calc. daily mean Temp
  ungroup() %>%
  rowwise() %>%
  mutate(GDD=max(min(TMEAN,30)-7,0)) %>% #Max Temp of 30 deg C; Base Temp of 10 deg C
  ungroup() %>%
  mutate(TDIFF=TMAX-TMIN) %>% #Calc. diurnal temperature difference
  select(date,month,location,TMAX,TMIN,TMEAN,TDIFF,SRAD,LIGHT,GDD,PRECIP)

# data %>%
#   ggplot() +
#   geom_line(aes(date,TDIFF,color=location))

sum <- data %>%
  filter((location=="Becker"&date%in%c(ymd("1900-05-01"):ymd("1900-09-15")))|
         (location=="Gembloux"&date%in%c(ymd("1900-04-20"):ymd("1900-09-20")))|
         (location=="Saint-Leonard"&date%in%c(ymd("1900-06-01"):ymd("1900-10-10")))|
         (location=="Balcarce"&date%in%c(ymd("1900-10-10"):ymd("1900-12-31"),ymd("1900-01-01"):ymd("1900-03-10")))) %>%
  group_by(location) %>%
  summarize(DAYS=n(),
            GDD=sum(GDD),
            TMEAN=mean(TMEAN),
            PRECIP=sum(PRECIP),
            SRAD=mean(SRAD),
            LIGHT=mean(LIGHT),
            TDIFF=mean(TDIFF),
            .groups = "drop")

write_csv(sum,"data/source/Weather Spark/Summarized/sum.csv")

sum.early <- data %>%
  filter((location=="Becker"&date%in%c((ymd("1900-05-01")+40):(ymd("1900-05-01")+70)))|
         (location=="Gembloux"&date%in%c((ymd("1900-04-20")+40):(ymd("1900-04-20")+70)))|
         (location=="Saint-Leonard"&date%in%c((ymd("1900-06-01")+40):(ymd("1900-06-01")+70)))|
         (location=="Balcarce"&date%in%c((ymd("1900-10-10")+40):(ymd("1900-10-10")+70)))) %>%
  group_by(location) %>%
  summarize(DAYS=n(),
            GDD=sum(GDD),
            TMAX=mean(TMAX),
            TMIN=mean(TMIN),
            TMEAN=mean(TMEAN),
            PRECIP=sum(PRECIP),
            SRAD=mean(SRAD),
            LIGHT=mean(LIGHT),
            TDIFF=mean(TDIFF),
            .groups = "drop")
