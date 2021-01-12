# initialization -----------------

library(tidyverse)
library(brms)
library(tidybayes)
library(gtable)
library(grid)
library(gridExtra)
library(egg)

# read in data -------------------
data <- read_csv("data/analysis/data_cndc.csv",col_types="cccccccdcdd")#; data_cndc = data
# data_cndc_index <- read_csv("data/analysis/data_cndc_index.csv",col_types="ccccccc"); #data_index = data_cndc_index

# read in model fit results ------------------

# model_old <- readRDS("brms/models/model_old.rds"); model_old
model <- readRDS("brms/models/model.rds"); model

# the fmin() function used in Stan isn't defined in R, so we need to create it so that when we try to use brms to make predictions, it knows what to do with the fmin()
fmin <- function(x,y){
  pmin(x,y)
}


# format data for results and figures ---------------

f.cndc.fit <- function(model){
  
  f.parm.fit <- function(model,parm){
    
    var1 <- paste("b_",parm,"_Intercept",sep="")
    var2 <- paste("r_location__",parm,sep="")
    var3 <- paste("r_location:variety__",parm,sep="")
    var4 <- paste("location_",parm,sep="")
    var5 <- paste("location:variety_",parm,sep="")
    
    d <- model %>%
      spread_draws(
        !!sym(var1),
        (!!sym(var2))[!!sym("location"),],
        (!!sym(var3))[!!sym("location:variety"),]
      ) %>%
      rowwise() %>%
      filter(is.na(str_match(`location:variety`,location))==F) %>%
      ungroup() %>%
      mutate(!!sym(var4) := !!sym(var1) + !!sym(var2)) %>%
      mutate(!!sym(var5) := !!sym(var1) + !!sym(var2) + !!sym(var3)) %>%
      mutate_at(vars(location,`location:variety`),as.character) 
    
    d <- d %>%
      mutate(variety=str_split(d$`location:variety`,"_",simplify=T)[,2]) %>% 
      rowwise() %>%
      mutate_at(vars(variety),~str_replace(.,"[.]"," ")) %>%
      ungroup()
    
    d <- d %>%
      arrange(location,variety) %>%
      mutate_at(vars(variety,location), ~as_factor(.)) %>%
      mutate_at(vars(variety,location), ~fct_inorder(.))
    
    return(d)
    
  }
  
  alpha1.fit <- f.parm.fit(model,"alpha1")
  alpha2.fit <- f.parm.fit(model,"alpha2")
  
  cndc.fit <- left_join(
    alpha1.fit,
    alpha2.fit,
    by = c(".chain", ".iteration", ".draw", "location", "variety", "location:variety")
  ) %>%
    relocate(all_of(c("location","variety","location:variety")),.after=.draw)
  
  
  
  return(cndc.fit)
  
}
cndc.fit <- f.cndc.fit(model)

f.cndc.fit.sum <- function(cndc.fit){
  
  cndc.fit.sum <- cndc.fit %>%
    left_join(
      expand.grid(
        .draw = cndc.fit %>% pull(.draw) %>% unique(),
        W = seq(1,40,0.1),stringsAsFactors=F
      ) %>%
        as_tibble() %>%
        arrange(as.numeric(.draw)),
      by=".draw"
    ) %>%
    mutate(`location:variety_N`=`location:variety_alpha1`*(W^(-`location:variety_alpha2`))) %>%
    mutate(`location_N`=`location_alpha1`*(W^(-`location_alpha2`))) %>%
    group_by(location,variety,location:variety,W) %>%
    summarize(qs = quantile(`location:variety_N`,c(0.05,0.50,0.95)), prob = c(0.05,0.50,0.95), .groups="drop") %>%
    pivot_wider(names_from=prob,
                names_prefix="location:variety_N_",
                values_from=qs)
  
  return(cndc.fit.sum)
  
}
cndc.fit.sum <- f.cndc.fit.sum(cndc.fit)

f.plateau.fit <- function(model,data){
  
  d0 <- data %>%
    select(index,location,variety) %>%
    left_join(
      left_join(
        model %>%
          spread_draws(b_Bmax_Intercept, r_index__Bmax[index,]) %>%
          mutate(index_Bmax = b_Bmax_Intercept + r_index__Bmax),
        model %>%
          spread_draws(b_Si_Intercept, r_index__Si[index,]) %>%
          mutate(index_Si = b_Si_Intercept + r_index__Si),
        by = c(".chain", ".iteration", ".draw", "index")) %>%
        select(.chain,.iteration,.draw,index,index_Bmax,index_Si) %>%
        rename(Bmax=index_Bmax,Si=index_Si) %>%
        mutate_at(vars(index),as.character),
      by=c("index"))
  
  d0$variety.name=str_replace(d0$variety," ",".")
  d0$`location:variety`=paste(d0$location,"_",d0$variety.name,sep="")
  d0$variety.name <- NULL
  d0 <- d0 %>% relocate(`location:variety`,.after=variety)
  
  d1 <- d0 %>% 
    left_join(
      left_join(
        model %>%
          spread_draws(b_alpha1_Intercept, `r_location__alpha1`[`location`,], `r_location:variety__alpha1`[`location:variety`,]) %>%
          rowwise() %>%
          filter(is.na(str_match(`location:variety`,location))==F) %>%
          ungroup() %>%
          mutate(`location:variety_alpha1` = b_alpha1_Intercept + r_location__alpha1 + `r_location:variety__alpha1`),
        model %>%
          spread_draws(b_alpha2_Intercept, `r_location__alpha2`[`location`,], `r_location:variety__alpha2`[`location:variety`,]) %>%
          rowwise() %>%
          filter(is.na(str_match(`location:variety`,location))==F) %>%
          ungroup() %>%
          mutate(`location:variety_alpha2` = b_alpha2_Intercept + r_location__alpha2 + `r_location:variety__alpha2`),
        by = c(".chain", ".iteration", ".draw", "location", "location:variety"="location:variety")) %>%
        select(.chain, .iteration, .draw, location, `location:variety`, `location:variety_alpha1`, `location:variety_alpha2`) %>%
        # rename(alpha1=`location:variety_alpha1`,alpha2=`location:variety_alpha2`) %>%
        mutate_at(vars(location,`location:variety`),as.character),
      by=c(".chain", ".iteration", ".draw", "location", "location:variety"="location:variety")) %>%
    mutate(Nc=`location:variety_alpha1`*(Bmax^(-`location:variety_alpha2`)))
  
  d2 <- d1 %>%
    relocate(.chain,.iteration,.draw,.before=index) %>%
    relocate(`location:variety_alpha1`,`location:variety_alpha2`,.before=Bmax) %>%
    arrange(location,variety) %>%
    mutate_at(vars(variety,location), ~as_factor(.)) %>%
    mutate_at(vars(variety,location), ~fct_inorder(.))
  
  d <- d2
  
  return(d)
  
}
plateau.fit <- f.plateau.fit(model,data)

f.Bmax.sum <- function(plateau.fit){
  
  Bmax.sum <- plateau.fit %>%
    group_by(location,variety,location:variety,index) %>%
    summarize(qs = quantile(`Bmax`,c(0.05,0.50,0.95)), prob = c(0.05,0.50,0.95), .groups="drop") %>%
    pivot_wider(names_from=prob,
                names_prefix="Bmax_",
                values_from=qs) %>%
    group_by(location,variety,location:variety) %>%
    summarize_at(vars(Bmax_0.05,Bmax_0.5,Bmax_0.95), max)
  
  return(Bmax.sum)
  
}
Bmax.sum <- f.Bmax.sum(plateau.fit)

f.plateau.fit.sum <- function(plateau.fit){
  
  plateau.fit.sum <- bind_rows(
    plateau.fit %>%
      mutate(W=Bmax,
             N=Nc) %>%
      group_by(location,variety,location:variety,index) %>%
      summarize_at(vars(W,N),~quantile(.,c(0.50)), prob = c(0.50), .groups="drop"),
    plateau.fit %>%
      mutate(W=Bmax,
             N=7.0) %>%
      group_by(location,variety,location:variety,index) %>%
      summarize_at(vars(W,N),~quantile(.,c(0.50)), prob = c(0.50), .groups="drop"),
    plateau.fit %>%
      mutate(W=Bmax + Si * (0-Nc),
             N=0) %>%
      group_by(location,variety,location:variety,index) %>%
      summarize_at(vars(W,N),~quantile(.,c(0.50)), prob = c(0.50), .groups="drop")
  )
  
  return(plateau.fit.sum)
  
}
plateau.fit.sum <- f.plateau.fit.sum(plateau.fit)

f.plot.format <- function(cndc.fit.sum,Bmax.sum,data){
  
  c <- cndc.fit.sum %>%
    left_join(Bmax.sum %>%
                select(location,variety,location:variety,Bmax_0.95),
              by=c("location","variety")) %>%
    filter(W <= Bmax_0.95) %>%
    select(-Bmax_0.95) 
  
  c <- c %>%
    filter(location%in%.location) %>%
    filter(variety%in%.variety) %>%
    mutate_at(vars(location,variety,`location:variety`), ~fct_drop(.))
  
  d <- data
  
  # d$variety.name=str_replace(d$variety," ",".")
  # d$`location:variety`=paste(d$location,":",d$variety.name,sep="")
  # d$variety.name <- NULL
  
  d$`location:variety`=paste(d$location,":",d$variety,sep="")
  
  d <- d %>% 
    relocate(`location:variety`,.after=variety) %>%
    arrange(location,variety) %>%
    mutate_at(vars(variety,location,`location:variety`), ~as_factor(.)) %>%
    mutate_at(vars(variety,location,`location:variety`), ~fct_inorder(.))
  
  d <- d %>%
    filter(location%in%.location) %>%
    filter(variety%in%.variety) %>%
    mutate_at(vars(location,variety,`location:variety`), ~fct_drop(.))
  
  p <- plateau.fit.sum %>%
    filter(location%in%.location) %>%
    filter(variety%in%.variety) %>%
    mutate_at(vars(location,variety,`location:variety`), ~fct_drop(.))
  
  # p$`location:variety`=str_replace(p$`location:variety`," ",".")
  
}

# figure 1 - plot of alpha parameter values ---------------

# parm = "alpha1"

f.fig1 <- function(cndc.fit,parm){
  
  # var1 <- paste("b_",parm,"_Intercept",sep="")
  # var2 <- paste("r_location__",parm,sep="")
  # var3 <- paste("r_location:variety__",parm,sep="")
  var4 <- paste("location_",parm,sep="")
  var5 <- paste("location:variety_",parm,sep="")
  
  limits.x <- case_when(
    parm == "alpha1" ~ c(4.0,5.5),
    parm == "alpha2" ~ c(0.1,0.7)
  )
  
  label.x <- case_when(
    parm == "alpha1" ~ expression(paste("parameter ", italic("a"), sep=" ")),
    parm == "alpha2" ~ expression(paste("parameter ", italic("b"), sep=" "))
  )
  
  # d <- model %>%
  #   spread_draws(
  #     !!sym(var1),
  #     (!!sym(var2))[!!sym("location"),],
  #     (!!sym(var3))[!!sym("location:variety"),]
  #     ) %>%
  #   rowwise() %>%
  #   filter(is.na(str_match(`location:variety`,location))==F) %>%
  #   ungroup() %>%
  #   mutate(!!sym(var4) := !!sym(var1) + !!sym(var2)) %>%
  #   mutate(!!sym(var5) := !!sym(var1) + !!sym(var2) + !!sym(var3)) %>%
  #   mutate_at(vars(location,`location:variety`),as.character) 
  # 
  # d <- d %>%
  #   mutate(variety=str_split(d$`location:variety`,"_",simplify=T)[,2]) %>% 
  #   rowwise() %>%
  #   mutate_at(vars(variety),~str_replace(.,"[.]"," ")) %>%
  #   ungroup()
  # 
  # d <- d %>%
  #   arrange(location,variety) %>%
  #   mutate_at(vars(variety,location), ~as_factor(.)) %>%
  #   mutate_at(vars(variety,location), ~fct_inorder(.))
  
  p1 <- cndc.fit %>%
    ggplot(aes(x = !!sym(var5), y = reorder(variety, desc(variety)))) +
    stat_halfeye() +
    facet_grid(location~., scales = "free_y", space = "free") +
    coord_cartesian(xlim=limits.x) +
    # theme_classic() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank()) 
    
  g1 <- ggplotGrob(p1)
  fg1 <- gtable_frame(g1, height = unit(3, "null")) #debug = TRUE,
  
  p2 <- cndc.fit %>%
    ggplot(aes(x = !!sym(var4), y = reorder(location, desc(location)))) +
    stat_halfeye() +
    coord_cartesian(xlim=limits.x) +
    labs(x=label.x) +
    # theme_classic() +
    theme(axis.title.y = element_blank())
  
  g2 <- ggplotGrob(p2)
  fg2 <- gtable_frame(g2, height = unit(1, "null"))
  
  # b <- rectGrob(gp = gpar(fill = "white"))
  
  fg <- rbind(fg1,fg2,size = "first")
  fg$widths <- unit.pmax(fg1$widths, fg2$widths)
  
  # grid.arrange(fg)
  # grid.newpage()
  # grid.rect(gp=gpar(fill="white"))
  # grid.draw(fg)
  
  return(fg)
  
}

fig1_a <- f.fig1(cndc.fit,"alpha1")
ggsave(filename="manuscript/images/figure1_a.pdf",plot=fig1_a,height=4,width=3,units="in",scale=1.5)
fig1_b <- f.fig1(cndc.fit,"alpha2")
ggsave(filename="manuscript/images/figure1_b.pdf",plot=fig1_b,height=4,width=3,units="in",scale=1.5)

# figure 2 - plot of model fit with data shown ------------------

# .location = "Minnesota"
# .variety = "Russet Burbank"

f.fig2 <- function(cndc.fit.sum,Bmax.sum,data,.location,.variety){
  
  c <- cndc.fit.sum %>%
    left_join(Bmax.sum %>%
                select(location,variety,location:variety,Bmax_0.95),
              by=c("location","variety")) %>%
    filter(W <= Bmax_0.95) %>%
    select(-Bmax_0.95) 
  
  c <- c %>%
    filter(location%in%.location) %>%
    filter(variety%in%.variety) %>%
    mutate_at(vars(location,variety,`location:variety`), ~fct_drop(.))
  
  d <- data
  
  # d$variety.name=str_replace(d$variety," ",".")
  # d$`location:variety`=paste(d$location,":",d$variety.name,sep="")
  # d$variety.name <- NULL
  
  d$`location:variety`=paste(d$location,":",d$variety,sep="")
  
  d <- d %>% 
    relocate(`location:variety`,.after=variety) %>%
    arrange(location,variety) %>%
    mutate_at(vars(variety,location,`location:variety`), ~as_factor(.)) %>%
    mutate_at(vars(variety,location,`location:variety`), ~fct_inorder(.))
  
  d <- d %>%
    filter(location%in%.location) %>%
    filter(variety%in%.variety) %>%
    mutate_at(vars(location,variety,`location:variety`), ~fct_drop(.))
  
  p <- plateau.fit.sum %>%
    filter(location%in%.location) %>%
    filter(variety%in%.variety) %>%
    mutate_at(vars(location,variety,`location:variety`), ~fct_drop(.))
  
  # p$`location:variety`=str_replace(p$`location:variety`," ",".")
  
  ggplot() +
    geom_line(data=c,aes(x=W,y=`location:variety_N_0.5`),linetype=1,alpha=1.0) +
    # geom_line(data=c,aes(x=W,y=`location:variety_N_0.05`),linetype=2,alpha=0.5) +
    # geom_line(data=c,aes(x=W,y=`location:variety_N_0.95`),linetype=2,alpha=0.5) +
    geom_line(data=p,aes(x=W,y=N),linetype=1,alpha=0.5) +
    geom_point(data=d,aes(x=W,y=N),alpha=0.33) +
    # facet_wrap(vars(`location:variety`,index)) +
    facet_wrap(vars(as.numeric(index))) +
    # facet_grid(`location:variety`~index) +
    labs(x = "W",
         y = "%N",
         title = paste(.location,.variety,sep=" - ")) + 
    coord_cartesian(xlim=c(0,NA),ylim=c(0,6.0)) +
    theme_classic()
  
}

fig2_a <- f.fig2(cndc.fit.sum,Bmax.sum,data,.location=c("Minnesota"),.variety=c("Russet Burbank"))
fig2_b <- f.fig2(cndc.fit.sum,Bmax.sum,data,.location=c("Minnesota"),.variety=c("Umatilla"))

# END --------------------