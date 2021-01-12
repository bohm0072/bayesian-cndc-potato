# initialization -----------------

library(tidyverse)
library(brms)
library(tidybayes)
library(gtable)
library(grid)
library(gridExtra)
library(egg)

# read in data -------------------
data <- read_csv("data/analysis/data_cndc.csv",col_types="cccccccdcdd"); data_cndc = data
# data_cndc_index <- read_csv("data/analysis/data_cndc_index.csv",col_types="ccccccc"); #data_index = data_cndc_index

# read in model fit results ------------------

model <- readRDS("brms/models/m0011.rds"); model

# the fmin() function used in Stan isn't defined in R, so we need to create it so that when we try to use brms to make predictions, it knows what to do with the fmin()
fmin <- function(x,y){
  pmin(x,y)
}

# figure 1 - plot of alpha parameter values ---------------

# parm = "alpha1"

f.fig1 <- function(model,parm){
  
  var1 <- paste("b_",parm,"_Intercept",sep="")
  var2 <- paste("r_location__",parm,sep="")
  var3 <- paste("r_location:variety__",parm,sep="")
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
  
  p1 <- d %>%
    ggplot(aes(x = !!sym(var5), y = reorder(variety, desc(variety)))) +
    stat_halfeye() +
    facet_grid(location~., scales = "free_y", space = "free") +
    coord_cartesian(xlim=limits.x) +
    # theme_classic() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank()) 
    
  g1 <- ggplotGrob(p1)
  fg1 <- gtable_frame(g1, height = unit(3, "null")) #debug = TRUE,
  
  p2 <- d %>%
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

fig1_a <- f.fig1(model,"alpha1")
ggsave(filename="manuscript/images/figure1_a.pdf",plot=fig1_a,height=4,width=3,units="in",scale=1.5)
fig1_b <- f.fig1(model,"alpha2")
ggsave(filename="manuscript/images/figure1_b.pdf",plot=fig1_b,height=4,width=3,units="in",scale=1.5)

# END --------------------