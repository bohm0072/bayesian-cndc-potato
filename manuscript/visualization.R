# initialization -----------------

library(tidyverse)
library(brms)
library(tidybayes)
library(gtable)
library(grid)
library(gridExtra)
library(egg)
library(ggExtra)
library(ggpubr)


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
    group_by(location,variety,`location:variety`,W) %>%
    summarize(qs = quantile(`location:variety_N`,c(0.05,0.50,0.95)), prob = c(0.05,0.50,0.95), .groups="drop") %>%
    pivot_wider(names_from=prob,
                names_prefix="N_",#"N_location:variety_",
                values_from=qs) %>%
    mutate_at(vars(N_0.05,N_0.5,N_0.95),as.numeric) %>%
    # mutate_at(vars(`N_location:variety_0.05`,`N_location:variety_0.5`,`N_location:variety_0.95`),as.numeric) %>%
    arrange(location,variety) %>%
    mutate_at(vars(variety,location,`location:variety`), ~as_factor(.)) %>%
    mutate_at(vars(variety,location,`location:variety`), ~fct_inorder(.))
  
  return(cndc.fit.sum)
  
}
cndc.fit.sum <- f.cndc.fit.sum(cndc.fit)

f.plateau.fit.sum <- function(plateau.fit){
  
  f.plateau.fit.quantile <- function(q){
    
    var1 <- paste("W",q,sep="_")
    var2 <- paste("N",q,sep="_")
    
    bind_rows(
      plateau.fit %>%
        mutate(W=Bmax,
               N=Nc) %>%
        group_by(location,variety,`location:variety`,index) %>%
        summarize_at(vars(W,N),~quantile(.,c(q)), prob = c(q), .groups="drop"),
        # summarize_at(vars(W,N),~quantile(.,c(0.50)), prob = c(0.50), .groups="drop"),
      plateau.fit %>%
        mutate(W=Bmax,
               N=7.0) %>%
        group_by(location,variety,`location:variety`,index) %>%
        summarize_at(vars(W,N),~quantile(.,c(q)), prob = c(q), .groups="drop"),
      plateau.fit %>%
        mutate(W=Bmax + Si * (0-Nc),
               N=0) %>%
        group_by(location,variety,`location:variety`,index) %>%
        summarize_at(vars(W,N),~quantile(.,c(q)), prob = c(q), .groups="drop")
    ) %>%
      ungroup() %>%
      mutate_at(vars(W,N),as.numeric) %>%
      rename(!!sym(var1):=W,!!sym(var2):=N) %>%
      # rename(W_0.5=W,N_0.5=N) %>%
      arrange(location,variety) %>%
      mutate_at(vars(variety,location,`location:variety`), ~as_factor(.)) %>%
      mutate_at(vars(variety,location,`location:variety`), ~fct_inorder(.))
    
  }
  
  plateau.fit.sum <- f.plateau.fit.quantile(0.5) %>%
    left_join(f.plateau.fit.quantile(0.05), c("location", "variety", "location:variety", "index")) %>%
    left_join(f.plateau.fit.quantile(0.95), c("location", "variety", "location:variety", "index"))
  
  return(plateau.fit.sum)
  
  # plateau.fit.sum <- bind_rows(
  #   plateau.fit %>%
  #     mutate(W=Bmax,
  #            N=Nc) %>%
  #     group_by(location,variety,`location:variety`,index) %>%
  #     summarize_at(vars(W,N),~quantile(.,c(0.50)), prob = c(0.50), .groups="drop"),
  #   plateau.fit %>%
  #     mutate(W=Bmax,
  #            N=7.0) %>%
  #     group_by(location,variety,`location:variety`,index) %>%
  #     summarize_at(vars(W,N),~quantile(.,c(0.50)), prob = c(0.50), .groups="drop"),
  #   plateau.fit %>%
  #     mutate(W=Bmax + Si * (0-Nc),
  #            N=0) %>%
  #     group_by(location,variety,`location:variety`,index) %>%
  #     summarize_at(vars(W,N),~quantile(.,c(0.50)), prob = c(0.50), .groups="drop")
  # ) %>%
  #   ungroup() %>%
  #   mutate_at(vars(W,N),as.numeric) %>%
  #   rename(W_0.5=W,N_0.5=N) %>%
  #   arrange(location,variety) %>%
  #   mutate_at(vars(variety,location,`location:variety`), ~as_factor(.)) %>%
  #   mutate_at(vars(variety,location,`location:variety`), ~fct_inorder(.))
  
}
plateau.fit.sum <- f.plateau.fit.sum(plateau.fit)

f.Bmax.sum <- function(plateau.fit){
  
  Bmax.sum <- plateau.fit %>%
    group_by(location,variety,`location:variety`,index) %>%
    summarize(qs = quantile(`Bmax`,c(0.05,0.50,0.95)), prob = c(0.05,0.50,0.95), .groups="drop") %>%
    pivot_wider(names_from=prob,
                names_prefix="Bmax_",
                values_from=qs) %>%
    group_by(location,variety,`location:variety`) %>%
    summarize_at(vars(Bmax_0.05,Bmax_0.5,Bmax_0.95), max) %>%
    ungroup() %>%
    arrange(location,variety) %>%
    mutate_at(vars(variety,location,`location:variety`), ~as_factor(.)) %>%
    mutate_at(vars(variety,location,`location:variety`), ~fct_inorder(.))
  
  return(Bmax.sum)
  
}
Bmax.sum <- f.Bmax.sum(plateau.fit)

f.plot.data <- function(data,cndc.fit.sum,plateau.fit.sum,Bmax.sum){
  
  c <- cndc.fit.sum %>%
    left_join(Bmax.sum %>%
                select(location,variety,`location:variety`,Bmax_0.95),
              by=c("location","variety","location:variety")) %>%
    filter(W <= Bmax_0.95) %>%
    select(-Bmax_0.95)
  
  d <- data %>%
    rowwise() %>%
    mutate(variety=str_replace(variety," ",".")) %>%
    mutate(`location:variety`=paste(location,"_",variety,sep="")) %>%
    mutate(variety=str_replace(variety,"[.]"," ")) %>%
    ungroup() %>%
    relocate(`location:variety`,.after=variety) %>%
    arrange(location,variety) %>%
    mutate_at(vars(variety,location,`location:variety`), ~as_factor(.)) %>%
    mutate_at(vars(variety,location,`location:variety`), ~fct_inorder(.))
  
  p <- plateau.fit.sum 
  
  out <- list(c=c,
              d=d,
              p=p)
  
  return(out)
  
}
plot.data <- f.plot.data(data,cndc.fit.sum,plateau.fit.sum,Bmax.sum)

# figure 1 - distribution of alpha parameter values for each parameter independently ------------------

# parm = "alpha1"

f.fig1 <- function(cndc.fit,parm){
  
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

# figure 2 - curve fits for each variety x location  ------------------

f.fig2 <- function(plot.data){
  
  c <- plot.data$c
  
  d <- plot.data$d
  
  p <- plot.data$p
  
  ggplot() +
    geom_ribbon(data=c,aes(ymin=N_0.05,ymax=N_0.95,x=W,group=`location:variety`),fill="#737373",alpha=0.66) + #color=NA,
    geom_point(data=d,aes(x=W,y=N),alpha=0.2,color="#FB9A99")+ #"#A6CEE3") +
    geom_line(data=c,aes(x=W,y=N_0.5,group=`location:variety`),linetype=1,alpha=1.0) +
    # geom_line(data=c,aes(x=W,y=N_0.05),linetype=1,color="#333333",size=0.01) +
    # geom_line(data=c,aes(x=W,y=N_0.95),linetype=1,color="#333333",size=0.01) +
    # geom_line(data=p,aes(x=W_0.5,y=N_0.5,group=index),linetype=1,alpha=0.05) +
    facet_wrap(vars(`location:variety`)) + #,scale="free_x") +
    # facet_wrap(vars(`location:variety`,index)) +
    # facet_wrap(vars(as.numeric(index))) +
    # facet_grid(`location:variety`~index) +
    labs(x = "W",
         y = "%N",
         title = paste(.location,.variety,sep=" - ")) + 
    coord_cartesian(xlim=c(0,NA),ylim=c(0,6.0)) +
    theme_classic() #+
    # scale_color_manual(values=c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#000000","#000000","#000000","#000000"))
  
}
fig2 <- f.fig2(plot.data)

# figure 3 - distribution of alpha parameters for each parameters simultaneously ----------------

# .location = "Minnesota"
# .variety = "Russet Burbank"
# .color = "#e41a1c"

f.fig3 <- function(cndc.fit,.location,.variety,.color){
  
  var1 <- expression(paste("parameter ", italic("a"), sep=" "))
  var2 <- expression(paste("parameter ", italic("b"), sep=" "))
  var3 <- paste(.location,.variety,sep=" - ")
  
  d <- cndc.fit %>%
    # select(-`location:variety`)
    filter(location %in% .location) %>%
    filter(variety %in% .variety)
  
  p1 <- ggplot(data = d, aes(x=`location:variety_alpha1`, y=`location:variety_alpha2`, color=`location:variety`)) +
    geom_point(alpha=0.01) +
    # geom_smooth(method="lm",formula=y~x) +
    geom_smooth(method="lm",formula=y~x,color="black",size=0.5) + #,color="black",linetype=2
    stat_regline_equation(color="black",size=2) +
    theme_classic() +
    scale_color_manual(values = .color) +
    # scale_color_brewer(palette = "Set1") +
    scale_x_continuous(limits=c(4.0,5.5)) + 
    scale_y_continuous(limits=c(0.01,0.7)) +
    guides(color="none") +
    labs(x=var1,
         y=var2,
         title=var3)
  
    # geom_point(data = d, aes(x=`location:variety_alpha1`, y=`location:variety_alpha2`), alpha=0.01, color="grey") +
    # facet_wrap(vars(`location:variety`)) +
    # scale_color_manual(values=c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999","#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666")) # +
    # coord_cartesian(xlim=c(4.0,5.5),ylim=c(0.1,0.7))
  
  p2 <- ggMarginal(p1, type="density", data = cndc.fit, groupColour = T)
  
  return(p2)
  
}
fig3_a <- f.fig3(cndc.fit,.location=c("Argentina"),.variety=c("Bannock Russet"),.color="#e41a1c")
ggsave(filename="manuscript/images/figure3_a.pdf",plot=fig3_a,height=1.5,width=1.5,units="in",scale=2.5)

fig3_b <- f.fig3(cndc.fit,.location=c("Argentina"),.variety=c("Gem Russet"),.color="#377eb8")
ggsave(filename="manuscript/images/figure3_b.pdf",plot=fig3_b,height=1.5,width=1.5,units="in",scale=2.5)

appx1_b <- f.appx1(plot.data,.location=c("Argentina"),.variety=c("Gem Russet"))
ggsave(filename="manuscript/images/appendix1_b.pdf",plot=appx1_b,height=5*(3.4/7),width=6,units="in",scale=1.5)

appx1_c <- f.appx1(plot.data,.location=c("Argentina"),.variety=c("Innovator"))
ggsave(filename="manuscript/images/appendix1_c.pdf",plot=appx1_c,height=5*(3.4/7),width=6,units="in",scale=1.5)

appx1_d <- f.appx1(plot.data,.location=c("Argentina"),.variety=c("Markies Russet"))
ggsave(filename="manuscript/images/appendix1_d.pdf",plot=appx1_d,height=5*(2.5/7),width=6,units="in",scale=1.5)

appx1_e <- f.appx1(plot.data,.location=c("Argentina"),.variety=c("Umatilla Russet"))
ggsave(filename="manuscript/images/appendix1_e.pdf",plot=appx1_e,height=5*(2.5/7),width=6,units="in",scale=1.5)

appx1_f <- f.appx1(plot.data,.location=c("Belgium"),.variety=c("Bintje"))
ggsave(filename="manuscript/images/appendix1_f.pdf",plot=appx1_f,height=5*(7/7),width=6,units="in",scale=1.5)

appx1_g <- f.appx1(plot.data,.location=c("Belgium"),.variety=c("Charlotte"))
ggsave(filename="manuscript/images/appendix1_g.pdf",plot=appx1_g,height=5*(3.4/7),width=6,units="in",scale=1.5)

appx1_h <- f.appx1(plot.data,.location=c("Canada"),.variety=c("Russet Burbank"))
ggsave(filename="manuscript/images/appendix1_h.pdf",plot=appx1_h,height=5*(4.3/7),width=6,units="in",scale=1.5)

appx1_i <- f.appx1(plot.data,.location=c("Canada"),.variety=c("Shepody"))
ggsave(filename="manuscript/images/appendix1_i.pdf",plot=appx1_i,height=5*(4.3/7),width=6,units="in",scale=1.5)

appx1_j <- f.appx1(plot.data,.location=c("Minnesota"),.variety=c("Clearwater"))
ggsave(filename="manuscript/images/appendix1_j.pdf",plot=appx1_j,height=5*(2.5/7),width=6,units="in",scale=1.5)

appx1_k <- f.appx1(plot.data,.location=c("Minnesota"),.variety=c("Dakota Russet"))
ggsave(filename="manuscript/images/appendix1_k.pdf",plot=appx1_k,height=5*(2.5/7),width=6,units="in",scale=1.5)

appx1_l <- f.appx1(plot.data,.location=c("Minnesota"),.variety=c("Easton"))
ggsave(filename="manuscript/images/appendix1_l.pdf",plot=appx1_l,height=5*(2.5/7),width=6,units="in",scale=1.5)

appx1_m <- f.appx1(plot.data,.location=c("Minnesota"),.variety=c("Russet Burbank"))
ggsave(filename="manuscript/images/appendix1_m.pdf",plot=appx1_m,height=5*(7/7),width=6,units="in",scale=1.5)

appx1_n <- f.appx1(plot.data,.location=c("Minnesota"),.variety=c("Umatilla"))
ggsave(filename="manuscript/images/appendix1_n.pdf",plot=appx1_n,height=5*(2.5/7),width=6,units="in",scale=1.5)

# appendix 1 - plateau model fit with point data for each date shown for each variety x location ------------------

# .location = "Minnesota"
# .variety = "Russet Burbank"

f.appx1 <- function(plot.data,.location,.variety){
  
  c <- plot.data$c %>%
    filter(location%in%.location) %>%
    filter(variety%in%.variety) %>%
    mutate_at(vars(location,variety,`location:variety`), ~fct_drop(.))
  
  d <- plot.data$d %>%
    filter(location%in%.location) %>%
    filter(variety%in%.variety) %>%
    mutate_at(vars(location,variety,`location:variety`), ~fct_drop(.))
  
  p <- plot.data$p %>%
    filter(location%in%.location) %>%
    filter(variety%in%.variety) %>%
    mutate_at(vars(location,variety,`location:variety`), ~fct_drop(.))
  
  ggplot() +
    geom_line(data=c,aes(x=W,y=N_0.5),linetype=1,alpha=1.0) +
    # geom_line(data=c,aes(x=W,y=N_0.05),linetype=2,alpha=0.5) +
    # geom_line(data=c,aes(x=W,y=N_0.95),linetype=2,alpha=0.5) +
    geom_line(data=p,aes(x=W_0.5,y=N_0.5,group=index),linetype=1,alpha=0.5) +
    # geom_line(data=p,aes(x=W_0.05,y=N_0.05,group=index),linetype=2,alpha=0.3) +
    # geom_line(data=p,aes(x=W_0.95,y=N_0.95,group=index),linetype=2,alpha=0.3) +
    geom_point(data=d,aes(x=W,y=N),alpha=0.33) +
    # facet_wrap(vars(`location:variety`,index)) +
    facet_wrap(vars(as.numeric(index)),ncol=8) +
    # facet_grid(`location:variety`~index) +
    labs(x = "W",
         y = "%N",
         title = paste(.location,.variety,sep=" - ")) + 
    coord_cartesian(xlim=c(0,NA),ylim=c(0,6.0)) +
    theme_classic()
  
}

appx1_a <- f.appx1(plot.data,.location=c("Argentina"),.variety=c("Bannock Russet"))
ggsave(filename="manuscript/images/appendix1_a.pdf",plot=appx1_a,height=5*(2.5/7),width=6,units="in",scale=1.5)

appx1_b <- f.appx1(plot.data,.location=c("Argentina"),.variety=c("Gem Russet"))
ggsave(filename="manuscript/images/appendix1_b.pdf",plot=appx1_b,height=5*(3.4/7),width=6,units="in",scale=1.5)

appx1_c <- f.appx1(plot.data,.location=c("Argentina"),.variety=c("Innovator"))
ggsave(filename="manuscript/images/appendix1_c.pdf",plot=appx1_c,height=5*(3.4/7),width=6,units="in",scale=1.5)

appx1_d <- f.appx1(plot.data,.location=c("Argentina"),.variety=c("Markies Russet"))
ggsave(filename="manuscript/images/appendix1_d.pdf",plot=appx1_d,height=5*(2.5/7),width=6,units="in",scale=1.5)

appx1_e <- f.appx1(plot.data,.location=c("Argentina"),.variety=c("Umatilla Russet"))
ggsave(filename="manuscript/images/appendix1_e.pdf",plot=appx1_e,height=5*(2.5/7),width=6,units="in",scale=1.5)

appx1_f <- f.appx1(plot.data,.location=c("Belgium"),.variety=c("Bintje"))
ggsave(filename="manuscript/images/appendix1_f.pdf",plot=appx1_f,height=5*(7/7),width=6,units="in",scale=1.5)

appx1_g <- f.appx1(plot.data,.location=c("Belgium"),.variety=c("Charlotte"))
ggsave(filename="manuscript/images/appendix1_g.pdf",plot=appx1_g,height=5*(3.4/7),width=6,units="in",scale=1.5)

appx1_h <- f.appx1(plot.data,.location=c("Canada"),.variety=c("Russet Burbank"))
ggsave(filename="manuscript/images/appendix1_h.pdf",plot=appx1_h,height=5*(4.3/7),width=6,units="in",scale=1.5)

appx1_i <- f.appx1(plot.data,.location=c("Canada"),.variety=c("Shepody"))
ggsave(filename="manuscript/images/appendix1_i.pdf",plot=appx1_i,height=5*(4.3/7),width=6,units="in",scale=1.5)

appx1_j <- f.appx1(plot.data,.location=c("Minnesota"),.variety=c("Clearwater"))
ggsave(filename="manuscript/images/appendix1_j.pdf",plot=appx1_j,height=5*(2.5/7),width=6,units="in",scale=1.5)

appx1_k <- f.appx1(plot.data,.location=c("Minnesota"),.variety=c("Dakota Russet"))
ggsave(filename="manuscript/images/appendix1_k.pdf",plot=appx1_k,height=5*(2.5/7),width=6,units="in",scale=1.5)

appx1_l <- f.appx1(plot.data,.location=c("Minnesota"),.variety=c("Easton"))
ggsave(filename="manuscript/images/appendix1_l.pdf",plot=appx1_l,height=5*(2.5/7),width=6,units="in",scale=1.5)

appx1_m <- f.appx1(plot.data,.location=c("Minnesota"),.variety=c("Russet Burbank"))
ggsave(filename="manuscript/images/appendix1_m.pdf",plot=appx1_m,height=5*(7/7),width=6,units="in",scale=1.5)

appx1_n <- f.appx1(plot.data,.location=c("Minnesota"),.variety=c("Umatilla"))
ggsave(filename="manuscript/images/appendix1_n.pdf",plot=appx1_n,height=5*(2.5/7),width=6,units="in",scale=1.5)


# END --------------------