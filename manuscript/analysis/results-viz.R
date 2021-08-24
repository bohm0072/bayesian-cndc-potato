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
data_cndc <- read_csv("data/analysis/data_cndc.csv",col_types="cccccccdcdd")

# read in model fit results ------------------

model <- readRDS("manuscript/models/model_070621.rds"); model

# the fmin() function used in Stan isn't defined in R, so we need to create it so that when we try to use brms to make predictions, it knows what to do with the fmin()
fmin <- function(x,y){
  pmin(x,y)
}


# format data for results and figures ---------------

f.plot.data <- function(data,model){
  
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
        arrange(variety,location) %>%
        mutate_at(vars(variety), ~as_factor(.)) %>%
        mutate_at(vars(variety), ~fct_inorder(.)) %>%
        arrange(location,variety) %>%
        mutate_at(vars(location), ~as_factor(.)) %>%
        mutate_at(vars(location), ~fct_inorder(.))
      
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
  
  # read in previous critical N dilution curve fits
  f.parm.orig.sum <- function(){
    
    tibble(
      location=c("Argentina","Argentina","Argentina","Argentina","Argentina","Belgium","Belgium","Canada","Canada","Minnesota","Minnesota","Minnesota","Minnesota","Minnesota"),
      variety=c("Bannock Russet","Gem Russet","Innovator","Markies Russet","Umatilla Russet","Bintje","Charlotte","Russet Burbank","Shepody","Clearwater","Dakota Russet","Easton","Russet Burbank","Umatilla Russet"),
      alpha1_orig=c(5.30,5.32,5.30,5.53,5.19,5.37,5.37,4.57,5.04,5.00,5.00,5.00,5.00,5.00),
      alpha2_orig=c(0.25,0.36,0.42,0.25,0.25,0.45,0.45,0.42,0.42,0.45,0.45,0.45,0.45,0.45)
    ) %>%
      mutate(variety.name=str_replace(variety," ",".")) %>%
      mutate(`location:variety`=paste(location,"_",variety.name,sep="")) %>%
      select(-variety.name) %>%
      relocate(`location:variety`,.after=variety) %>%
      arrange(variety,location) %>%
      mutate_at(vars(variety), ~as_factor(.)) %>%
      mutate_at(vars(variety), ~fct_inorder(.)) %>%
      arrange(location,variety) %>%
      mutate_at(vars(location,`location:variety`), ~as_factor(.)) %>%
      mutate_at(vars(location,`location:variety`), ~fct_inorder(.))
    
  }
  parm.orig.sum <- f.parm.orig.sum()
  
  f.parm.fit.sum <- function(cndc.fit){
    
    parm.fit.sum <- left_join(
      cndc.fit %>%
        group_by(location,variety,`location:variety`) %>%
        summarize(qs = quantile(`location:variety_alpha1`,c(0.05,0.50,0.95)), prob = c(0.05,0.50,0.95), .groups="drop") %>%
        pivot_wider(names_from=prob,
                    names_prefix="alpha1_",
                    values_from=qs) %>%
        mutate_at(vars(alpha1_0.05,alpha1_0.5,alpha1_0.95),as.numeric),
      cndc.fit %>%
        group_by(location,variety,`location:variety`) %>%
        summarize(qs = quantile(`location:variety_alpha2`,c(0.05,0.50,0.95)), prob = c(0.05,0.50,0.95), .groups="drop") %>%
        pivot_wider(names_from=prob,
                    names_prefix="alpha2_",
                    values_from=qs)%>%
        mutate_at(vars(alpha2_0.05,alpha2_0.5,alpha2_0.95),as.numeric),
      by = c("location","variety","location:variety")
    ) %>%
      arrange(variety,location) %>%
      mutate_at(vars(variety), ~as_factor(.)) %>%
      mutate_at(vars(variety), ~fct_inorder(.)) %>%
      arrange(location,variety) %>%
      mutate_at(vars(location,`location:variety`), ~as_factor(.)) %>%
      mutate_at(vars(location,`location:variety`), ~fct_inorder(.))
    
    return(parm.fit.sum)
    
  }
  parm.fit.sum <- f.parm.fit.sum(cndc.fit)
  
  f.parm.fit.sum2 <- function(cndc.fit){
    
    parm.fit.sum2 <- left_join(
      cndc.fit %>%
        group_by(location) %>%
        summarize(qs = quantile(`location_alpha1`,c(0.05,0.50,0.95)), prob = c(0.05,0.50,0.95), .groups="drop") %>%
        pivot_wider(names_from=prob,
                    names_prefix="alpha1_",
                    values_from=qs) %>%
        mutate_at(vars(alpha1_0.05,alpha1_0.5,alpha1_0.95),as.numeric),
      cndc.fit %>%
        group_by(location) %>%
        summarize(qs = quantile(`location_alpha2`,c(0.05,0.50,0.95)), prob = c(0.05,0.50,0.95), .groups="drop") %>%
        pivot_wider(names_from=prob,
                    names_prefix="alpha2_",
                    values_from=qs)%>%
        mutate_at(vars(alpha2_0.05,alpha2_0.5,alpha2_0.95),as.numeric),
      by = c("location")
    ) %>%
      arrange(location) %>%
      mutate_at(vars(location), ~as_factor(.)) %>%
      mutate_at(vars(location), ~fct_inorder(.))
    
    return(parm.fit.sum2)
    
  }
  parm.fit.sum2 <- f.parm.fit.sum2(cndc.fit)
  
  f.cndc.draws <- function(cndc.fit){
    
    cndc.draws <- cndc.fit %>%
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
      mutate(`location_N`=`location_alpha1`*(W^(-`location_alpha2`)))
    
    return(cndc.draws)
    
  }
  cndc.draws <- f.cndc.draws(cndc.fit)
  
  f.cndc.fit.sum <- function(cndc.draws){
    
    cndc.fit.sum <- cndc.draws %>%
      group_by(location,variety,`location:variety`,W) %>%
      summarize(qs = quantile(`location:variety_N`,c(0.05,0.50,0.95)), prob = c(0.05,0.50,0.95), .groups="drop") %>%
      pivot_wider(names_from=prob,
                  names_prefix="N_",#"N_location:variety_",
                  values_from=qs) %>%
      mutate_at(vars(N_0.05,N_0.5,N_0.95),as.numeric) %>%
      # mutate_at(vars(`N_location:variety_0.05`,`N_location:variety_0.5`,`N_location:variety_0.95`),as.numeric) %>%
      arrange(variety,location) %>%
      mutate_at(vars(variety), ~as_factor(.)) %>%
      mutate_at(vars(variety), ~fct_inorder(.)) %>%
      arrange(location,variety) %>%
      mutate_at(vars(location,`location:variety`), ~as_factor(.)) %>%
      mutate_at(vars(location,`location:variety`), ~fct_inorder(.))
    
    return(cndc.fit.sum)
    
  }
  cndc.fit.sum <- f.cndc.fit.sum(cndc.draws)
  
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
          mutate_at(vars(location,`location:variety`),as.character),
        by=c(".chain", ".iteration", ".draw", "location", "location:variety"="location:variety")) %>%
      mutate(Nc=`location:variety_alpha1`*(Bmax^(-`location:variety_alpha2`)))
    
    d2 <- d1 %>%
      relocate(.chain,.iteration,.draw,.before=index) %>%
      relocate(`location:variety_alpha1`,`location:variety_alpha2`,.before=Bmax) %>%
      arrange(variety,location) %>%
      mutate_at(vars(variety), ~as_factor(.)) %>%
      mutate_at(vars(variety), ~fct_inorder(.)) %>%
      arrange(location,variety) %>%
      mutate_at(vars(location), ~as_factor(.)) %>%
      mutate_at(vars(location), ~fct_inorder(.))
    
    d <- d2
    
    return(d)
    
  }
  plateau.fit <- f.plateau.fit(model,data)
  
  f.plateau.fit.sum <- function(plateau.fit){
    
    #q=0.5
    f.plateau.fit.quantile <- function(q){
      
      var1 <- paste("W",q,sep="_")
      var2 <- paste("N",q,sep="_")
      
      bind_rows(
        plateau.fit %>%
          mutate(W=Bmax,
                 N=Nc) %>%
          group_by(location,variety,`location:variety`,index) %>%
          summarize_at(vars(W,N),~quantile(.,c(q),na.rm=T), prob = c(q), .groups="drop"),
        # summarize_at(vars(W,N),~quantile(.,c(0.50)), prob = c(0.50), .groups="drop"),
        plateau.fit %>%
          mutate(W=Bmax,
                 N=7.0) %>%
          group_by(location,variety,`location:variety`,index) %>%
          summarize_at(vars(W,N),~quantile(.,c(q),na.rm=T), prob = c(q), .groups="drop"),
        plateau.fit %>%
          mutate(W=Bmax + Si * (0-Nc),
                 N=0) %>%
          group_by(location,variety,`location:variety`,index) %>%
          summarize_at(vars(W,N),~quantile(.,c(q),na.rm=T), prob = c(q), .groups="drop")
      ) %>%
        ungroup() %>%
        mutate_at(vars(W,N),as.numeric) %>%
        rename(!!sym(var1):=W,!!sym(var2):=N) %>%
        # rename(W_0.5=W,N_0.5=N) %>%
        arrange(variety,location) %>%
        mutate_at(vars(variety), ~as_factor(.)) %>%
        mutate_at(vars(variety), ~fct_inorder(.)) %>%
        arrange(location,variety) %>%
        mutate_at(vars(location,`location:variety`), ~as_factor(.)) %>%
        mutate_at(vars(location,`location:variety`), ~fct_inorder(.))
      
    }
    
    plateau.fit.sum <- f.plateau.fit.quantile(0.5) %>%
      left_join(f.plateau.fit.quantile(0.05), c("location", "variety", "location:variety", "index")) %>%
      left_join(f.plateau.fit.quantile(0.95), c("location", "variety", "location:variety", "index"))
    
    return(plateau.fit.sum)
    
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
      arrange(variety,location) %>%
      mutate_at(vars(variety), ~as_factor(.)) %>%
      mutate_at(vars(variety), ~fct_inorder(.)) %>%
      arrange(location,variety) %>%
      mutate_at(vars(location,`location:variety`), ~as_factor(.)) %>%
      mutate_at(vars(location,`location:variety`), ~fct_inorder(.))
    
    return(Bmax.sum)
    
  }
  Bmax.sum <- f.Bmax.sum(plateau.fit)
  
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
    arrange(variety,location) %>%
    mutate_at(vars(variety), ~as_factor(.)) %>%
    mutate_at(vars(variety), ~fct_inorder(.)) %>%
    arrange(location,variety) %>%
    mutate_at(vars(location,`location:variety`), ~as_factor(.)) %>%
    mutate_at(vars(location,`location:variety`), ~fct_inorder(.))
  
  # p <- plateau.fit.sum 
   
  r <- cndc.draws %>%
    select(.draw,location,variety,`location:variety`,W,`location:variety_N`) %>%
    rename(N=`location:variety_N`) %>%
    left_join(Bmax.sum %>%
                select(location,variety,`location:variety`,Bmax_0.95),
              by=c("location","variety","location:variety")) %>%
    filter(W <= Bmax_0.95) %>%
    select(-Bmax_0.95)
  
  out <- list(c=c,
              d=d,
              p=plateau.fit.sum,
              r=r,
              cndc.fit=cndc.fit,
              parm.fit.sum=parm.fit.sum,
              parm.fit.sum2=parm.fit.sum2,
              parm.orig.sum=parm.orig.sum)
  
  return(out)
  
}
plot.data <- f.plot.data(data_cndc,model)

#Free Unsued Memory
gc()

# set color scale -----------------

plot.colors.1 <- c("#490809","#931012","#CA1619","#C55962","#E74B5E","#2C652A","#33a02c","#CC6600","#ff7f00","#15527A","#0177A2","#1F78B4","#01B3F4","#62B0E4")
plot.colors.2 <- c("#CA1619","#33a02c","#ff7f00","#1F78B4")

# figure 3 table - distribution of alpha parameter values for each parameter independently #####

plot.data$tab1 <- bind_rows(plot.data$parm.fit.sum,
                  plot.data$parm.fit.sum2 %>%
                    mutate(variety=".",
                           `location:variety`=".")) %>%
  arrange(location,variety)

# figure 3 - distribution of alpha parameter values for each parameter independently ------------------

# parm = "alpha1"
# .colors1 = plot.colors.1
# .colors2 = plot.colors.2
# alpha1_0.5 <- 5.00
# alpha1_0.05 <- 4.75
# alpha1_0.95 <- 5.25

f.fig1 <- function(plot.data,parm,.colors1,.colors2){
  
  var4 <- paste("location_",parm,sep="")
  var5 <- paste("location:variety_",parm,sep="")
  
  coord.x <- case_when(
    parm == "alpha1" ~ c(4.0,6.75),
    parm == "alpha2" ~ c(0.1,1.25)
  )
  
  breaks.x <- case_when(
    parm == "alpha1" ~ c(4.0,4.5,5.0,5.5),
    parm == "alpha2" ~ c(0.1,0.3,0.5,0.7)
  )
  
  minor.breaks.x <- case_when(
    parm == "alpha1" ~ c(4.25,4.75,5.25,5.5),
    parm == "alpha2" ~ c(0.0,0.2,0.4,0.6)
  )
  
  label.x <- case_when(
    parm == "alpha1" ~ expression(paste("parameter ", italic("a"), sep=" ")),
    parm == "alpha2" ~ expression(paste("parameter ", italic("b"), sep=" "))
  )
  
  parm.pos.x <- case_when(
    parm == "alpha1" ~ 6.75,
    parm == "alpha2" ~ 1.28
  )
  
  parm.lab.x <- case_when(
    parm == "alpha1" ~ "alpha1_lab",
    parm == "alpha2" ~ "alpha2_lab"
  )
  
  t <- plot.data$tab1 %>%
    mutate(alpha1_lab=paste(format(round(alpha1_0.5,2),nsmall=2),' (',format(round(alpha1_0.05,2),nsmall=2),', ',format(round(alpha1_0.95,2),nsmall=2),')',sep=''),
           alpha2_lab=paste(format(round(alpha2_0.5,3),nsmall=3),' (',format(round(alpha2_0.05,3),nsmall=3),', ',format(round(alpha2_0.95,3),nsmall=3),')',sep=''))
  
  p1 <- ggplot() +
    stat_halfeye(data=plot.data$cndc.fit, aes(x = !!sym(var5), y = reorder(variety, desc(variety)), fill=`location:variety`),.width = c(0.1, 0.90),size=2) +
    geom_text(data=filter(t,variety!="."), aes(x=parm.pos.x,y=variety,label=!!sym(parm.lab.x)),nudge_y=0.1,hjust=1,vjust=0,size=3,color="#4d4d4d") + #label="5.00 (4.75,5.25)"
    facet_grid(location~., scales = "free_y", space = "free") +
    coord_cartesian(xlim=coord.x) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    guides(fill = "none") +
    scale_fill_manual(values=.colors1) +
    scale_x_continuous(breaks=breaks.x,
                       minor_breaks=minor.breaks.x)
    
  g1 <- ggplotGrob(p1)
  fg1 <- gtable_frame(g1, height = unit(3, "null")) #debug = TRUE,
  
  p2 <- ggplot() +
    stat_halfeye(data=plot.data$cndc.fit, aes(x = !!sym(var4), y = reorder(location, desc(location)), fill=location),.width = c(0.1, 0.90),size=2) +
    geom_text(data=filter(t,variety=="."), aes(x=parm.pos.x,y=location,label=!!sym(parm.lab.x)),nudge_y=0.1,hjust=1,vjust=0,size=3,color="#4d4d4d") + #label="5.00 (4.75,5.25)"
    coord_cartesian(xlim=coord.x) +
    labs(x=label.x) +
    theme_bw() +
    theme(axis.title.y = element_blank()) +
    guides(fill = "none") +
    scale_fill_manual(values=.colors2) +
    scale_x_continuous(breaks=breaks.x,
                       minor_breaks=minor.breaks.x)
  
  g2 <- ggplotGrob(p2)
  fg2 <- gtable_frame(g2, height = unit(1, "null"))
  
  fg <- rbind(fg1,fg2,size = "first")
  fg$widths <- unit.pmax(fg1$widths, fg2$widths)
  
  return(fg)
  
}

fig1_a <- f.fig1(plot.data,"alpha1",plot.colors.1,plot.colors.2)
# ggsave(filename="manuscript/images/figure1_a.pdf",plot=fig1_a,height=4.5,width=3,units="in",scale=1.3)
# ggsave(filename="manuscript/images/figure1_a.png",plot=fig1_a,height=4.5,width=3,units="in",scale=1.3,dpi=1000)
fig1_b <- f.fig1(plot.data,"alpha2",plot.colors.1,plot.colors.2)
# ggsave(filename="manuscript/images/figure1_b.pdf",plot=fig1_b,height=4.5,width=3,units="in",scale=1.3)
# ggsave(filename="manuscript/images/figure1_b.png",plot=fig1_b,height=4.5,width=3,units="in",scale=1.3,dpi=1000)

fig1.layout <- rbind(c(1,2))

fig1 <- grid.arrange(fig1_a,fig1_b,
                     layout_matrix=fig1.layout)

ggsave(filename="manuscript/images/figure3.pdf",plot=fig1,height=4.5,width=6.0,units="in",scale=1.3,device=cairo_pdf) #height=4.5,width=6,scale=1.3
ggsave(filename="manuscript/images/figure3.png",plot=fig1,height=4.5,width=6,units="in",scale=1.3,dpi=1000)

# figure 4 - distribution of alpha parameters for each parameters simultaneously ----------------

# .location = "Argentina"
# .variety = "Innovator"
# .color = "#e41a1c"

f.fig2 <- function(plot.data,.location,.variety,.color){
  
  var3 <- paste("       ",.variety,sep="")
  
  d <- plot.data$cndc.fit %>%
    filter(location %in% .location) %>%
    filter(variety %in% .variety)
  
  p1 <- ggplot(data = d, aes(x=`location:variety_alpha1`, y=`location:variety_alpha2`, color=`location:variety`)) +
    geom_point(alpha=0.02) +
    stat_cor(color="black",size=2,
             aes(label = ..r.label..),
             label.x=4.0,
             label.y=0.78,
             r.accuracy=0.01,
             hjust=0) +
    theme_classic() +
    theme(text=element_text(size=8),
          plot.title=element_text(size=8),
          axis.title=element_blank(),
          plot.caption = element_text(hjust=0,size=7), 
          plot.title.position = "plot", 
          plot.caption.position =  "plot") +
    scale_color_manual(values = .color) +
    scale_x_continuous(limits=c(4.0,5.5),breaks=c(4.0,4.5,5.0,5.5)) +
    scale_y_continuous(limits=c(0.0,0.8),breaks=c(0.1,0.3,0.5,0.7)) +
    guides(color="none") +
    labs(caption=var3)
  
  p2 <- ggMarginal(p1, type="density", data = plot.data$cndc.fit, groupColour = T)
  
  return(p2)
  
}

fig2.list <- list(
  location=c("Argentina","Argentina","Argentina","Argentina","Argentina","Belgium","Belgium","Canada","Canada","Minnesota","Minnesota","Minnesota","Minnesota","Minnesota"),
  variety=c("Bannock Russet","Gem Russet","Innovator","Markies Russet","Umatilla Russet","Bintje","Charlotte","Russet Burbank","Shepody","Clearwater","Dakota Russet","Easton","Russet Burbank","Umatilla Russet"),
  color=plot.colors.1
)

fig2.sub <- pmap(fig2.list,~f.fig2(plot.data,
                                   .location=..1,
                                   .variety=..2,
                                   .color=..3))

# f.fig2(plot.data,fig2.list$location[[14]],fig2.list$variety[[14]],fig2.list$color[[14]])

# fig2.sub[[1]]

f.fig2.lab.facet <- function(.location){
  
  ggplot() +
    geom_text(aes(x=0,y=0,label=.location),angle=270,size=2.5) +
    theme_classic() +
    # theme_grey() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_rect(fill="#d2d2d2",color="black"))
  
}
f.fig2.lab.axis.y <- function(){
  
  ggplot() +
    geom_text(aes(x=0,y=0),label=expression(paste("parameter ", italic("b"), sep=" ")),parse=T,angle=90,size=3) +
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank())
  
}
f.fig2.lab.axis.x <- function(){
  
  ggplot() +
    geom_text(aes(x=0,y=0),label=expression(paste("parameter ", italic("a"), sep=" ")),parse=T,angle=0,size=3) +
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank())
  
}

fig2.layout <- rbind(c(19,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,15),
                     c(19,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,15),
                     c(19,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,15),
                     c(19,6,6,6,7,7,7,16,NA,NA,8,8,8,9,9,9,17),
                     c(19,6,6,6,7,7,7,16,NA,NA,8,8,8,9,9,9,17),
                     c(19,6,6,6,7,7,7,16,NA,NA,8,8,8,9,9,9,17),
                     c(19,10,10,10,11,11,11,12,12,12,13,13,13,14,14,14,18),
                     c(19,10,10,10,11,11,11,12,12,12,13,13,13,14,14,14,18),
                     c(19,10,10,10,11,11,11,12,12,12,13,13,13,14,14,14,18),
                     c(NA,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20))

fig2 <- grid.arrange(fig2.sub[[1]],fig2.sub[[2]],fig2.sub[[3]],fig2.sub[[4]],fig2.sub[[5]],
                     fig2.sub[[6]],fig2.sub[[7]],fig2.sub[[8]],fig2.sub[[9]],
                     fig2.sub[[10]],fig2.sub[[11]],fig2.sub[[12]],fig2.sub[[13]],fig2.sub[[14]],
                     f.fig2.lab.facet("Argentina"),f.fig2.lab.facet("Belgium"),f.fig2.lab.facet("Canada"),f.fig2.lab.facet("Minnesota"),
                     f.fig2.lab.axis.y(),f.fig2.lab.axis.x(),
                     layout_matrix=fig2.layout)

ggsave(filename="manuscript/images/figure4.pdf",plot=fig2,height=4,width=6,units="in",scale=1,device=cairo_pdf)
ggsave(filename="manuscript/images/figure4.png",plot=fig2,height=4,width=6,units="in",scale=1,dpi=1000)


# figure 5 - curve fits for each variety x location  ------------------

# .location = "Argentina"
# .variety = "Innovator"
# .color = "#e41a1c"

f.fig3 <- function(plot.data,.location,.variety,.color){
  
  # var1 <- "W"
  # var2 <- "%N"
  var3 <- paste("     ",.variety,sep="")
  
  c <- plot.data$c %>%
    filter(location %in% .location) %>%
    filter(variety %in% .variety)
  
  d <- plot.data$d %>%
    filter(location %in% .location) %>%
    filter(variety %in% .variety)
  
  p <- plot.data$p %>%
    filter(location %in% .location) %>%
    filter(variety %in% .variety)
  
  n <- d %>%
    drop_na() %>%
    summarize(n=n()) %>%
    mutate(x=max(c$W)) %>%
    mutate(lab=paste("n=",n,sep=""))
  
  i <- d %>%
    select(index) %>%
    distinct() %>%
    drop_na() %>%
    summarize(i=n()) %>%
    mutate(x=max(c$W)) %>%
    mutate(lab=paste("i=",i,sep=""))
  
  ggplot() +
    geom_line(data=p,aes(x=W_0.5,y=N_0.5,group=index),color="grey",linetype=1,alpha=0.20) + #,color=`location:variety`
    geom_point(data=d,aes(x=W,y=N,color=`location:variety`),alpha=0.20) + 
    geom_line(data=c,aes(x=W,y=N_0.5,group=`location:variety`),linetype=1,alpha=1.0) +
    geom_text(data=n,aes(x=x,y=6.4,label=lab),size=2,hjust=1,vjust=1) +
    geom_text(data=i,aes(x=x,y=5.4,label=lab),size=2,hjust=1,vjust=1) +
    labs(caption=var3) +
    coord_cartesian(xlim=c(0,NA),ylim=c(0,6.5)) +
    theme_classic() +
    theme(text=element_text(size=8),
          plot.title=element_text(size=8),
          axis.title=element_blank(),
          plot.caption = element_text(hjust=0,size=7), 
          plot.title.position = "plot", 
          plot.caption.position =  "plot") +
    guides(color="none") +
    scale_color_manual(values=.color)  +
    scale_x_continuous(n.breaks=4)
   
}

fig3.list <- list(
  location=c("Argentina","Argentina","Argentina","Argentina","Argentina","Belgium","Belgium","Canada","Canada","Minnesota","Minnesota","Minnesota","Minnesota","Minnesota"),
  variety=c("Bannock Russet","Gem Russet","Innovator","Markies Russet","Umatilla Russet","Bintje","Charlotte","Russet Burbank","Shepody","Clearwater","Dakota Russet","Easton","Russet Burbank","Umatilla Russet"),
  color=plot.colors.1
)

fig3.sub <- pmap(fig3.list,~f.fig3(plot.data,
                                   .location=..1,
                                   .variety=..2,
                                   .color=..3))

f.fig3.lab.facet <- function(.location){
  
  ggplot() +
    geom_text(aes(x=0,y=0,label=.location),angle=270,size=2.5) +
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_rect(fill="#d2d2d2",color="black"))
  
}
f.fig3.lab.axis.y <- function(){
  
  ggplot() +
    geom_text(aes(x=0,y=0),label=expression("%N [g N 100 g"^-1*"]"),angle=90,size=3,parse=T) +
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank())
  
}
f.fig3.lab.axis.x <- function(){
  
  ggplot() +
    geom_text(aes(x=0,y=0),label=expression("Biomass [Mg ha"^-1*"]"),angle=0,size=3,parse=T) +
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank())
  
}

fig3.layout <- rbind(c(19,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,15),
                     c(19,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,15),
                     c(19,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,15),
                     c(19,6,6,6,7,7,7,16,NA,NA,8,8,8,9,9,9,17),
                     c(19,6,6,6,7,7,7,16,NA,NA,8,8,8,9,9,9,17),
                     c(19,6,6,6,7,7,7,16,NA,NA,8,8,8,9,9,9,17),
                     c(19,10,10,10,11,11,11,12,12,12,13,13,13,14,14,14,18),
                     c(19,10,10,10,11,11,11,12,12,12,13,13,13,14,14,14,18),
                     c(19,10,10,10,11,11,11,12,12,12,13,13,13,14,14,14,18),
                     c(NA,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20))

fig3 <- grid.arrange(fig3.sub[[1]],fig3.sub[[2]],fig3.sub[[3]],fig3.sub[[4]],fig3.sub[[5]],
                     fig3.sub[[6]],fig3.sub[[7]],fig3.sub[[8]],fig3.sub[[9]],
                     fig3.sub[[10]],fig3.sub[[11]],fig3.sub[[12]],fig3.sub[[13]],fig3.sub[[14]],
                     f.fig3.lab.facet("Argentina"),f.fig3.lab.facet("Belgium"),f.fig3.lab.facet("Canada"),f.fig3.lab.facet("Minnesota"),
                     f.fig3.lab.axis.y(),f.fig3.lab.axis.x(),
                     layout_matrix=fig3.layout)

ggsave(filename="manuscript/images/figure5.pdf",plot=fig3,height=4,width=6,units="in",scale=1.0,device=cairo_pdf)
ggsave(filename="manuscript/images/figure5.png",plot=fig3,height=4,width=6,units="in",scale=1.0,dpi=1000)

# table 6 - fitted boundary nls curves -----------------

# .location = "Minnesota"
# .variety = "Russet Burbank"

f.tab2 <- function(plot.data,.location,.variety){
  
  t1 <- plot.data$tab1 %>%
    filter(variety!=".") %>%
    filter(location %in% .location) %>%
    filter(variety %in% .variety)
  
  c <- plot.data$c %>%
    filter(location %in% .location) %>%
    filter(variety %in% .variety)
  
  nls_0.95 <- nls(data=c,
                  formula=N_0.95~alpha1*W^(-alpha2),
                  start=list(alpha1=4.94,
                             alpha2=0.40))
  
  nls_0.5 <- nls(data=c,
                 formula=N_0.5~alpha1*W^(-alpha2),
                 start=list(alpha1=4.94,
                            alpha2=0.40))
  
  nls_0.05 <- nls(data=c,
                  formula=N_0.05~alpha1*W^(-alpha2),
                  start=list(alpha1=4.94,
                             alpha2=0.40))
  
  nls_alpha1_0.05 <- summary(nls_0.05)$parameters[1,1]
  nls_alpha2_0.05 <- summary(nls_0.05)$parameters[2,1]
  nls_alpha1_0.5 <- summary(nls_0.5)$parameters[1,1]
  nls_alpha2_0.5 <- summary(nls_0.5)$parameters[2,1]
  nls_alpha1_0.95 <- summary(nls_0.95)$parameters[1,1]
  nls_alpha2_0.95 <- summary(nls_0.95)$parameters[2,1]
  
  t <- t1 %>%
    mutate(nls_alpha1_0.05=nls_alpha1_0.05,
           nls_alpha2_0.05=nls_alpha2_0.05,
           nls_alpha1_0.5=nls_alpha1_0.5,
           nls_alpha2_0.5=nls_alpha2_0.5,
           nls_alpha1_0.95=nls_alpha1_0.95,
           nls_alpha2_0.95=nls_alpha2_0.95)
  
  # a1 = "nls_alpha1_0.05"
  # a2 = "nls_alpha2_0.05"
  # name = "nls_est_lo"
  
  f.t2.select <- function(t,a1,a2,name){
    
    # .t <- t %>%
    #   select(!!sym(a1),!!sym(a2)) %>%
    #   pivot_longer(cols=c(!!sym(a1),!!sym(a2)),values_to=name) %>%
    #   select(-name)
    
    # .t <- bind_rows(
    #   .t %>%
    #     slice(1) %>%
    #     mutate_all(~format(round(.,2),nsmall=2)),
    #   .t %>%
    #     slice(2) %>%
    #     mutate_all(~format(round(.,3),nsmall=3))
    # )
    
    .t <- t %>%
      select(!!sym(a1),!!sym(a2))

    .t <- .t %>%
      mutate(across(.cols=c(!!sym(a1)),.fns=~format(round(.,2),nsmall=2))) %>%
      mutate(across(.cols=c(!!sym(a2)),.fns=~format(round(.,3),nsmall=3)))
    
    return(.t)  
    
  }
  
  # t2.list <- list(
  #   a1 = c("alpha1_0.05","nls_alpha1_0.05","alpha1_0.5","nls_alpha1_0.95","alpha1_0.95"),
  #   a2 = c("alpha2_0.95","nls_alpha2_0.05","alpha2_0.5","nls_alpha2_0.95","alpha2_0.05"),
  #   name = c("conserv_est_lo","nls_est_lo","median","nls_est_hi","conserv_est_hi")
  # )
  
  t2.list <- list(
    a1 = c("nls_alpha1_0.05","alpha1_0.5","nls_alpha1_0.95"),
    a2 = c("nls_alpha2_0.05","alpha2_0.5","nls_alpha2_0.95"),
    name = c("nls_est_lo","median","nls_est_hi")
  )
  
  t2 <- pmap(t2.list,~f.t2.select(t,
                                  a1=..1,
                                  a2=..2,
                                  name=..3)) %>% bind_cols()
  
  t2 <- bind_cols(
    c %>%
      select(location,variety,`location:variety`) %>%
      distinct(),
    t2)
  
  return(t2)
  
}

tab2.list <- list(
  location=c("Argentina","Argentina","Argentina","Argentina","Argentina","Belgium","Belgium","Canada","Canada","Minnesota","Minnesota","Minnesota","Minnesota","Minnesota"),
  variety=c("Bannock Russet","Gem Russet","Innovator","Markies Russet","Umatilla Russet","Bintje","Charlotte","Russet Burbank","Shepody","Clearwater","Dakota Russet","Easton","Russet Burbank","Umatilla Russet")
)

tab2 <- pmap(tab2.list,~f.tab2(plot.data,
                               .location=..1,
                               .variety=..2)) %>% bind_rows()

write_csv(tab2,"manuscript/tables/table6.csv")

# figure 6a - curve fits for each variety x location with uncertainty  ------------------

# .location = "Argentina"
# .variety = "Innovator"

f.fig4 <- function(plot.data,.location,.variety){
  
  var1 <- "W"
  # var2 <- "%N"
  var3 <- paste(.location,.variety,sep=" - ")
  var4 <- "%N Difference"
  var5 <- paste("     ",.variety,sep="")
  
  c <- plot.data$c %>%
    filter(location %in% .location) %>%
    filter(variety %in% .variety)
  
  parm <- plot.data$parm.fit.sum %>%
    filter(location %in% .location) %>%
    filter(variety %in% .variety)
  
  nls_0.95 <- nls(data=c,
                  formula=N_0.95~alpha1*W^(-alpha2),
                  start=list(alpha1=4.94,
                             alpha2=0.40))
  
  nls_0.5 <- nls(data=c,
                 formula=N_0.5~alpha1*W^(-alpha2),
                 start=list(alpha1=4.94,
                            alpha2=0.40))
  
  nls_0.05 <- nls(data=c,
                  formula=N_0.05~alpha1*W^(-alpha2),
                  start=list(alpha1=4.94,
                             alpha2=0.40))
  
  nls_alpha1_0.05 <- summary(nls_0.05)$parameters[1,1]
  nls_alpha2_0.05 <- summary(nls_0.05)$parameters[2,1]
  nls_alpha1_0.5 <- summary(nls_0.5)$parameters[1,1]
  nls_alpha2_0.5 <- summary(nls_0.5)$parameters[2,1]
  nls_alpha1_0.95 <- summary(nls_0.95)$parameters[1,1]
  nls_alpha2_0.95 <- summary(nls_0.95)$parameters[2,1]
  
  c <- c %>%
    rowwise() %>%
    mutate(N_0.05_nls = nls_alpha1_0.05*W^(-nls_alpha2_0.05),
           N_0.5_nls = nls_alpha1_0.5*W^(-nls_alpha2_0.5),
           N_0.95_nls = nls_alpha1_0.95*W^(-nls_alpha2_0.95)) %>%
    ungroup()
  
  c <- c %>%
    left_join(parm, by = c("location","variety","location:variety")) %>%
    rowwise() %>%
    mutate(N_0.05_est = alpha1_0.05*W^(-alpha2_0.95),
           N_0.5_est = alpha1_0.5*W^(-alpha2_0.5),
           N_0.95_est = alpha1_0.95*W^(-alpha2_0.05)) %>%
    ungroup() %>%
    select(-c(alpha1_0.05,alpha1_0.5,alpha1_0.95,alpha2_0.05,alpha2_0.5,alpha2_0.95))
  
  c <- c %>%
    mutate(N_cred_lo=N_0.05-N_0.5,
           N_cred_up=N_0.95-N_0.5,
           N_cred_nls_lo=N_0.05_nls-N_0.5,
           N_cred_nls_up=N_0.95_nls-N_0.5,
           N_cred_est_lo=N_0.05_est-N_0.5,
           N_cred_est_up=N_0.95_est-N_0.5)
  
  l <- c %>%
    summarize_all(last) %>%
    mutate(lab_median=paste("(",format(round(parm$alpha1_0.5,2),nsmall=2),", ",format(round(parm$alpha2_0.5,3),nsmall=3),")",sep=""),
           lab_cred_nls_lo=paste("(",format(round(nls_alpha1_0.05,2),nsmall=2),", ",format(round(nls_alpha2_0.05,3),nsmall=3),")",sep=""),
           lab_cred_nls_up=paste("(",format(round(nls_alpha1_0.95,2),nsmall=2),", ",format(round(nls_alpha2_0.95,3),nsmall=3),")",sep=""),
           lab_cred_est_lo=paste("(",format(round(parm$alpha1_0.05,2),nsmall=2),", ",format(round(parm$alpha2_0.95,3),nsmall=3),")",sep=""),
           lab_cred_est_up=paste("(",format(round(parm$alpha1_0.95,2),nsmall=2),", ",format(round(parm$alpha2_0.05,3),nsmall=3),")",sep=""))
  
  ggplot() +
    geom_ribbon(data=c,aes(x=W,ymin=N_cred_lo,ymax=N_cred_up),alpha=0.20) + #,fill="#737373"
    geom_line(data=c,aes(x=W,y=0,group=`location:variety`),linetype=1,alpha=1.0) +
    geom_line(data=c,aes(x=W,y=N_cred_nls_lo,group=`location:variety`),linetype=3,alpha=1.0,size=0.2) +
    geom_line(data=c,aes(x=W,y=N_cred_nls_up,group=`location:variety`),linetype=3,alpha=1.0,size=0.2) +
    geom_line(data=c,aes(x=W,y=N_cred_est_lo,group=`location:variety`),linetype=2,alpha=1.0,size=0.2) +
    geom_line(data=c,aes(x=W,y=N_cred_est_up,group=`location:variety`),linetype=2,alpha=1.0,size=0.2) +
    theme_classic() +
    theme(text=element_text(size=8),
          plot.title=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          plot.caption = element_text(hjust=0,size=7), 
          plot.title.position = "plot", 
          plot.caption.position =  "plot") +
    labs(x=var1,
         y=var4,
         title=var3,
         caption=var5) +
    scale_y_continuous(limits=c(-0.7,0.7),breaks=c(-0.6,-0.3,0,0.3,0.6)) +
    scale_x_continuous(limits=c(0,NA),n.breaks=4)
  
}

fig4.list <- list(
  location=c("Argentina","Argentina","Argentina","Argentina","Argentina","Belgium","Belgium","Canada","Canada","Minnesota","Minnesota","Minnesota","Minnesota","Minnesota"),
  variety=c("Bannock Russet","Gem Russet","Innovator","Markies Russet","Umatilla Russet","Bintje","Charlotte","Russet Burbank","Shepody","Clearwater","Dakota Russet","Easton","Russet Burbank","Umatilla Russet")
)

fig4.sub <- pmap(fig4.list,~f.fig4(plot.data,
                                   .location=..1,
                                   .variety=..2))

f.fig4.lab.facet <- function(.location){
  
  ggplot() +
    geom_text(aes(x=0,y=0,label=.location),angle=270,size=2.5) +
    theme_classic() +
    # theme_grey() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_rect(fill="#d2d2d2",color="black"))
  
}
f.fig4.lab.axis.y <- function(){
  
  ggplot() +
    geom_text(aes(x=0,y=0),label=expression("∆ %N"[c]*" [g N 100 g"^-1*"]"),angle=90,size=3,parse=T) + #expression("%N Difference [g N 100 g"^-1*"]")
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank())
  
}
f.fig4.lab.axis.x <- function(){
  
  ggplot() +
    geom_text(aes(x=0,y=0),label=expression("Biomass [Mg ha"^-1*"]"),angle=0,size=3,parse=T) +
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank())
  
}

fig4.layout <- rbind(c(19,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,15),
                     c(19,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,15),
                     c(19,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,15),
                     c(19,6,6,6,7,7,7,16,NA,NA,8,8,8,9,9,9,17),
                     c(19,6,6,6,7,7,7,16,NA,NA,8,8,8,9,9,9,17),
                     c(19,6,6,6,7,7,7,16,NA,NA,8,8,8,9,9,9,17),
                     c(19,10,10,10,11,11,11,12,12,12,13,13,13,14,14,14,18),
                     c(19,10,10,10,11,11,11,12,12,12,13,13,13,14,14,14,18),
                     c(19,10,10,10,11,11,11,12,12,12,13,13,13,14,14,14,18),
                     c(NA,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20))

fig4 <- grid.arrange(fig4.sub[[1]],fig4.sub[[2]],fig4.sub[[3]],fig4.sub[[4]],fig4.sub[[5]],
                     fig4.sub[[6]],fig4.sub[[7]],fig4.sub[[8]],fig4.sub[[9]],
                     fig4.sub[[10]],fig4.sub[[11]],fig4.sub[[12]],fig4.sub[[13]],fig4.sub[[14]],
                     f.fig4.lab.facet("Argentina"),f.fig4.lab.facet("Belgium"),f.fig4.lab.facet("Canada"),f.fig4.lab.facet("Minnesota"),
                     f.fig4.lab.axis.y(),f.fig4.lab.axis.x(),
                     layout_matrix=fig4.layout)

ggsave(filename="manuscript/images/figure6a.pdf",plot=fig4,height=4,width=6,units="in",scale=1.0,device=cairo_pdf)
ggsave(filename="manuscript/images/figure6a.png",plot=fig4,height=4,width=6,units="in",scale=1.0,dpi=1000)

# figure 6b - curve fits with uncertainty and spaghetti plot for single variety x location  ------------------

# .location = "Minnesota"
# .variety = "Russet Burbank"

f.fig4.callout <- function(plot.data,.location,.variety,.plot){
  
  var1 <- "W"
  # var2 <- "%N"
  var3 <- paste(.location,.variety,sep=" - ")
  var4 <- "%N Difference"
  var5 <- paste("     ",.variety,sep="")
  
  set.seed(999999)
  
  draw.list <- sample(x = c(1:(plot.data$r %>%select(.draw) %>% distinct() %>% count() %>% pull())),
                      size = 15,
                      replace = F)
  
  r <- plot.data$r %>%
    filter(location %in% .location) %>%
    filter(variety %in% .variety) %>%
    filter(.draw %in% draw.list)
  
  c <- plot.data$c %>%
    filter(location %in% .location) %>%
    filter(variety %in% .variety)
  
  parm <- plot.data$parm.fit.sum %>%
    filter(location %in% .location) %>%
    filter(variety %in% .variety)
  
  nls_0.95 <- nls(data=c,
                  formula=N_0.95~alpha1*W^(-alpha2),
                  start=list(alpha1=4.94,
                             alpha2=0.40))
  
  nls_0.5 <- nls(data=c,
                 formula=N_0.5~alpha1*W^(-alpha2),
                 start=list(alpha1=4.94,
                            alpha2=0.40))
  
  nls_0.05 <- nls(data=c,
                  formula=N_0.05~alpha1*W^(-alpha2),
                  start=list(alpha1=4.94,
                             alpha2=0.40))
  
  nls_alpha1_0.05 <- summary(nls_0.05)$parameters[1,1]
  nls_alpha2_0.05 <- summary(nls_0.05)$parameters[2,1]
  nls_alpha1_0.5 <- summary(nls_0.5)$parameters[1,1]
  nls_alpha2_0.5 <- summary(nls_0.5)$parameters[2,1]
  nls_alpha1_0.95 <- summary(nls_0.95)$parameters[1,1]
  nls_alpha2_0.95 <- summary(nls_0.95)$parameters[2,1]
  
  c <- c %>%
    rowwise() %>%
    mutate(N_0.05_nls = nls_alpha1_0.05*W^(-nls_alpha2_0.05),
           N_0.5_nls = nls_alpha1_0.5*W^(-nls_alpha2_0.5),
           N_0.95_nls = nls_alpha1_0.95*W^(-nls_alpha2_0.95)) %>%
    ungroup()
  
  c <- c %>%
    left_join(parm, by = c("location","variety","location:variety")) %>%
    rowwise() %>%
    mutate(N_0.05_est = alpha1_0.05*W^(-alpha2_0.95),
           N_0.5_est = alpha1_0.5*W^(-alpha2_0.5),
           N_0.95_est = alpha1_0.95*W^(-alpha2_0.05)) %>%
    ungroup() %>%
    select(-c(alpha1_0.05,alpha1_0.5,alpha1_0.95,alpha2_0.05,alpha2_0.5,alpha2_0.95))
  
  c <- c %>%
    mutate(N_cred_lo=N_0.05-N_0.5,
           N_cred_up=N_0.95-N_0.5,
           N_cred_nls_lo=N_0.05_nls-N_0.5,
           N_cred_nls_up=N_0.95_nls-N_0.5,
           N_cred_est_lo=N_0.05_est-N_0.5,
           N_cred_est_up=N_0.95_est-N_0.5)
  
  r <- r %>%
    rename(N_draw=N) %>%
    left_join(c %>%
                select(location,variety,`location:variety`,W,N_0.5),
              by = c("location", "variety", "location:variety", "W")) %>%
    mutate(N_draw=N_draw-N_0.5) %>%
    select(-N_0.5)
  
  l <- c %>%
    summarize_all(last) %>%
    mutate(lab_median=paste("(",format(round(parm$alpha1_0.5,2),nsmall=2),", ",format(round(parm$alpha2_0.5,3),nsmall=3),")",sep=""),
           lab_cred_nls_lo=paste("(",format(round(nls_alpha1_0.05,2),nsmall=2),", ",format(round(nls_alpha2_0.05,3),nsmall=3),")",sep=""),
           lab_cred_nls_up=paste("(",format(round(nls_alpha1_0.95,2),nsmall=2),", ",format(round(nls_alpha2_0.95,3),nsmall=3),")",sep=""),
           lab_cred_est_lo=paste("(",format(round(parm$alpha1_0.05,2),nsmall=2),", ",format(round(parm$alpha2_0.95,3),nsmall=3),")",sep=""),
           lab_cred_est_up=paste("(",format(round(parm$alpha1_0.95,2),nsmall=2),", ",format(round(parm$alpha2_0.05,3),nsmall=3),")",sep=""))
  
  g <- ggplot() +
    geom_ribbon(data=c,aes(x=W,ymin=N_cred_lo,ymax=N_cred_up),alpha=0.20) + #,fill="#404040" #,fill="#737373"
    geom_line(data=c,aes(x=W,y=0,group=`location:variety`),linetype=1,alpha=1.0) +
    theme_classic() +
    theme(text=element_text(size=8),
          plot.title=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          plot.caption = element_text(hjust=0,size=7), 
          plot.title.position = "plot", 
          plot.caption.position =  "plot") +
    scale_y_continuous(limits=c(-0.45,0.45),breaks=c(-0.6,-0.3,0,0.3,0.6)) +
    scale_x_continuous(limits=c(0,NA),breaks=c(0,5,10,15))
  
  if(.plot=="draws"){
    g +
      geom_line(data=r,aes(x=W,y=N_draw,group=.draw),alpha=0.80,color="#ca0020",size=0.20) +
      geom_line(data=c,aes(x=W,y=0,group=`location:variety`),linetype=1,alpha=1.0)
  }else if(.plot=="estimate"){
    g +
      geom_line(data=c,aes(x=W,y=N_cred_nls_lo,group=`location:variety`),linetype=3,alpha=1.0,size=0.2) +
      geom_line(data=c,aes(x=W,y=N_cred_nls_up,group=`location:variety`),linetype=3,alpha=1.0,size=0.2)
  }else if(.plot=="conservative"){
    g +
      geom_line(data=c,aes(x=W,y=N_cred_est_lo,group=`location:variety`),linetype=2,alpha=1.0,size=0.2) +
      geom_line(data=c,aes(x=W,y=N_cred_est_up,group=`location:variety`),linetype=2,alpha=1.0,size=0.2)
  }else if(.plot=="estimate+conservative"){
    g +
      geom_line(data=c,aes(x=W,y=N_cred_nls_lo,group=`location:variety`),linetype=3,alpha=1.0,size=0.2) +
      geom_line(data=c,aes(x=W,y=N_cred_nls_up,group=`location:variety`),linetype=3,alpha=1.0,size=0.2) +
      geom_line(data=c,aes(x=W,y=N_cred_est_lo,group=`location:variety`),linetype=2,alpha=1.0,size=0.2) +
      geom_line(data=c,aes(x=W,y=N_cred_est_up,group=`location:variety`),linetype=2,alpha=1.0,size=0.2)
  }else if(.plot=="all"){
    g +
      geom_line(data=r,aes(x=W,y=N_draw,group=.draw),alpha=0.80,color="#ca0020",size=0.20) +
      geom_line(data=c,aes(x=W,y=0,group=`location:variety`),linetype=1,alpha=1.0) +
      geom_line(data=c,aes(x=W,y=N_cred_nls_lo,group=`location:variety`),linetype=3,alpha=1.0,size=0.2) +
      geom_line(data=c,aes(x=W,y=N_cred_nls_up,group=`location:variety`),linetype=3,alpha=1.0,size=0.2) +
      geom_line(data=c,aes(x=W,y=N_cred_est_lo,group=`location:variety`),linetype=2,alpha=1.0,size=0.2) +
      geom_line(data=c,aes(x=W,y=N_cred_est_up,group=`location:variety`),linetype=2,alpha=1.0,size=0.2)
  }else{
    g
  }
  
}

fig4.callout.list <- list(
  .plot=c("all","draws","estimate","conservative","estimate+conservative")
)

fig4.callout.sub <- pmap(fig4.callout.list,~f.fig4.callout(plot.data,
                                                           .location="Minnesota",
                                                           .variety="Russet Burbank",
                                                           .plot=..1))

f.fig4.callout.lab.facet <- function(.location){
  
  ggplot() +
    geom_text(aes(x=0,y=0,label=.location),angle=270,size=2.5) +
    theme_classic() +
    # theme_grey() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_rect(fill="#d2d2d2",color="black"))
  
}
f.fig4.callout.lab.axis.y <- function(){
  
  ggplot() +
    geom_text(aes(x=0,y=0),label=expression("∆ %N"[c]*" [g N 100 g"^-1*"]"),angle=90,size=3,parse=T) +
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank())
  
}
f.fig4.callout.lab.axis.x <- function(){
  
  ggplot() +
    geom_text(aes(x=0,y=0),label=expression("Biomass [Mg ha"^-1*"]"),angle=0,size=3,parse=T) +
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank())
  
}

# fig4.callout.layout <- rbind(c(4,1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,6),
#                              c(4,1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,6),
#                              c(4,1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,6),
#                              c(4,1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,6),
#                              c(4,1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,6),
#                              c(NA,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,NA))
# 
# fig4.callout <- grid.arrange(fig4.callout.sub[[2]],fig4.callout.sub[[3]],fig4.callout.sub[[4]],
#                             f.fig4.callout.lab.axis.y(),f.fig4.callout.lab.axis.x(), f.fig4.callout.lab.facet("Minnesota x Russet Burbank"),
#                              layout_matrix=fig4.callout.layout)

fig4.callout.layout <- rbind(c(NA,NA,3,1,1,1,1,1,NA,2,2,2,2,2,5,NA,NA),
                             c(NA,NA,3,1,1,1,1,1,NA,2,2,2,2,2,5,NA,NA),
                             c(NA,NA,3,1,1,1,1,1,NA,2,2,2,2,2,5,NA,NA),
                             c(NA,NA,3,1,1,1,1,1,NA,2,2,2,2,2,5,NA,NA),
                             c(NA,NA,3,1,1,1,1,1,NA,2,2,2,2,2,5,NA,NA),
                             c(NA,NA,NA,4,4,4,4,4,4,4,4,4,4,4,NA,NA,NA))

fig4.callout <- grid.arrange(fig4.callout.sub[[2]],fig4.callout.sub[[5]],
                            f.fig4.callout.lab.axis.y(),f.fig4.callout.lab.axis.x(), f.fig4.callout.lab.facet("Minnesota x Russet Burbank"),
                             layout_matrix=fig4.callout.layout)

ggsave(filename="manuscript/images/figure6b.pdf",plot=fig4.callout,height=2.2,width=6.0,units="in",scale=1.0,device=cairo_pdf)
ggsave(filename="manuscript/images/figure6b.png",plot=fig4.callout,height=2.2,width=6.0,units="in",scale=1.0,dpi=1000)

# fig4.callout.layout <- rbind(c(2,1,1,1,1,1,1),
#                              c(2,1,1,1,1,1,1),
#                              c(2,1,1,1,1,1,1),
#                              c(2,1,1,1,1,1,1),
#                              c(2,1,1,1,1,1,1),
#                              c(2,1,1,1,1,1,1),
#                              c(NA,3,3,3,3,3,3))
# 
# fig4.callout <- grid.arrange(fig4.callout.sub[[1]],
#                              f.fig4.callout.lab.axis.y(),f.fig4.callout.lab.axis.x(),
#                              layout_matrix=fig4.callout.layout)
# 
# ggsave(filename="manuscript/images/figure4_callout.pdf",plot=fig4.callout,height=3.0,width=3.0,units="in",scale=1.0)
# ggsave(filename="manuscript/images/figure4_callout.png",plot=fig4.callout,height=3.0,width=3.0,units="in",scale=1.0,dpi=1000)


# figure 7 - comparing curves to each other fits -----------------

# .location_comp = "Belgium"
# .variety_comp = "Bintje"
# .location_ref = "Minnesota"
# .variety_ref = "Russet Burbank"

f.fig5 <- function(plot.data,.location_ref,.variety_ref,.location_comp,.variety_comp){
  
  var1 <- "W"
  var2 <- "%N - Diff"
  var3 <- paste(.location_ref,.variety_ref,sep=" - ")
  var4 <- .variety_comp
  var5 <- .location_comp
  
  r <- left_join(
    plot.data$r %>%
      filter(location %in% .location_ref) %>%
      filter(variety %in% .variety_ref) %>%
      rename(N_ref=N) %>%
      rename(location_ref=location,
             variety_ref=variety,
             `location:variety_ref`=`location:variety`),
    plot.data$r %>%
      filter(location %in% .location_comp) %>%
      filter(variety %in% .variety_comp) %>%
      rename(N_comp=N) %>%
      rename(location_comp=location,
             variety_comp=variety,
             `location:variety_comp`=`location:variety`),
    by=c(".draw","W")
  ) %>%
    drop_na()
  
  r <- r %>%
    mutate(N_diff = N_comp - N_ref) %>% 
    group_by(location_ref, variety_ref, `location:variety_ref`, 
             location_comp, variety_comp, `location:variety_comp`, 
             W) %>% 
    # median_qi(N_diff, .width = 0.9) %>% 
    summarize(qs = quantile(`N_diff`,c(0.05,0.50,0.95)), prob = c(0.05,0.50,0.95), .groups="drop") %>%
    pivot_wider(names_from=prob,
                names_prefix="N_diff_",#"N_location:variety_",
                values_from=qs) %>%
    mutate_at(vars(N_diff_0.05,N_diff_0.5,N_diff_0.95),as.numeric) %>%
    mutate(N_class = N_diff_0.05 < 0 & N_diff_0.95 > 0)
  
  r_fix <- r %>%
    slice(1) %>%
    mutate_at(vars(W,N_diff_0.05,N_diff_0.5,N_diff_0.95,N_class),~NA)
  
  r_plot <- r %>%
    bind_rows(r_fix %>% mutate(N_class=TRUE)) %>%
    bind_rows(r_fix %>% mutate(N_class=FALSE))
  
  r_range <- r %>% 
    filter(N_class==TRUE) %>%
    # group_by(location_ref,variety,`location:variety`) %>%
    group_by(location_ref, variety_ref, `location:variety_ref`, 
             location_comp, variety_comp, `location:variety_comp`) %>%
    summarize(range_min=min(W),range_max=max(W),.groups="drop")
  
  if(lengths(r_range[,1])==0){
    r_range <- r %>%
      select(c("location_ref", "variety_ref", "location:variety_ref", 
               "location_comp", "variety_comp", "location:variety_comp")) %>%
      distinct() %>%
      mutate(range_min=NA,
             range_max=NA)
  }

  g <- ggplot() +
    geom_ribbon(data=r_plot,aes(x=W,ymin=N_diff_0.05,ymax=N_diff_0.95),alpha=0.20) + #,fill="#737373"
    geom_point(data=r_plot,aes(x=W,y=N_diff_0.5,group=`location:variety_comp`,color=N_class),alpha=1.0,size=0.2) + #linetype=1,
    geom_line(data=r_plot,aes(x=W,y=0,group=`location:variety_comp`),linetype=1,alpha=1.0) +
    theme_classic() +
    labs(x=var1,
         y=var2,
         title=var3,
         caption=var4) +
    theme(text=element_text(size=8),
          plot.title=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          plot.caption = element_text(hjust=0,size=7), 
          plot.title.position = "plot", 
          plot.caption.position =  "plot") +
    guides(color="none") +
    scale_color_manual(values=c("#ca0020","#0571b0")) +
    scale_x_continuous(limits=c(0,NA),breaks=c(0,5,10,15,20,25,30))
    
    if(var5%in%c("Minnesota","Belgium")) {
      g + 
        scale_y_continuous(limits=c(-0.33,0.33),n.breaks=5) +
        geom_text(data=r_range,aes(x=1,y=0.33,label=paste0("[",format(round(range_min,1),nsmall=1),", ",format(round(range_max,1),nsmall=1),"]")),size=2.0,hjust="inward",vjust=1)
    } else if(var5=="Canada"){
      g + 
        scale_y_continuous(limits=c(-0.5,1.0),n.breaks=5) +
        geom_text(data=r_range,aes(x=1,y=1.0,label=paste0("[",format(round(range_min,1),nsmall=1),", ",format(round(range_max,1),nsmall=1),"]")),size=2.0,hjust="inward",vjust=1)
    } else if(var5=="Argentina"){
      g + 
        scale_y_continuous(limits=c(-0.5,3.0),n.breaks=5) +
        geom_text(data=r_range,aes(x=1,y=3.0,label=paste0("[",format(round(range_min,1),nsmall=1),", ",format(round(range_max,1),nsmall=1),"]")),size=2.0,hjust="inward",vjust=1)
    } else {
      g
    }
  
}

fig5.list <- list(
  location_ref=c("Minnesota"),
  variety_ref=c("Russet Burbank"),
  location_comp=c("Argentina","Argentina","Argentina","Argentina","Argentina","Belgium","Belgium","Canada","Canada","Minnesota","Minnesota","Minnesota","Minnesota","Minnesota"),
  variety_comp=c("Bannock Russet","Gem Russet","Innovator","Markies Russet","Umatilla Russet","Bintje","Charlotte","Russet Burbank","Shepody","Clearwater","Dakota Russet","Easton","Russet Burbank","Umatilla Russet")
)

fig5.sub <- pmap(fig5.list,~f.fig5(plot.data,
                                   .location_ref=..1,
                                   .variety_ref=..2,
                                   .location_comp=..3,
                                   .variety_comp=..4))

f.fig5.lab.facet <- function(.location){
  
  ggplot() +
    geom_text(aes(x=0,y=0,label=.location),angle=270,size=2.5) +
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_rect(fill="#d2d2d2",color="black"))
  
}
f.fig5.lab.axis.y <- function(){
  
  ggplot() +
    geom_text(aes(x=0,y=0),label=expression("∆ %N"[c]*" [g N 100 g"^-1*"]"),angle=90,size=3,parse=T) +
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank())
  
}
f.fig5.lab.axis.x <- function(){
  
  ggplot() +
    geom_text(aes(x=0,y=0),label=expression("Biomass [Mg ha"^-1*"]"),angle=0,size=3,parse=T) +
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank())
  
}
f.fig5.reference <- function(.variety){
  
  ggplot() +
    geom_text(aes(x=0,y=0),label="Reference",angle=0,size=3,parse=T) +
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          text=element_text(size=8),
          plot.title=element_blank(),
          plot.caption = element_text(hjust=0,size=7),
          plot.title.position = "plot",
          plot.caption.position =  "plot") +
    labs(caption=.variety)
  
}

fig5.layout <- rbind(c(19,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,15),
                     c(19,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,15),
                     c(19,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,15),
                     c(19,6,6,6,7,7,7,16,NA,NA,8,8,8,9,9,9,17),
                     c(19,6,6,6,7,7,7,16,NA,NA,8,8,8,9,9,9,17),
                     c(19,6,6,6,7,7,7,16,NA,NA,8,8,8,9,9,9,17),
                     c(19,10,10,10,11,11,11,12,12,12,13,13,13,14,14,14,18),
                     c(19,10,10,10,11,11,11,12,12,12,13,13,13,14,14,14,18),
                     c(19,10,10,10,11,11,11,12,12,12,13,13,13,14,14,14,18),
                     c(NA,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20))

fig5 <- grid.arrange(fig5.sub[[1]],fig5.sub[[2]],fig5.sub[[3]],fig5.sub[[4]],fig5.sub[[5]],
                     fig5.sub[[6]],fig5.sub[[7]],fig5.sub[[8]],fig5.sub[[9]],
                     fig5.sub[[10]],fig5.sub[[11]],fig5.sub[[12]],f.fig5.reference("Russet Burbank"),fig5.sub[[14]],
                     f.fig5.lab.facet("Argentina"),f.fig5.lab.facet("Belgium"),f.fig5.lab.facet("Canada"),f.fig5.lab.facet("Minnesota"),
                     f.fig5.lab.axis.y(),f.fig5.lab.axis.x(),
                     layout_matrix=fig5.layout)

ggsave(filename="manuscript/images/figure7.pdf",plot=fig5,height=4,width=6,units="in",scale=1.0,device=cairo_pdf)
ggsave(filename="manuscript/images/figure7.png",plot=fig5,height=4,width=6,units="in",scale=1.0,dpi=1000)

# figure 7 table - tabular values for comparing curves to each other fits ------------

# .location_comp = "Belgium"
# .variety_comp = "Bintje"
# .location_ref = "Minnesota"
# .variety_ref = "Russet Burbank"

f.fig5.tab <- function(plot.data,.location_ref,.variety_ref,.location_comp,.variety_comp){
  
  var1 <- "W"
  var2 <- "%N - Diff"
  var3 <- paste(.location_ref,.variety_ref,sep=" - ")
  var4 <- .variety_comp
  var5 <- .location_comp
  
  r <- left_join(
    plot.data$r %>%
      filter(location %in% .location_ref) %>%
      filter(variety %in% .variety_ref) %>%
      rename(N_ref=N) %>%
      rename(location_ref=location,
             variety_ref=variety,
             `location:variety_ref`=`location:variety`),
    plot.data$r %>%
      filter(location %in% .location_comp) %>%
      filter(variety %in% .variety_comp) %>%
      rename(N_comp=N) %>%
      rename(location_comp=location,
             variety_comp=variety,
             `location:variety_comp`=`location:variety`),
    by=c(".draw","W")
  ) %>%
    drop_na()
  
  r <- r %>%
    mutate(N_diff = N_comp - N_ref) %>% 
    group_by(location_ref, variety_ref, `location:variety_ref`, 
             location_comp, variety_comp, `location:variety_comp`, 
             W) %>% 
    # median_qi(N_diff, .width = 0.9) %>% 
    summarize(qs = quantile(`N_diff`,c(0.05,0.50,0.95)), prob = c(0.05,0.50,0.95), .groups="drop") %>%
    pivot_wider(names_from=prob,
                names_prefix="N_diff_",#"N_location:variety_",
                values_from=qs) %>%
    mutate_at(vars(N_diff_0.05,N_diff_0.5,N_diff_0.95),as.numeric) %>%
    mutate(N_class = N_diff_0.05 < 0 & N_diff_0.95 > 0)
  
  r_fix <- r %>%
    slice(1) %>%
    mutate_at(vars(W,N_diff_0.05,N_diff_0.5,N_diff_0.95,N_class),~NA)
  
  r_plot <- r %>%
    bind_rows(r_fix %>% mutate(N_class=TRUE)) %>%
    bind_rows(r_fix %>% mutate(N_class=FALSE))
  
  r_range <- r %>% 
    filter(N_class==TRUE) %>%
    # group_by(location_ref,variety,`location:variety`) %>%
    group_by(location_ref, variety_ref, `location:variety_ref`, 
             location_comp, variety_comp, `location:variety_comp`) %>%
    summarize(range_min=min(W),range_max=max(W),.groups="drop")
  
  if(lengths(r_range[,1])==0){
    r_range <- r %>%
      select(c("location_ref", "variety_ref", "location:variety_ref", 
               "location_comp", "variety_comp", "location:variety_comp")) %>%
      distinct() %>%
      mutate(range_min=NA,
             range_max=NA)
  }
  
  out <- r %>%
    filter(N_class==FALSE) %>%
    mutate(across(N_diff_0.5,abs)) %>%
    group_by(location_ref, variety_ref, `location:variety_ref`, 
             location_comp, variety_comp, `location:variety_comp`) %>%
    summarize(across(N_diff_0.5,max),.groups="drop") 
  
  if(lengths(out[,1])==0){
    out <- r %>%
      select(c("location_ref", "variety_ref", "location:variety_ref", 
               "location_comp", "variety_comp", "location:variety_comp")) %>%
      distinct() %>%
      mutate(N_diff_0.5=NA)
  }
  
  out <- out %>%
    bind_cols(r_range %>% select(-c("location_ref", "variety_ref", "location:variety_ref", 
                                    "location_comp", "variety_comp", "location:variety_comp"))) %>%
    rename(N_max_diff=N_diff_0.5,
           W_nodiff_min=range_min,
           W_nodiff_max=range_max)
  
  return(out)
  
}

fig5.tab.list <- list(
  location_ref=c("Minnesota"),
  variety_ref=c("Russet Burbank"),
  location_comp=c("Argentina","Argentina","Argentina","Argentina","Argentina","Belgium","Belgium","Canada","Canada","Minnesota","Minnesota","Minnesota","Minnesota","Minnesota"),
  variety_comp=c("Bannock Russet","Gem Russet","Innovator","Markies Russet","Umatilla Russet","Bintje","Charlotte","Russet Burbank","Shepody","Clearwater","Dakota Russet","Easton","Russet Burbank","Umatilla Russet")
)

fig5.tab <- pmap(fig5.tab.list,~f.fig5.tab(plot.data,
                                           .location_ref=..1,
                                           .variety_ref=..2,
                                           .location_comp=..3,
                                           .variety_comp=..4)) %>% bind_rows()

write_csv(fig5.tab,"manuscript/tables/figure7_table.csv")

# figure 8 - comparing these curves to previous curve fits -----------------

# .location = "Belgium"
# .variety = "Bintje"

f.fig6 <- function(plot.data,.location,.variety){
  
  var1 <- "W"
  var2 <- "%N - Diff"
  var3 <- paste(.location,.variety,sep=" - ")
  var4 <- .variety
  
  c_ref <- plot.data$c %>%
    filter(location %in% .location) %>%
    filter(variety %in% .variety) %>%
    rename(N_ref_0.05=N_0.05,
           N_ref_0.5=N_0.5,
           N_ref_0.95=N_0.95)
  
  parm_ref <- plot.data$parm.fit.sum %>%
    filter(location %in% .location) %>%
    filter(variety %in% .variety)
  
  parm_orig <- plot.data$parm.orig.sum %>%
    filter(location %in% .location) %>%
    filter(variety %in% .variety)
  
  c_orig <- plot.data$c %>%
    filter(location %in% .location) %>%
    filter(variety %in% .variety) %>%
    select(location,variety,`location:variety`,W) %>%
    mutate(alpha1_orig=parm_orig$alpha1_orig,
           alpha2_orig=parm_orig$alpha2_orig) %>%
    mutate(N_orig=alpha1_orig*W^(-alpha2_orig))
  
  c <- left_join(
    c_ref,
    c_orig %>%
      select(W,N_orig),
    by="W"
  ) %>% 
    mutate(N_ref_norm=0,
           N_ref_up=N_ref_0.95-N_ref_0.5,
           N_ref_lo=N_ref_0.05-N_ref_0.5,
           N_orig_norm=N_orig-N_ref_0.5) %>%
    mutate(N_class=case_when(
      (N_orig_norm <= N_ref_up) & (N_orig_norm >= N_ref_lo) ~ T,
      (N_orig_norm > N_ref_up) ~ F,
      (N_orig_norm < N_ref_lo) ~ F,
      TRUE ~ NA)) %>%
    na.omit()
  
  c_range <- c %>% 
    filter(N_class==TRUE) %>%
    group_by(location,variety,`location:variety`) %>%
    summarize(range_min=min(W),range_max=max(W),.groups="drop")
  
  if(lengths(c_range[,1])==0){
    c_range <- c %>%
      select(location,variety,`location:variety`) %>%
      distinct() %>%
      mutate(range_min=NA,
             range_max=NA)
  }
  
  g <- ggplot() +
    geom_ribbon(data=c,aes(x=W,ymin=N_ref_lo,ymax=N_ref_up),alpha=0.20) + #,fill="#737373"
    geom_point(data=c,aes(x=W,y=N_orig_norm,group=`location:variety`,color=N_class),alpha=1.0,size=0.2) + #linetype=1,
    geom_line(data=c,aes(x=W,y=0,group=`location:variety`),linetype=1,alpha=1.0) +
    geom_text(data=c_range,aes(x=1,y=1.5,label=paste0("[",format(round(range_min,1),nsmall=1),", ",format(round(range_max,1),nsmall=1),"]")),size=2.0,hjust="inward",vjust=1) +
    theme_classic() +
    labs(x=var1,
         y=var2,
         title=var3,
         caption=var4) +
    theme(text=element_text(size=8),
          plot.title=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          plot.caption = element_text(hjust=0,size=7), 
          plot.title.position = "plot", 
          plot.caption.position =  "plot") +
    guides(color="none") +
    scale_color_manual(values=c("#ca0020","#0571b0")) +
    scale_y_continuous(limits=c(-1.5,1.5),n.breaks=5) +
    scale_x_continuous(limits=c(0,NA),n.breaks=4)
  
}

fig6.list <- list(
  location=c("Argentina","Argentina","Argentina","Argentina","Argentina","Belgium","Belgium","Canada","Canada"), #,"Minnesota","Minnesota","Minnesota","Minnesota","Minnesota"),
  variety=c("Bannock Russet","Gem Russet","Innovator","Markies Russet","Umatilla Russet","Bintje","Charlotte","Russet Burbank","Shepody") #,"Clearwater","Dakota Russet","Easton","Russet Burbank","Umatilla Russet"),
)

fig6.sub <- pmap(fig6.list,~f.fig6(plot.data,
                                   .location=..1,
                                   .variety=..2))

f.fig6.lab.facet <- function(.location){
  
  ggplot() +
    geom_text(aes(x=0,y=0,label=.location),angle=270,size=2.5) +
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_rect(fill="#d2d2d2",color="black"))
  
}
f.fig6.lab.axis.y <- function(){
  
  ggplot() +
    geom_text(aes(x=0,y=0),label=expression("∆%N"[c]*" [g N 100 g"^-1*"]"),angle=90,size=3,parse=T) +
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank())
  
}
f.fig6.lab.axis.x <- function(){
  
  ggplot() +
    geom_text(aes(x=0,y=0),label=expression("Biomass [Mg ha"^-1*"]"),angle=0,size=3,parse=T) +
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank())
  
}

fig6.layout <- rbind(c(19,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,15),
                     c(19,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,15),
                     c(19,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,15),
                     c(19,6,6,6,7,7,7,16,NA,NA,8,8,8,9,9,9,17),
                     c(19,6,6,6,7,7,7,16,NA,NA,8,8,8,9,9,9,17),
                     c(19,6,6,6,7,7,7,16,NA,NA,8,8,8,9,9,9,17),
                     # c(19,10,10,10,11,11,11,12,12,12,13,13,13,14,14,14,18),
                     # c(19,10,10,10,11,11,11,12,12,12,13,13,13,14,14,14,18),
                     # c(19,10,10,10,11,11,11,12,12,12,13,13,13,14,14,14,18),
                     c(NA,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20))

fig6 <- grid.arrange(fig6.sub[[1]],fig6.sub[[2]],fig6.sub[[3]],fig6.sub[[4]],fig6.sub[[5]],
                     fig6.sub[[6]],fig6.sub[[7]],fig6.sub[[8]],fig6.sub[[9]],
                     # fig6.sub[[10]],fig6.sub[[11]],fig6.sub[[12]],fig6.sub[[13]],fig6.sub[[14]],
                     f.fig6.lab.facet("Argentina"),f.fig6.lab.facet("Belgium"),f.fig6.lab.facet("Canada"),#f.fig6.lab.facet("Minnesota"),
                     f.fig6.lab.axis.y(),f.fig6.lab.axis.x(),
                     layout_matrix=fig6.layout)

ggsave(filename="manuscript/images/figure8.pdf",plot=fig6,height=2.8,width=6,units="in",scale=1.0,device=cairo_pdf)
ggsave(filename="manuscript/images/figure8.png",plot=fig6,height=2.8,width=6,units="in",scale=1.0,dpi=1000)

# figure 8 table - tabular values for the previous curve fit comparison --------------------

# .location = "Canada"
# .variety = "Shepody"

f.fig6.tab <- function(plot.data,.location,.variety){
  
  var1 <- "W"
  var2 <- "%N - Diff"
  var3 <- paste(.location,.variety,sep=" - ")
  var4 <- .variety
  
  c_ref <- plot.data$c %>%
    filter(location %in% .location) %>%
    filter(variety %in% .variety) %>%
    rename(N_ref_0.05=N_0.05,
           N_ref_0.5=N_0.5,
           N_ref_0.95=N_0.95)
  
  parm_ref <- plot.data$parm.fit.sum %>%
    filter(location %in% .location) %>%
    filter(variety %in% .variety)
  
  parm_orig <- plot.data$parm.orig.sum %>%
    filter(location %in% .location) %>%
    filter(variety %in% .variety)
  
  c_orig <- plot.data$c %>%
    filter(location %in% .location) %>%
    filter(variety %in% .variety) %>%
    select(location,variety,`location:variety`,W) %>%
    mutate(alpha1_orig=parm_orig$alpha1_orig,
           alpha2_orig=parm_orig$alpha2_orig) %>%
    mutate(N_orig=alpha1_orig*W^(-alpha2_orig))
  
  c <- left_join(
    c_ref,
    c_orig %>%
      select(W,N_orig),
    by="W"
  ) %>% 
    mutate(N_ref_norm=0,
           N_ref_up=N_ref_0.95-N_ref_0.5,
           N_ref_lo=N_ref_0.05-N_ref_0.5,
           N_orig_norm=N_orig-N_ref_0.5) %>%
    mutate(N_class=case_when(
      (N_orig_norm <= N_ref_up) & (N_orig_norm >= N_ref_lo) ~ T,
      (N_orig_norm > N_ref_up) ~ F,
      (N_orig_norm < N_ref_lo) ~ F,
      TRUE ~ NA)) %>%
    na.omit()
  
  c_range <- c %>% 
    filter(N_class==TRUE) %>%
    group_by(location,variety,`location:variety`) %>%
    summarize(range_min=min(W),range_max=max(W),.groups="drop")
  
  if(lengths(c_range[,1])==0){
    c_range <- c %>%
      select(location,variety,`location:variety`) %>%
      distinct() %>%
      mutate(range_min=NA,
             range_max=NA)
  }
  
  out <- c %>%
    filter(N_class==FALSE) %>%
    mutate(across(N_orig_norm,abs)) %>%
    group_by(location,variety,location:variety) %>%
    summarize(across(N_orig_norm,max),.groups="drop") %>%
    bind_cols(c_range %>% select(-c("location","variety","location:variety"))) %>%
    rename(N_max_diff=N_orig_norm,
           W_nodiff_min=range_min,
           W_nodiff_max=range_max)
  
  return(out)
  
}

fig6.tab.list <- list(
  location=c("Argentina","Argentina","Argentina","Argentina","Argentina","Belgium","Belgium","Canada","Canada"), #,"Minnesota","Minnesota","Minnesota","Minnesota","Minnesota"),
  variety=c("Bannock Russet","Gem Russet","Innovator","Markies Russet","Umatilla Russet","Bintje","Charlotte","Russet Burbank","Shepody") #,"Clearwater","Dakota Russet","Easton","Russet Burbank","Umatilla Russet"),
)

fig6.tab <- pmap(fig6.tab.list,~f.fig6.tab(plot.data,
                                           .location=..1,
                                           .variety=..2)) %>% bind_rows()

write_csv(fig6.tab,"manuscript/tables/figure8_table.csv")


# figure S1 - plateau model fit with point data for each date shown for each variety x location ------------------

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
    mutate_at(vars(location,variety,`location:variety`), ~fct_drop(.)) %>%
    left_join(plot.data$tab1 %>%
                select(location,variety,`location:variety`,alpha1_0.5,alpha2_0.5),
              by = c("location", "variety", "location:variety")) %>%
    mutate(Nc=alpha1_0.5*(W^(-alpha2_0.5))) %>%
    mutate(Nclass=case_when(
      N>Nc ~ "Surplus",
      N<Nc ~ "Deficit",
      N==Nc ~ "Optimal")
    ) %>%
    filter(is.na(Nclass)!=T)
  
  p <- plot.data$p %>%
    filter(location%in%.location) %>%
    filter(variety%in%.variety) %>%
    mutate_at(vars(location,variety,`location:variety`), ~fct_drop(.))
  
  d.sum.index <- d %>%
    group_by(index,Nclass) %>%
    count(Nclass) %>%
    pivot_wider(names_from="Nclass",values_from=n)
  
  d.sum.group <- d %>%
    group_by(group,Nclass) %>%
    count(Nclass) %>%
    pivot_wider(names_from="Nclass",values_from=n)
  
  d.colors <- c("#ca0020","#0571b0")
  
  ggplot() +
    geom_line(data=c,aes(x=W,y=N_0.5),linetype=1,alpha=1.0) +
    geom_line(data=p,aes(x=W_0.5,y=N_0.5,group=index),linetype=1,alpha=0.5) +
    geom_point(data=d,aes(x=W,y=N,color=Nclass),alpha=0.33) +
    geom_text(data=d.sum.index,aes(x=0,y=0,label=Deficit),color=d.colors[[1]],size=3,hjust=0,vjust=0) + # paste0("Deficit: ",Deficit)
    geom_text(data=d.sum.index,aes(x=0,y=6,label=Surplus),color=d.colors[[2]],size=3,hjust=0,vjust=1) + 
    facet_wrap(vars(as.numeric(index)),ncol=8) +
    labs(x = expression("Biomass [Mg ha"^-1*"]"),
         y = expression("%N [g N 100 g"^-1*"]"),
         title = paste(.location,.variety,sep=" - "),
         color = "N Status\nClassification",
         subtitle = paste0("Deficit: ",d.sum.group$Deficit,", ","Surplus: ",d.sum.group$Surplus)) + 
    guides(color=guide_legend()) +
    coord_cartesian(xlim=c(0,NA),ylim=c(0,6.0)) +
    theme_classic() +
    theme(legend.position = "right") +
    scale_color_manual(values=d.colors)
  
}

f.appx1.join <- function(){
  
  f.appx1_title <- function(){
    
    appx1_title <- "Figure S1: Linear-plateau curves fitted for each experimental sampling date"
    appx1_caption <- "Fitted critical N dilution curve for each variety within location shown for each level of index with the fitted linear-plateau 
curve and experimental data observed for that date. Experimental data points are classified depending on whether the
N concentration [%N] for that given level of biomass is less than the critical N concentration [%Nc] (i.e., Deficit) or is 
greater than %Nc (i.e., Surplus). The total number of experimental observations classified as Deficit or Surplus is 
summarized for each level of index and also for the each variety within location."
    
    ggplot() +
      labs(title = appx1_title,
           subtitle = appx1_caption) +
      theme_classic() +
      theme(axis.line = element_blank())
    
  }
  appx1_title <- f.appx1_title()
  g.appx1_title <- ggplotGrob(appx1_title)
  fg.appx1_title <- gtable_frame(g.appx1_title, height = unit(5*(0.33/7), "null"), width = unit(6, "null"))
  
  appx1_a <- f.appx1(plot.data,.location=c("Argentina"),.variety=c("Bannock Russet"))
  g.appx1_a <- ggplotGrob(appx1_a)
  fg.appx1_a <- gtable_frame(g.appx1_a, height = unit(5*(2/7), "null"), width = unit(6, "null"))
  ggsave(filename="manuscript/images/figureS1_a.png",plot=fg.appx1_a,height=5*(2/7)*1.40,width=6,scale=1.5,limitsize=F,dpi=1000)
  ggsave(filename="manuscript/images/figureS1_a.pdf",plot=fg.appx1_a,height=5*(2/7)*1.40,width=6,scale=1.5,limitsize=F,device=cairo_pdf)

  appx1_b <- f.appx1(plot.data,.location=c("Argentina"),.variety=c("Gem Russet"))
  g.appx1_b <- ggplotGrob(appx1_b)
  fg.appx1_b <- gtable_frame(g.appx1_b, height = unit(5*(3/7), "null"), width = unit(6, "null"))
  ggsave(filename="manuscript/images/figureS1_b.png",plot=fg.appx1_b,height=5*(3/7)*1.25,width=6,scale=1.5,limitsize=F,dpi=1000)
  ggsave(filename="manuscript/images/figureS1_b.pdf",plot=fg.appx1_b,height=5*(3/7)*1.25,width=6,scale=1.5,limitsize=F,device=cairo_pdf)
  
  appx1_c <- f.appx1(plot.data,.location=c("Argentina"),.variety=c("Innovator"))
  g.appx1_c <- ggplotGrob(appx1_c)
  fg.appx1_c <- gtable_frame(g.appx1_c, height = unit(5*(3/7), "null"), width = unit(6, "null"))
  ggsave(filename="manuscript/images/figureS1_c.png",plot=fg.appx1_c,height=5*(3/7)*1.25,width=6,scale=1.5,limitsize=F,dpi=1000)
  ggsave(filename="manuscript/images/figureS1_c.pdf",plot=fg.appx1_c,height=5*(3/7)*1.25,width=6,scale=1.5,limitsize=F,device=cairo_pdf)
  
  appx1_d <- f.appx1(plot.data,.location=c("Argentina"),.variety=c("Markies Russet"))
  g.appx1_d <- ggplotGrob(appx1_d)
  fg.appx1_d <- gtable_frame(g.appx1_d, height = unit(5*(2/7), "null"), width = unit(6, "null"))
  ggsave(filename="manuscript/images/figureS1_d.png",plot=fg.appx1_d,height=5*(2/7)*1.40,width=6,scale=1.5,limitsize=F,dpi=1000)
  ggsave(filename="manuscript/images/figureS1_d.pdf",plot=fg.appx1_d,height=5*(2/7)*1.40,width=6,scale=1.5,limitsize=F,device=cairo_pdf)
  
  appx1_e <- f.appx1(plot.data,.location=c("Argentina"),.variety=c("Umatilla Russet"))
  g.appx1_e <- ggplotGrob(appx1_e)
  fg.appx1_e <- gtable_frame(g.appx1_e, height = unit(5*(2/7), "null"), width = unit(6, "null"))
  ggsave(filename="manuscript/images/figureS1_e.png",plot=fg.appx1_e,height=5*(2/7)*1.40,width=6,scale=1.5,limitsize=F,dpi=1000)
  ggsave(filename="manuscript/images/figureS1_e.pdf",plot=fg.appx1_e,height=5*(2/7)*1.40,width=6,scale=1.5,limitsize=F,device=cairo_pdf)
  
  appx1_f <- f.appx1(plot.data,.location=c("Belgium"),.variety=c("Bintje"))
  g.appx1_f <- ggplotGrob(appx1_f)
  fg.appx1_f <- gtable_frame(g.appx1_f, height = unit(5*(7/7), "null"), width = unit(6, "null"))
  ggsave(filename="manuscript/images/figureS1_f.png",plot=fg.appx1_f,height=5*(7/7)*1.05,width=6,scale=1.5,limitsize=F,dpi=1000)
  ggsave(filename="manuscript/images/figureS1_f.pdf",plot=fg.appx1_f,height=5*(7/7)*1.05,width=6,scale=1.5,limitsize=F,device=cairo_pdf)
  
  appx1_g <- f.appx1(plot.data,.location=c("Belgium"),.variety=c("Charlotte"))
  g.appx1_g <- ggplotGrob(appx1_g)
  fg.appx1_g <- gtable_frame(g.appx1_g, height = unit(5*(3/7), "null"), width = unit(6, "null"))
  ggsave(filename="manuscript/images/figureS1_g.png",plot=fg.appx1_g,height=5*(3/7)*1.25,width=6,scale=1.5,limitsize=F,dpi=1000)
  ggsave(filename="manuscript/images/figureS1_g.pdf",plot=fg.appx1_g,height=5*(3/7)*1.25,width=6,scale=1.5,limitsize=F,device=cairo_pdf)
  
  appx1_h <- f.appx1(plot.data,.location=c("Canada"),.variety=c("Russet Burbank"))
  g.appx1_h <- ggplotGrob(appx1_h)
  fg.appx1_h <- gtable_frame(g.appx1_h, height = unit(5*(4/7), "null"), width = unit(6, "null"))
  ggsave(filename="manuscript/images/figureS1_h.png",plot=fg.appx1_h,height=5*(4/7)*1.20,width=6,scale=1.5,limitsize=F,dpi=1000)
  ggsave(filename="manuscript/images/figureS1_h.pdf",plot=fg.appx1_h,height=5*(4/7)*1.20,width=6,scale=1.5,limitsize=F,device=cairo_pdf)
  
  appx1_i <- f.appx1(plot.data,.location=c("Canada"),.variety=c("Shepody"))
  g.appx1_i <- ggplotGrob(appx1_i)
  fg.appx1_i <- gtable_frame(g.appx1_i, height = unit(5*(4/7), "null"), width = unit(6, "null"))
  ggsave(filename="manuscript/images/figureS1_i.png",plot=fg.appx1_i,height=5*(4/7)*1.20,width=6,scale=1.5,limitsize=F,dpi=1000)
  ggsave(filename="manuscript/images/figureS1_i.pdf",plot=fg.appx1_i,height=5*(4/7)*1.20,width=6,scale=1.5,limitsize=F,device=cairo_pdf)
  
  appx1_j <- f.appx1(plot.data,.location=c("Minnesota"),.variety=c("Clearwater"))
  g.appx1_j <- ggplotGrob(appx1_j)
  fg.appx1_j <- gtable_frame(g.appx1_j, height = unit(5*(2/7), "null"), width = unit(6, "null"))
  ggsave(filename="manuscript/images/figureS1_j.png",plot=fg.appx1_j,height=5*(2/7)*1.40,width=6,scale=1.5,limitsize=F,dpi=1000)
  ggsave(filename="manuscript/images/figureS1_j.pdf",plot=fg.appx1_j,height=5*(2/7)*1.40,width=6,scale=1.5,limitsize=F,device=cairo_pdf)
  
  appx1_k <- f.appx1(plot.data,.location=c("Minnesota"),.variety=c("Dakota Russet"))
  g.appx1_k <- ggplotGrob(appx1_k)
  fg.appx1_k <- gtable_frame(g.appx1_k, height = unit(5*(2/7), "null"), width = unit(6, "null"))
  ggsave(filename="manuscript/images/figureS1_k.png",plot=fg.appx1_k,height=5*(2/7)*1.40,width=6,scale=1.5,limitsize=F,dpi=1000)
  ggsave(filename="manuscript/images/figureS1_k.pdf",plot=fg.appx1_k,height=5*(2/7)*1.40,width=6,scale=1.5,limitsize=F,device=cairo_pdf)
  
  appx1_l <- f.appx1(plot.data,.location=c("Minnesota"),.variety=c("Easton"))
  g.appx1_l <- ggplotGrob(appx1_l)
  fg.appx1_l <- gtable_frame(g.appx1_l, height = unit(5*(2/7), "null"), width = unit(6, "null"))
  ggsave(filename="manuscript/images/figureS1_l.png",plot=fg.appx1_l,height=5*(2/7)*1.40,width=6,scale=1.5,limitsize=F,dpi=1000)
  ggsave(filename="manuscript/images/figureS1_l.pdf",plot=fg.appx1_l,height=5*(2/7)*1.40,width=6,scale=1.5,limitsize=F,device=cairo_pdf)
  
  appx1_m <- f.appx1(plot.data,.location=c("Minnesota"),.variety=c("Russet Burbank"))
  g.appx1_m <- ggplotGrob(appx1_m)
  fg.appx1_m <- gtable_frame(g.appx1_m, height = unit(5*(7/7), "null"), width = unit(6, "null"))
  ggsave(filename="manuscript/images/figureS1_m.png",plot=fg.appx1_m,height=5*(7/7)*1.05,width=6,scale=1.5,limitsize=F,dpi=1000)
  ggsave(filename="manuscript/images/figureS1_m.pdf",plot=fg.appx1_m,height=5*(7/7)*1.05,width=6,scale=1.5,limitsize=F,device=cairo_pdf)
  
  appx1_n <- f.appx1(plot.data,.location=c("Minnesota"),.variety=c("Umatilla Russet"))
  g.appx1_n <- ggplotGrob(appx1_n)
  fg.appx1_n <- gtable_frame(g.appx1_n, height = unit(5*(2/7), "null"), width = unit(6, "null"))
  ggsave(filename="manuscript/images/figureS1_n.png",plot=fg.appx1_n,height=5*(2/7)*1.40,width=6,scale=1.5,limitsize=F,dpi=1000)
  ggsave(filename="manuscript/images/figureS1_n.pdf",plot=fg.appx1_n,height=5*(2/7)*1.40,width=6,scale=1.5,limitsize=F,device=cairo_pdf)
  
  fg <- rbind(fg.appx1_title,fg.appx1_a,fg.appx1_b,fg.appx1_c,fg.appx1_d,fg.appx1_e,fg.appx1_f,fg.appx1_g,fg.appx1_h,fg.appx1_i,fg.appx1_j,fg.appx1_k,fg.appx1_l,fg.appx1_m,fg.appx1_n,size = "first")
  fg$widths <- unit.pmax(fg.appx1_title$widths, fg.appx1_a$widths, fg.appx1_b$widths, fg.appx1_c$widths, fg.appx1_d$widths, fg.appx1_e$widths, fg.appx1_f$widths, fg.appx1_g$widths, fg.appx1_h$widths, fg.appx1_i$widths, fg.appx1_j$widths, fg.appx1_k$widths, fg.appx1_l$widths, fg.appx1_m$widths, fg.appx1_n$widths)
  
  return(fg)
}

appx1 <- f.appx1.join()

# ggsave(filename="manuscript/images/figureS1.pdf",plot=appx1,height=40,width=6,scale=1.5,limitsize=F,device=cairo_pdf)
# ggsave(filename="manuscript/images/figureS1.png",plot=appx1,height=40,width=6,scale=1.5,limitsize=F,dpi=1000)

# figure S2 - full comparing curves to each other fits -----------------

# .location = "Minnesota"
# .variety = "Russet Burbank"

f.appx2 <- function(plot.data,.location,.variety){
  
  var1 <- expression("Biomass [Mg ha"^-1*"]")
  var2 <- expression("∆%N"[c]*" [g N 100 g"^-1*"]")
  var3 <- paste(.location,.variety,sep=" - ")
  var4 <- paste(.location,str_replace(.variety," ","."),sep="_")
  
  r <- left_join(
    plot.data$r %>%
      filter(location %in% .location) %>%
      filter(variety %in% .variety) %>%
      # filter(location %in% .location_ref) %>%
      # filter(variety %in% .variety_ref) %>%
      rename(N_ref=N) %>%
      rename(location_ref=location,
             variety_ref=variety,
             `location:variety_ref`=`location:variety`),
    plot.data$r %>%
      # filter(location %in% .location_comp) %>%
      # filter(variety %in% .variety_comp) %>%
      rename(N_comp=N) %>%
      rename(location_comp=location,
             variety_comp=variety,
             `location:variety_comp`=`location:variety`),
    by=c(".draw","W")
  ) %>%
    filter(`location:variety_ref`!=`location:variety_comp`) %>%
    drop_na()
  
  r <- r %>%
    mutate(N_diff = N_comp - N_ref) %>% 
    group_by(location_ref, variety_ref, `location:variety_ref`, 
             location_comp, variety_comp, `location:variety_comp`, 
             W) %>% 
    summarize(qs = quantile(`N_diff`,c(0.05,0.50,0.95)), prob = c(0.05,0.50,0.95), .groups="drop") %>%
    pivot_wider(names_from=prob,
                names_prefix="N_diff_",
                values_from=qs) %>%
    mutate_at(vars(N_diff_0.05,N_diff_0.5,N_diff_0.95),as.numeric) %>%
    mutate(N_class = N_diff_0.05 < 0 & N_diff_0.95 > 0)
  
  r_range <- r %>%
    filter(N_class==TRUE) %>%
    group_by(location_comp,variety_comp,`location:variety_comp`) %>%
    summarize(range_min=min(W),range_max=max(W),.groups="drop")
  
  r_fix <- r %>%
    slice(1) %>%
    mutate_at(vars(W,N_diff_0.05,N_diff_0.5,N_diff_0.95,N_class),~NA)
  
  r_plot <- r %>%
    mutate(N_class = !N_class) %>%
    bind_rows(r_fix %>% mutate(N_class=TRUE)) %>%
    bind_rows(r_fix %>% mutate(N_class=FALSE))
  
  g <- ggplot() +
    geom_ribbon(data=r_plot,aes(x=W,ymin=N_diff_0.05,ymax=N_diff_0.95),alpha=0.20) + #,fill="#737373"
    geom_point(data=r_plot,aes(x=W,y=N_diff_0.5,group=`location:variety_comp`,color=N_class),alpha=1.0,size=0.2) + #linetype=1,
    geom_line(data=r_plot,aes(x=W,y=0,group=`location:variety_comp`),linetype=1,alpha=1.0) +
    geom_text(data=r_range,aes(x=1,y=3,label=paste0("[",format(round(range_min,1),nsmall=1),", ",format(round(range_max,1),nsmall=1),"]")),size=2.5,hjust="inward",vjust=1) +
    theme_classic() +
    # facet_wrap(vars(`location:variety_comp`),scales="free_x",ncol=7) + 
    facet_wrap(vars(location_comp,variety_comp),scales="free_x",ncol=7) + 
    labs(x=var1,
         y=var2,
         title=var3,
         color="Significant\nDifference") +
    # guides(color="none") +
    scale_color_manual(values=c("#0571b0","#ca0020")) +
    scale_y_continuous(n.breaks=4) +
    coord_cartesian(ylim=c(-3.0,3.0)) +
    scale_x_continuous(limits=c(0,NA),n.breaks=4) #,minor_breaks=c(5,15,25) ,breaks=c(0,10,20,30)
  
}

appx2.list <- list(
  location=c("Argentina","Argentina","Argentina","Argentina","Argentina","Belgium","Belgium","Canada","Canada","Minnesota","Minnesota","Minnesota","Minnesota","Minnesota"),
  variety=c("Bannock Russet","Gem Russet","Innovator","Markies Russet","Umatilla Russet","Bintje","Charlotte","Russet Burbank","Shepody","Clearwater","Dakota Russet","Easton","Russet Burbank","Umatilla Russet")
)

appx2_sub <- pmap(appx2.list,~f.appx2(plot.data,
                                      .location=..1,
                                      .variety=..2))

f.appx2_title <- function(){
  
  appx2_title <- "         Figure S2: Pairwise comparisons of differences in critical N concertation for variety within location"
  appx2_caption <- "           Pairwise comparison of the difference in critical N concentration values [∆%Nc] between the critical N concentration [%Nc]
           for a given reference curve (i.e., title) and %Nc for all levels of variety within location evaluated in the present study
           (i.e., facets). The grey shaded region represents the 90% credible region (lower bound, 5% quantile; upper bound, 
           95% quantile) for ∆%Nc. The colored points represent the median value for ∆%Nc at a given biomass level where blue or
           red color respectively indicate that credible region for ∆%Nc does or does not contain zero. The solid black line represents
           a constant value of zero (i.e., %Nc for reference curve). The range of biomass values for which ∆%Nc is not significantly
           different (i.e., credible region contains zero) is given in brackets."
  
  ggplot() +
    labs(title = appx2_title,
         subtitle = appx2_caption) +
    theme_classic() +
    theme(axis.line = element_blank())
  
}
appx2_title <- f.appx2_title()
# g.appx1_title <- ggplotGrob(appx1_title)
# fg.appx1_title <- gtable_frame(g.appx1_title, height = unit(5*(0.33/7), "null"), width = unit(6, "null"))


# appx2.layout <- rbind(1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,9,10,10,10,11,11,11,12,12,12,13,13,13,14,14,14,15,15,15)
appx2.layout <- rbind(1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13,14,14,15,15)

appx2 <- grid.arrange(appx2_title,
                      appx2_sub[[1]],appx2_sub[[2]],appx2_sub[[3]],appx2_sub[[4]],appx2_sub[[5]],
                      appx2_sub[[6]],appx2_sub[[7]],appx2_sub[[8]],appx2_sub[[9]],appx2_sub[[10]],
                      appx2_sub[[11]],appx2_sub[[12]],appx2_sub[[13]],appx2_sub[[14]],
                      layout_matrix=appx2.layout)

# ggsave(filename="manuscript/images/figureS2.pdf",plot=appx2,scale=1.5,height=35,width=6,limitsize=F,device=cairo_pdf)
# ggsave(filename="manuscript/images/figureS2.png",plot=appx2,scale=1.5,height=35,width=6,limitsize=F,dpi=1000)

ggsave(filename="manuscript/images/figureS2_a.pdf",plot=appx2_sub[[1]],scale=1.5,height=2.5,width=6,limitsize=F,device=cairo_pdf)
ggsave(filename="manuscript/images/figureS2_a.png",plot=appx2_sub[[1]],scale=1.5,height=2.5,width=6,limitsize=F,dpi=1000)

ggsave(filename="manuscript/images/figureS2_b.pdf",plot=appx2_sub[[2]],scale=1.5,height=2.5,width=6,limitsize=F,device=cairo_pdf)
ggsave(filename="manuscript/images/figureS2_b.png",plot=appx2_sub[[2]],scale=1.5,height=2.5,width=6,limitsize=F,dpi=1000)

ggsave(filename="manuscript/images/figureS2_c.pdf",plot=appx2_sub[[3]],scale=1.5,height=2.5,width=6,limitsize=F,device=cairo_pdf)
ggsave(filename="manuscript/images/figureS2_c.png",plot=appx2_sub[[3]],scale=1.5,height=2.5,width=6,limitsize=F,dpi=1000)

ggsave(filename="manuscript/images/figureS2_d.pdf",plot=appx2_sub[[4]],scale=1.5,height=2.5,width=6,limitsize=F,device=cairo_pdf)
ggsave(filename="manuscript/images/figureS2_d.png",plot=appx2_sub[[4]],scale=1.5,height=2.5,width=6,limitsize=F,dpi=1000)

ggsave(filename="manuscript/images/figureS2_e.pdf",plot=appx2_sub[[5]],scale=1.5,height=2.5,width=6,limitsize=F,device=cairo_pdf)
ggsave(filename="manuscript/images/figureS2_e.png",plot=appx2_sub[[5]],scale=1.5,height=2.5,width=6,limitsize=F,dpi=1000)

ggsave(filename="manuscript/images/figureS2_f.pdf",plot=appx2_sub[[6]],scale=1.5,height=2.5,width=6,limitsize=F,device=cairo_pdf)
ggsave(filename="manuscript/images/figureS2_f.png",plot=appx2_sub[[6]],scale=1.5,height=2.5,width=6,limitsize=F,dpi=1000)

ggsave(filename="manuscript/images/figureS2_g.pdf",plot=appx2_sub[[7]],scale=1.5,height=2.5,width=6,limitsize=F,device=cairo_pdf)
ggsave(filename="manuscript/images/figureS2_g.png",plot=appx2_sub[[7]],scale=1.5,height=2.5,width=6,limitsize=F,dpi=1000)

ggsave(filename="manuscript/images/figureS2_h.pdf",plot=appx2_sub[[8]],scale=1.5,height=2.5,width=6,limitsize=F,device=cairo_pdf)
ggsave(filename="manuscript/images/figureS2_h.png",plot=appx2_sub[[8]],scale=1.5,height=2.5,width=6,limitsize=F,dpi=1000)

ggsave(filename="manuscript/images/figureS2_i.pdf",plot=appx2_sub[[9]],scale=1.5,height=2.5,width=6,limitsize=F,device=cairo_pdf)
ggsave(filename="manuscript/images/figureS2_i.png",plot=appx2_sub[[9]],scale=1.5,height=2.5,width=6,limitsize=F,dpi=1000)

ggsave(filename="manuscript/images/figureS2_j.pdf",plot=appx2_sub[[10]],scale=1.5,height=2.5,width=6,limitsize=F,device=cairo_pdf)
ggsave(filename="manuscript/images/figureS2_j.png",plot=appx2_sub[[10]],scale=1.5,height=2.5,width=6,limitsize=F,dpi=1000)

ggsave(filename="manuscript/images/figureS2_k.pdf",plot=appx2_sub[[11]],scale=1.5,height=2.5,width=6,limitsize=F,device=cairo_pdf)
ggsave(filename="manuscript/images/figureS2_k.png",plot=appx2_sub[[11]],scale=1.5,height=2.5,width=6,limitsize=F,dpi=1000)

ggsave(filename="manuscript/images/figureS2_l.pdf",plot=appx2_sub[[12]],scale=1.5,height=2.5,width=6,limitsize=F,device=cairo_pdf)
ggsave(filename="manuscript/images/figureS2_l.png",plot=appx2_sub[[12]],scale=1.5,height=2.5,width=6,limitsize=F,dpi=1000)

ggsave(filename="manuscript/images/figureS2_m.pdf",plot=appx2_sub[[13]],scale=1.5,height=2.5,width=6,limitsize=F,device=cairo_pdf)
ggsave(filename="manuscript/images/figureS2_m.png",plot=appx2_sub[[13]],scale=1.5,height=2.5,width=6,limitsize=F,dpi=1000)

ggsave(filename="manuscript/images/figureS2_n.pdf",plot=appx2_sub[[14]],scale=1.5,height=2.5,width=6,limitsize=F,device=cairo_pdf)
ggsave(filename="manuscript/images/figureS2_n.png",plot=appx2_sub[[14]],scale=1.5,height=2.5,width=6,limitsize=F,dpi=1000)


# END --------------------
