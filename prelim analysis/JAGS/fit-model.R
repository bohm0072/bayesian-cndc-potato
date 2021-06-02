# Fit hierarchical model for preliminary analysis of Bayesian methods for CNDC
# Using framework from Makowski et al. (2020) - doi:10.1016/j.eja.2020.126076

##### Initialization #####
library(tidyverse)
library(lubridate)
library(stringr)

# library(rjags)
# renv::remove("rjags")

##### Import Data #####

# This csv file contains both unpublished data from Rosen Lab and previously published data from Giletto et. al (2020), Appendix A and B - doi:10.1016/j.eja.2020.126114
# The data is formatted to both 1) work with the framework from Makowski et al. (2020) - doi:10.1016/j.eja.2020.126076, and 2) provide necessary metadata to allow for parameterizing multiple models in the same workflow

data <- read_csv("data.csv",col_types="cccdcdd") 

# Metadata

# study: name of study author
  # Bohman: unpublished data from Rosen Lab
  # Giletto: data Giletto et. al (2020), Appendix A - doi:10.1016/j.eja.2020.126114

# location: this is nested w/in study
# note: we may want to fit and compare curves by location
  # Bohman:
    # Minnesota
  # Giletto:
    # Canada
    # Argentina

# variety: this is nested w/in location and study
# note: we may want to fit and compare curves by variety
# note: some varieties are in common (e.g., Russet Burbank)
# note: some varieties are the same but with different names (e.g., Umatilla and Umatilla Russet)
  # Bohman, Minnesota:
    # Russet Burbank
    # Clearwater
    # Umatilla
    # Dakota Russet
    # Easton
  # Giletto, Canada:
    # Russet Burbank
    # Shepody
  # Giletto, Argentina:
    # Innovator
    # Gem Russet
    # Umatilla Russet
    # Bannock Russet
    # Markies Russet

# rate_n_kgha: these are the treatment levels used for a particular study, location, and experiment year 
# note: these levels are continuous and won't be explicitly used to compare models
# note: these levels are implicitly considered in the hierarchical model (e.g. as "replicates") 

# Date: these are the dates on which samples were collected and are explictly used in the hierarchical model fit
# note: the value of Date for the Giletto study are not true dates (based on DAP/DAE and origin of Jan 1 and year of study)

# W: this is the value of dry wt. wholeplant biomass in units of Mg d.w. ha-1

# N: this is the value of whole plant N concentration in units of g N 100 g-1 d.w.

##### Define Models #####

# These are the models as defined in Makowski et al. (2020), Appendix D

f.makowski_model_1 <- function(){

  modelstring="
model{

	for (i in 1:Q)
	{
		W[i]~dnorm(mu[i],tau_b)	
		N[i]~dnorm(Nc[Date[i]],tau_n)
		mu[i]<-min(Bmax[Date[i]],Bmax[Date[i]]+S[Date[i]]*(N[i]-Nc[Date[i]]))
	}

	for (j in 1:K)
	{
		Nc[j]=A1*Bmax[j]^(-A2)
		Bmax[j]~dnorm(Mu_Bmax,Prec_Bmax)T(0,)
		S[j]~dnorm(Mu_S,Prec_S)T(0,)
			}

			#Weakly informative
			Mu_Bmax~dnorm(10,10)
			Mu_S~dnorm(5,2)
			A1~dunif(3,12)
			A2~dunif(0.2,0.8)

			Prec_Bmax~dgamma(0.001,0.001)
			Prec_S~dgamma(0.001,0.001)
			tau_b~dgamma(0.001,0.001)
			tau_n~dgamma(0.001,0.001)
			
}
"

writeLines(modelstring,con="model.txt")

}
f.makowski_model_2 <- function(){
  
  modelstring="
model{

	for (i in 1:Q)
	{
		W[i]~dnorm(mu[i],tau_b)	
		N[i]~dnorm(Nc[Date[i]]+Theta[Date[i]],tau_n)
		mu[i]<-min(Bmax[Date[i]],Bmax[Date[i]]+S[Date[i]]*(N[i]-Nc[Date[i]]))
	}

	for (j in 1:K)
	{
		Nc[j]=A1*Bmax[j]^(-A2)
		Bmax[j]~dnorm(Mu_Bmax,Prec_Bmax)T(0,)
		S[j]~dnorm(Mu_S,Prec_S)T(0,)
		Theta[j]~dnorm(0,tau_t)
			}

			#Weakly informative
			Mu_Bmax~dnorm(6,0.1)
			Mu_S~dnorm(0,0.1)
			A1~dunif(2,6)
			A2~dunif(0,0.5)
			
			Prec_Bmax~dgamma(0.001,0.001)
			Prec_S~dgamma(0.001,0.001)
			tau_b~dgamma(0.001,0.001)
			tau_n~dgamma(0.001,0.001)
			tau_t~dgamma(0.001,0.001)
}
"
writeLines(modelstring,con="model.txt")

}

##### Define Fit Function #####

# An implementation of the method used by Makowski et al. (2020)
# Initially designed to work with purrr::map
# Flexible to accept new JAGS models (e.g., using dataset ID (study X location X variety) as an additional hierarchical structure)

f.fit <- function(model,data,study,location,variety){
  
  if (model=="makowski_model_1"){
    f.makowski_model_1()
  } else if (model=="makowski_model_2"){
    f.makowski_model_2()
  } else {
    writeLines("",con="model.txt")
  }
  
  d <- data %>%
    filter(study%in%study) %>% #=="Bohman") %>% #
    filter(location%in%location) %>% #=="Minnesota") %>% #
    filter(variety%in%variety) #=="Russet Burbank") #
  
  d <- d %>%
    mutate(Date=factor(Date)) %>%
    mutate(Date=as.numeric(Date))
  
  Q <- length(d$Date)
  K <- length(unique(d$Date))
  
  model <- jags.model('model.txt',data=list('W'=d$W,'N'=d$N,'Date'=d$Date,'Q'=Q,'K'=K),
                      n.chains=3, n.adapt=10000)
  
  output <- jags.samples(model,variable.names=c("A1","A2"),n.iter=1000L,thin=1L)
  
  parms <- tibble(A1=as.vector(output$A1),
                A2=as.vector(output$A2))
  
  out <- list(model=model,
              output=output,
              parms=parms)
  
  return(out)
    
}

##### Example Fit and Visualization #####

# Example fit using makowski_model_1 for Bohman data and only Russet Burbank
# Using priors as defined in the original makowski_model_1 (appropriate?...)

fit <- f.fit(model = "makowski_model_1",
             data = data,
             study = "Bohman",
             location = "Minnesota",
             variety = "Russet Burbank"
             )

# Based on conventional nls fit methods, expected value of A1=5.00 and A2=0.45
# Based on bootstrapped conventional nls fit methods, expected range of A1=[4.40,6.14] and A2=[0.36,0.56]

ggplot(fit$parms,aes(A1,A2)) +
  geom_point() +
  theme_classic()

##### END #####