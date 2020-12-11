# Analyzing uncertainty in critical nitrogen dilution curves
# Makowski et al. (2020) - doi:10.1016/j.eja.2020.126076

# Appendix D

##### Header #####

###Data

# Q=total number of biomass observations

# K=number of dates

# W=column of data including biomass observations

# N=column of data including observations of nitrogen concentration

# Date=column with the indices identifying the different dates of the dataset



###Model parameters

# Nc=Critical nitrogen concentration

# Bmax=maximum biomass value in a specific date

# S=slope of the linear-plus-plateau function

# W=biomass increase per unit of nitrogen concentration

# A1 and A2 = parameters of the critical N curve

# tau_b and tau_n = 1/residual variances for biomass and nitrogen content observations

# Mu_Bmax,Prec_Bmax = parameters defining the between-date variability of Bmax

# Mu_S,Prec_S = parameters defining the between-date variability of S

##### Model 1 #####

Q<-length(Date)
K<-length(unique(Date))

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
			Mu_Bmax~dnorm(6,0.1)
			Mu_S~dnorm(0,0.1)
			A1~dunif(2,6)
			A2~dunif(0,0.5)

			#Informative prior C3
			#A1~dnorm(4.89,7.72)T(4,5.5)
			#ZA2~dbeta(2.12,2.12)
			#A2=(0.4-0.3)*ZA2+0.3
			#ZMu_Bmax~dbeta(2.31,2.31)
			#Mu_Bmax=(15-1)*ZMu_Bmax+1
			#Mu_S~dnorm(0,0.1)

			#Informative prior C4
			#ZA1~dbeta(2.03,1.5)
			#A1=(4-3)*ZA1+3
			#ZA2~dbeta(2.12,2.12)
			#A2=(0.4-0.3)*ZA2+0.3
			#ZMu_Bmax~dbeta(2.31,2.31)
			#Mu_Bmax=(15-1)*ZMu_Bmax+1
			#Mu_S~dnorm(0,0.1)

			Prec_Bmax~dgamma(0.001,0.001)
			Prec_S~dgamma(0.001,0.001)
			tau_b~dgamma(0.001,0.001)
			tau_n~dgamma(0.001,0.001)
			
}
"

writeLines(modelstring,con="model.txt")

model<-jags.model('model.txt',data=list('W'=W,'N'=N,'Date'=Date,'Q'=Q,'K'=K),
	n.chains=3, n.adapt=10000)

##### Model 2 #####

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

			#Informative prior C3
			#A1~dnorm(4.89,7.72)T(4,5.5)
			#ZA2~dbeta(2.12,2.12)
			#A2=(0.4-0.3)*ZA2+0.3
			#ZMu_Bmax~dbeta(2.31,2.31)
			#Mu_Bmax=(15-1)*ZMu_Bmax+1
			#Mu_S~dnorm(0,0.1)

			#Informative prior C4
			#ZA1~dbeta(2.03,1.5)
			#A1=(4-3)*ZA1+3
			#ZA2~dbeta(2.12,2.12)
			#A2=(0.4-0.3)*ZA2+0.3
			#ZMu_Bmax~dbeta(2.31,2.31)
			#Mu_Bmax=(15-1)*ZMu_Bmax+1
			#Mu_S~dnorm(0,0.1)

			Prec_Bmax~dgamma(0.001,0.001)
			Prec_S~dgamma(0.001,0.001)
			tau_b~dgamma(0.001,0.001)
			tau_n~dgamma(0.001,0.001)
			tau_t~dgamma(0.001,0.001)
}
"
writeLines(modelstring,con="model.txt")
model<-jags.model('model.txt',data=list('W'=W,'N'=N,'Date'=Date,'Q'=Q,'K'=K),
                  n.chains=3, n.adapt=10000)

##### END #####


