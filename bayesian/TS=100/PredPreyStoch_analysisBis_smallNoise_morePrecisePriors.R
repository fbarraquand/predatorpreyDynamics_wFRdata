### F. Barraquand 21/04/2015 - Code for analyzing noisy predator-prey system (incl. noise on the functional response)
### Updated 01/06/2017. 
### Case with "small noise" on the functional response
### FB 03/01/2019 Added plots for the paired posterior distributions
### FB 14/01/2019 Corrected FR to nondelayed version
### FB 31/01/2019 -- more precise priors on (C,D)

### Log of previous edits to other versions ######################################################################################

### Edited 25/05/2015 Add of more diagnostic plots. Problems with the functional response model, I should probably
# analyze it separarely to see if all components can be identified - this looks quite unclear. 
###### Edit again 25/05. actually no pb with the FR model, there's just that with the amount of noise, the Beta model for the FR
# is probably not identifiable and the FR half-sat constant is estimated at zero, which means a constant amount. 
# There was by the way initially an error in the simulation, FR[t+1]<-C*(1-B)*N[t]/(D+N[t]) was FR[t+1]<-C*(1-B)*N[1]/(D+N[1]) 
# So the FR was rightfully more constant!! Corrected now. 
# we can probably force a handling-time or half-sat value but even with an informative prior it will be estimated low. 
# Keep in mind for the gyr and ptarmigan

### Predator prey-only, fit of functional response (28/04/2015) but direct estimation (same model simulated as fitted). 
# small previous problem with Beverton-Holt formulation. Now Getz formulation, better. 

### For the functional response, the beta model is perhaps good for simulations, but I feel like we need something else
# to fit the model. Perhaps a logNormal. 

##################################################################################################################################

rm(list=ls())
graphics.off()

library("R2jags")      # Load R2jags package

### Parameters for simulation 

n.years<-100  	# Number of years 
N1<-1			# Initial pop size
P1<-0.1
K<-1			# threshold dd 
beta<-1			# density-dependence exponent
rmax_V<-2			# Max AVERAGE growth rate (thus not a true max...)
rmax_P<-0.5
sigma2.proc<-0.05		# worked well with 0.005 as well
# Process sigma on the log-scale, use the Peretti et al. value. 0.005


### FR and predator parameters
C<-2.5
D<-1
Q<-10

### Simulation of data
#set.seed(42) 
set.seed(41)
#set.seed(40)

y<-N<-P<-FR<-numeric(n.years)
N[1]<-N1
P[1]<-P1
FR[1]<-C*N[1]/(D+N[1])

rV<-rnorm(n.years-1,rmax_V,sqrt(sigma2.proc))
rP<-rnorm(n.years-1,rmax_P,sqrt(sigma2.proc))
FRnoise<-rnorm(n.years-1,0,sqrt(sigma2.proc))

for (t in 1:(n.years-1)){

  N[t+1]<-N[t]*(exp(rV[t])/(1+(N[t]/K)^beta))*exp(-FR[t]*P[t]/N[t])
  P[t+1]<-P[t]*exp(rP[t])/(1+P[t]*Q/N[t])
  FR[t+1]<-(C*N[t+1]/(D+N[t+1])) + FRnoise[t+1] #after for update
}
## Plotting time series of abundances and FR
par(mfrow=c(2,2))
plot(1:n.years,N,type="b")
plot(1:n.years,P,type="b")
#curve(dbeta(x,a,b),from=0, to=1)
plot(N,FR)
plot(log(N),log(P))
    
#seq(0,1,0.01)
# Bundle data
jags.data <- list(T=n.years,logN=log(N),logP=log(P),FR=FR)

sink("predprey.txt")
cat("
    model {
    
    # Priors and constraints
    logN[1] ~ dnorm(0,0.01) # Prior on initial pop size on the log scale
    logP[1] ~ dnorm(0,0.01) # Prior on initial pop size on the log scale

    # Priors prey population dynamics
    r_V ~ dnorm(1,0.001) # below the truth, rather flat prior
    K_V ~ dunif(0.2,10)
    sigma_V ~ dunif(0.01,5) # rather vague 
    sigma2_V<-pow(sigma_V, 2)
    tau_V<-pow(sigma_V,-2)

 
     #Priors predator population dynamics
    Q ~ dgamma(0.1,0.1)
    r_P ~ dnorm(1,0.1)
    sigma_P ~ dunif(0.01,2) # rather vague 
    sigma2_P<-pow(sigma_P, 2)
    tau_P<-pow(sigma_P,-2)
    
    #Priors predation parameters 
    tau_FR ~ dgamma(.01,.01)
    C~dgamma(2,1) #C~dgamma(2,0.1) # beware, curve(dgamma(x,2,1),xlim=c(0,10)) in R
    D~dgamma(2,1) #D~dgamma(2,0.1)
    # check model Leslie to see how she specified priors...

    # Likelihood
    # state process

    for (t in 1:(T-1)){        

    FRUpdate[t] <- C*N[t]/(D+N[t]) #functional response equation, including noise
    FR[t] ~  dnorm(FRUpdate[t],tau_FR) #small trick to use FR data

    logNupdate[t] <- logN[t] + r_V -log(1+N[t]/K_V) -FR[t]*exp(logP[t])/N[t]
    logN[t+1] ~ dnorm(logNupdate[t],tau_V)
    N[t]<-exp(logN[t])
    # for some reason, log(1+(exp(r_V)-1)*N[t]/K_V) was not working

    logP[t+1]~ dnorm(logPupdate[t],tau_P)
    logPupdate[t] <- logP[t] + r_P - log(1+exp(logP[t])*Q/exp(logN[t]) )  
    
    }
    
    }
    ",fill=TRUE)
sink()


# Initial values
inits <- function () {
  list(sigma_V=runif(1,0.1,2), sigma_P=runif(1,0.1,2), r_V=runif(1,0.1,2),r_P=runif(1,0.1,2), K_V=runif(1,0.2,8), Q=runif(1,0,5),tau_FR=runif(1,1,10),C=runif(1,10,100),D=runif(1,0.01,0.1))}

# Parameters monitored
#parameters<-c("r_V","K_V","r_P","Q","sigma2_V","sigma2_P","a","b","C","D","logNupdate","logPupdate","FR")
parameters<-c("r_V","K_V","r_P","Q","sigma2_V","sigma2_P","tau_FR","C","D")


# MCMC settings
nc <- 3 #number of chains
nb <- 14000 # “burn in”
#ni <- 14000# “number of iterations” # that's for a symmetric distrib...
ni<-34000
nt <- 10 # “thinning”

# run model
out <- jags(jags.data, inits, parameters, "predprey.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, working.directory = getwd())
print(out, dig = 2)

CEb<-out$BUGSoutput$mean$C
DEb<-out$BUGSoutput$mean$D

# logN1<-out$BUGSoutput$mean$logNupdate
# logP1<-out$BUGSoutput$mean$logPupdate
# Better to output logN!

# 
# par(mfrow=c(2,2))
# plot(1:(n.years-1),logN1,type="o")
# plot(1:(n.years-1),logP1,type="o")

### Fit functional response
par(mfrow=c(1,1))

fr_fit<-nls(FR~CE*N/(DE+N),start=list(CE=1,DE=1))
CE<-coef(fr_fit)[1]
DE<-coef(fr_fit)[2]
plot(N,FR,ylim=c(0,max(FR,na.rm=T)))
curve(CE*x/(DE+x),add=TRUE)
curve(CEb*x/(DEb+x),col="blue",add=TRUE)

### Now try to fit a model without the FR data. 
sink("predprey_without_sepFR.txt")
cat("
    model {
    
    # Priors and constraints
    logN[1] ~ dnorm(0,0.01) # Prior on initial pop size on the log scale
    logP[1] ~ dnorm(0,0.01) # Prior on initial pop size on the log scale
    
    # Priors prey population dynamics
    r_V ~ dnorm(1,0.001) # below the truth, rather flat prior
    K_V ~ dunif(0.2,10)
    sigma_V ~ dunif(0.01,5) # rather vague 
    sigma2_V<-pow(sigma_V, 2)
    tau_V<-pow(sigma_V,-2)
    
    
    #Priors predator population dynamics
    Q ~ dgamma(0.1,0.1)
    r_P ~ dnorm(1,0.1)
    sigma_P ~ dunif(0.01,2) # rather vague 
    sigma2_P<-pow(sigma_P, 2)
    tau_P<-pow(sigma_P,-2)
    
    #Priors predation parameters 
    
    C~dgamma(2,1) # C~dgamma(2,0.1)
    D~dgamma(2,1) # C~dgamma(2,0.1)
    
    # check model Leslie to see how she specified priors...
    
    # Likelihood
    # state process
    
    for (t in 1:(T-1)){        

    
    logNupdate[t] <- logN[t] + r_V -log(1+N[t]/K_V) -C*exp(logP[t])/(D+N[t])
    logN[t+1] ~ dnorm(logNupdate[t],tau_V)
    N[t]<-exp(logN[t])
    
    # for some reason, log(1+(exp(r_V)-1)*N[t]/K_V) was not working
    logP[t+1]~ dnorm(logPupdate[t],tau_P)
    logPupdate[t] <- logP[t] + r_P - log(1+exp(logP[t])*Q/exp(logN[t]) )  
    
    }
    
    }
    ",fill=TRUE)
sink()


# Initial values
inits <- function () {
  list(sigma_V=runif(1,0.1,2), sigma_P=runif(1,0.1,2), r_V=runif(1,0.1,2),r_P=runif(1,0.1,2), K_V=runif(1,0.2,8), Q=runif(1,0,5),C=runif(1,10,100),D=runif(1,0.01,0.1))}

# Parameters monitored
#parameters<-c("r_V","K_V","r_P","Q","sigma2_V","sigma2_P","a","b","C","D","logNupdate","logPupdate","FR")
parameters<-c("r_V","K_V","r_P","Q","sigma2_V","sigma2_P","C","D")


# MCMC settings
nc <- 3 #number of chains
nb <- 14000 # “burn in”
#ni <- 14000# “number of iterations” # that's for a symmetric distrib...
ni<-34000
nt <- 10 # “thinning”

# run model
out2 <- jags(jags.data, inits, parameters, "predprey_without_sepFR.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, working.directory = getwd())
print(out2, dig = 2)
print(out,dig=2) # 

# Inference for Bugs model at "predprey.txt", fit using jags,
# 3 chains, each with 34000 iterations (first 14000 discarded), n.thin = 10
# n.sims = 6000 iterations saved
# mu.vect sd.vect   2.5%    25%    50%    75%  97.5% Rhat n.eff
# C           2.58    0.12   2.35   2.49   2.57   2.66   2.84 1.00  4300
# D           1.24    0.32   0.66   1.01   1.22   1.44   1.91 1.00  4000
# K_V         1.65    0.62   0.61   1.21   1.58   2.04   3.04 1.04    63
# Q          11.16    2.89   6.43   9.09  10.88  12.87  17.65 1.01   420
# r_P         0.54    0.10   0.35   0.47   0.53   0.61   0.76 1.01   420
# r_V         1.71    0.32   1.20   1.48   1.67   1.87   2.46 1.04    61
# sigma2_P    0.05    0.01   0.03   0.04   0.05   0.05   0.06 1.00  6000
# sigma2_V    0.05    0.01   0.04   0.05   0.05   0.06   0.07 1.00  4000
# tau_FR     17.48    2.49  12.98  15.71  17.42  19.09  22.52 1.00  3300
# deviance  -23.39    4.23 -29.80 -26.53 -23.98 -20.95 -13.52 1.00   740


# For C ~ dgamma(2,0.1), D ~ dgamma(2,0.1)
# Inference for Bugs model at "predprey_without_sepFR.txt", fit using jags,
# 3 chains, each with 34000 iterations (first 14000 discarded), n.thin = 10
# n.sims = 6000 iterations saved
# mu.vect sd.vect   2.5%    25%    50%    75%  97.5% Rhat n.eff
# C          13.29    7.53   2.84   7.78  11.89  17.18  32.23 1.00   980
# D          26.18   14.63   5.90  15.76  23.44  33.56  61.68 1.00  1400
# K_V         2.14    0.81   0.80   1.56   2.04   2.63   3.95 1.01   190
# Q          11.17    2.91   6.36   9.11  10.86  12.85  17.67 1.00   580
# r_P         0.54    0.11   0.34   0.47   0.53   0.61   0.76 1.00   650
# r_V         1.53    0.30   1.03   1.31   1.50   1.70   2.25 1.01   230
# sigma2_P    0.05    0.01   0.03   0.04   0.05   0.05   0.06 1.00  2500
# sigma2_V    0.05    0.01   0.04   0.05   0.05   0.06   0.07 1.00  6000
# deviance  -21.84    3.77 -27.19 -24.57 -22.48 -19.84 -12.67 1.00  6000
# 
# 
# DIC info (using the rule, pD = var(deviance)/2)
# pD = 7.1 and DIC = -14.7
# DIC is an estimate of expected predictive error (lower deviance is better).
# > print(out,dig=2) #
# Inference for Bugs model at "predprey.txt", fit using jags,
# 3 chains, each with 34000 iterations (first 14000 discarded), n.thin = 10
# n.sims = 6000 iterations saved
# mu.vect sd.vect   2.5%    25%    50%    75%  97.5% Rhat n.eff
# C           2.64    0.14   2.39   2.54   2.63   2.73   2.95 1.01   400
# D           1.40    0.37   0.74   1.15   1.37   1.62   2.21 1.01   400
# K_V         1.63    0.64   0.52   1.18   1.58   2.00   3.03 1.07    65
# Q          11.22    2.99   6.13   9.14  10.90  12.97  17.86 1.01   450
# r_P         0.54    0.11   0.34   0.47   0.54   0.61   0.77 1.01   440
# r_V         1.73    0.35   1.21   1.49   1.67   1.91   2.61 1.05    75
# sigma2_P    0.05    0.01   0.03   0.04   0.05   0.05   0.06 1.00  4800
# sigma2_V    0.05    0.01   0.04   0.05   0.05   0.06   0.07 1.00  1700
# tau_FR     17.42    2.49  13.05  15.70  17.28  19.04  22.57 1.00  4200
# deviance  -22.87    4.55 -29.75 -26.21 -23.56 -20.31 -12.09 1.00   530



#### Results for   
#C~dgamma(1,1) # 
#D~dgamma(2,1) # 

# Behaviour of the posterior seems dominated by the prior
# > print(out2, dig = 2)
# Inference for Bugs model at "predprey_without_sepFR.txt", fit using jags,
# 3 chains, each with 34000 iterations (first 14000 discarded), n.thin = 10
# n.sims = 6000 iterations saved
# mu.vect sd.vect   2.5%    25%    50%    75%  97.5% Rhat n.eff
# C           1.35    0.98   0.05   0.58   1.17   1.93   3.69 1.00  4600
# D           1.96    1.38   0.24   0.93   1.65   2.65   5.55 1.00  1300
# K_V         1.85    0.79   0.57   1.35   1.74   2.24   3.74 1.02   200
# Q          11.47    3.03   6.43   9.37  11.19  13.16  18.15 1.00   670
# r_P         0.55    0.11   0.35   0.48   0.55   0.62   0.77 1.00   740
# r_V         1.57    0.37   0.98   1.33   1.52   1.73   2.50 1.02   250
# sigma2_P    0.05    0.01   0.03   0.04   0.05   0.05   0.06 1.00  6000
# sigma2_V    0.06    0.01   0.04   0.05   0.05   0.06   0.07 1.00  6000
# deviance  -20.52    3.92 -26.24 -23.31 -21.20 -18.39 -10.84 1.00  1100
# 
# For each parameter, n.eff is a crude measure of effective sample size,
# and Rhat is the potential scale reduction factor (at convergence, Rhat=1).
# 
# DIC info (using the rule, pD = var(deviance)/2)
# pD = 7.7 and DIC = -12.8
# DIC is an estimate of expected predictive error (lower deviance is better).
# > print(out,dig=2) #
# Inference for Bugs model at "predprey.txt", fit using jags,
# 3 chains, each with 34000 iterations (first 14000 discarded), n.thin = 10
# n.sims = 6000 iterations saved
# mu.vect sd.vect   2.5%    25%    50%    75%  97.5% Rhat n.eff
# C           2.58    0.13   2.35   2.49   2.57   2.66   2.85 1.01   410
# D           1.24    0.33   0.67   1.01   1.21   1.44   1.95 1.00   480
# K_V         1.66    0.59   0.60   1.24   1.61   2.01   2.97 1.04   160
# Q          11.39    2.93   6.41   9.33  11.20  13.16  17.90 1.00  1800
# r_P         0.55    0.11   0.35   0.47   0.55   0.62   0.77 1.00  2000
# r_V         1.70    0.31   1.22   1.49   1.65   1.86   2.49 1.03   190
# sigma2_P    0.05    0.01   0.03   0.04   0.05   0.05   0.06 1.00  5900
# sigma2_V    0.05    0.01   0.04   0.05   0.05   0.06   0.07 1.00  6000
# tau_FR     17.43    2.51  12.84  15.70  17.27  19.07  22.76 1.00  6000
# deviance  -23.48    4.23 -29.79 -26.58 -24.17 -20.95 -13.55 1.00   730
# 
# For each parameter, n.eff is a crude measure of effective sample size,
# and Rhat is the potential scale reduction factor (at convergence, Rhat=1).
# 
# DIC info (using the rule, pD = var(deviance)/2)
# pD = 8.9 and DIC = -14.5
# DIC is an estimate of expected predictive error (lower deviance is better).

# 01/06/2017 // Olivier says DIC is crap in this (and other?) context, avoid this...  
#plot(as.mcmc(out2)) 
#plot(as.mcmc(out)) 

#save.image('/media/frederic/DATA/Simuls_wOlivier/predatorpreyFRdata/simTS100.RData')

### plot densities
library(mcmcplots)
denplot(out,c("r_V","K_V","r_P","Q","sigma2_V","sigma2_P","C","D"))
denplot(out2,c("r_V","K_V","r_P","Q","sigma2_V","sigma2_P","C","D"))

### Trace plots
traplot(out,c("r_V","K_V","r_P","Q","sigma2_V","sigma2_P","C","D"))
traplot(out2,c("r_V","K_V","r_P","Q","sigma2_V","sigma2_P","C","D"))


