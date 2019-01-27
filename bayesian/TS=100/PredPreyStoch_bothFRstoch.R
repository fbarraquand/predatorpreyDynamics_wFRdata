### F. Barraquand 21/04/2015 - Code for analyzing noisy predator-prey system (incl. noise on the functional response)
### Updated 01/06/2017. 
### Case with "small noise" on the functional response
### FB 03/01/2019 Added plots for the paired posterior distributions
### FB 14/01/2019 Corrected FR to nondelayed version
### FB 27/01/2019 Both FR are stochastic here, only the data differs

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
    C~dgamma(.01,.01) # uninformative priors OK for that one
    D~dgamma(.01,.01)
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


# Output summary statistics
jags.sum<-out$BUGSoutput$summary
write.table(x=jags.sum,file="JAGSsummary_PredPreyStoch_analysisBis.txt")
# MCMC Output
#pdf("Output_MCMC__PredPreyStoch_analysisBis.pdf")
out.mcmc<-as.mcmc(out)
plot(out.mcmc) 
#dev.off()

CEb<-out$BUGSoutput$mean$C
DEb<-out$BUGSoutput$mean$D


## Removes FR data
jags.data <- list(T=n.years,logN=log(N),logP=log(P))
  
# run model
out2 <- jags(jags.data, inits, parameters, "predprey.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, working.directory = getwd())
print(out2, dig = 2)
print(out,dig=2) # 

# 
# > print(out2, dig = 2)
# Inference for Bugs model at "predprey.txt", fit using jags,
# 3 chains, each with 34000 iterations (first 14000 discarded), n.thin = 10
# n.sims = 6000 iterations saved
# mu.vect sd.vect   2.5%    25%    50%    75%  97.5% Rhat n.eff
# C           0.00    0.01   0.00   0.00   0.00   0.00   0.01 1.01   160
# D           0.00    0.00   0.00   0.00   0.00   0.00   0.00 1.22    17
# K_V         1.91    0.82   0.67   1.29   1.81   2.39   3.78 1.06    37
# Q          11.21    2.92   6.34   9.15  10.92  12.94  17.81 1.00  1200
# r_P         0.54    0.11   0.35   0.47   0.54   0.61   0.77 1.00  1200
# r_V         1.48    0.35   0.93   1.23   1.43   1.70   2.25 1.06    37
# sigma2_P    0.05    0.01   0.03   0.04   0.05   0.05   0.06 1.00  4800
# sigma2_V    0.05    0.01   0.03   0.05   0.05   0.06   0.07 1.01   850
# tau_FR     17.62   36.92   0.17   1.08   4.11  15.65 124.88 1.01   420
# deviance  -25.29   15.42 -67.66 -26.17 -21.64 -18.31 -10.76 1.04   440
# 
# For each parameter, n.eff is a crude measure of effective sample size,
# and Rhat is the potential scale reduction factor (at convergence, Rhat=1).
# 
# DIC info (using the rule, pD = var(deviance)/2)
# pD = 118.4 and DIC = 93.1
# DIC is an estimate of expected predictive error (lower deviance is better).

# ----------------------------------------------------------------------- #

# > print(out,dig=2) #
# Inference for Bugs model at "predprey.txt", fit using jags,
# 3 chains, each with 34000 iterations (first 14000 discarded), n.thin = 10
# n.sims = 6000 iterations saved
# mu.vect sd.vect   2.5%    25%    50%    75%  97.5% Rhat n.eff
# C           2.56    0.14   2.32   2.46   2.55   2.65   2.86 1.00   680
# D           1.19    0.37   0.57   0.93   1.16   1.41   1.98 1.00   960
# K_V         1.71    0.61   0.75   1.28   1.63   2.07   3.13 1.01   210
# Q          11.26    2.93   6.46   9.21  10.95  13.05  17.58 1.01   310
# r_P         0.54    0.11   0.35   0.47   0.54   0.61   0.76 1.01   300
# r_V         1.67    0.28   1.18   1.46   1.65   1.84   2.28 1.01   200
# sigma2_P    0.05    0.01   0.03   0.04   0.05   0.05   0.06 1.00  6000
# sigma2_V    0.05    0.01   0.04   0.05   0.05   0.06   0.07 1.00  3900
# tau_FR     17.41    2.52  12.75  15.69  17.32  18.99  22.60 1.00  5800
# deviance  -23.10    4.48 -29.83 -26.40 -23.74 -20.45 -12.76 1.00   520
# 
# For each parameter, n.eff is a crude measure of effective sample size,
# and Rhat is the potential scale reduction factor (at convergence, Rhat=1).
# 
# DIC info (using the rule, pD = var(deviance)/2)
# pD = 10.0 and DIC = -13.1
# DIC is an estimate of expected predictive error (lower deviance is better).

