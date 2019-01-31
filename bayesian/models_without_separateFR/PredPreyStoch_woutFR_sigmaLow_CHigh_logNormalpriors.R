### F. Barraquand 21/04/2015 - Code for analyzing noisy predator-prey system (incl. noise on the functional response) 
### Predator prey-only, but no fit of functional response (28/04/2015)
### Correction predator density multiplication -- 26/12/2018
### New priors for C and D -- 26/12/2018
### Lognormal priors with dlnorm -- 31/01/2019

# Looks like the FR parameters are the most difficult to estimate indirectly, which is interesting for a 
# joint model fit. 

# small previous problem with Beverton-Holt formulation. Now Getz formulation, better. 

rm(list=ls())
graphics.off()

library("R2jags")      # Load R2jags package

### Parameters for simulation of Hassell model

n.years<-100  	# Number of years - 25 first, perhaps use 50 or 100 / worked very well with almost no process error on the real scale
N1<-1			# Initial pop size
P1<-0.1
K<-1			# threshold dd 
beta<-1			# density-dependence exponent
rmax_V<-2			# Max AVERAGE growth rate (thus not a true max...)
rmax_P<-0.5
sigma2.proc<-0.05 #1		# Process sigma on the log-scale

C<-2.5 #0.5
D<-1
Q<-10

### Simulation of data
set.seed(42) 
y<-N<-P<-numeric(n.years)
N[1]<-N1
P[1]<-P1

rV<-rnorm(n.years-1,rmax_V,sqrt(sigma2.proc))
rP<-rnorm(n.years-1,rmax_P,sqrt(sigma2.proc))
for (t in 1:(n.years-1)){
  
  N[t+1]<-N[t]*(exp(rV[t])/(1+(N[t]/K)^beta))*exp(-C*P[t]/(D+N[t]))
  P[t+1]<-P[t]*exp(rP[t])/(1+P[t]*Q/N[t])
}
## Plotting time series
par(mfrow=c(2,1))
plot(1:n.years,N,type="b")
plot(1:n.years,P,type="b")

# Bundle data
jags.data <- list(T=n.years,logN=log(N),logP=log(P))

sink("ssm.predprey1.txt")
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
    C ~ dlnorm(0,0.5) #C~dgamma(.01,.01)
    D ~ dlnorm(0,0.5) #D~gamma(.01,.01)


    # Likelihood
    # state process

    for (t in 1:(T-1)){    
    
    logN[t+1] ~ dnorm(logNupdate[t],tau_V)
    logNupdate[t] <- logN[t] + r_V -log(1+N[t]/K_V) -C*exp(logP[t])/(D+N[t])
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
  list(sigma_V=runif(1,0.1,2), sigma_P=runif(1,0.1,2), r_V=runif(1,0.1,2),r_P=runif(1,0.1,2), K_V=runif(1,0.2,10), Q=runif(1,0,5),logC=rnorm(1,0,1),logD=rnorm(1,0,1))}
#C=runif(1,10,100),D=runif(1,0.01,0.1)

# Parameters monitored
parameters<-c("r_V","K_V","r_P","Q","sigma2_V","sigma2_P","C","D") #"logN","logP"

# MCMC settings
nc <- 3 #number of chains
nb <- 10000 # “burn in”
#ni <- 14000# “number of iterations” # that's for a symmetric distrib...
ni<-20000 #34000
nt <- 10 # “thinning”

# run model
out <- jags(jags.data, inits, parameters, "ssm.predprey1.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, working.directory = getwd())
print(out, dig = 2)

### Diagnostics
library(mcmcplots)
### Trace plots
traplot(out,c("r_V","K_V","r_P","Q","sigma2_V","sigma2_P","C","D"))
### plot densities
denplot(out,c("r_V","K_V","r_P","Q","sigma2_V","sigma2_P","C","D"))

# Inference for Bugs model at "ssm.predprey1.txt", fit using jags,
# 3 chains, each with 20000 iterations (first 10000 discarded), n.thin = 10
# n.sims = 3000 iterations saved
# mu.vect sd.vect   2.5%    25%    50%    75%  97.5% Rhat n.eff
# C           1.07    1.07   0.06   0.34   0.74   1.47   3.84 1.00  1100
# D           2.50    4.01   0.07   0.44   1.08   2.76  13.93 1.00  1200
# K_V         1.39    0.65   0.27   0.95   1.34   1.76   2.83 1.12    41
# Q          12.24    2.78   7.43  10.32  11.98  13.97  18.45 1.01   420
# r_P         0.55    0.09   0.38   0.49   0.55   0.62   0.75 1.01   380
# r_V         1.78    0.48   1.12   1.47   1.69   1.98   3.14 1.08    51
# sigma2_P    0.04    0.01   0.03   0.04   0.04   0.04   0.05 1.00  2400
# sigma2_V    0.06    0.01   0.04   0.05   0.06   0.06   0.08 1.00   810
# deviance  -29.55    3.64 -34.63 -32.27 -30.16 -27.54 -20.89 1.00   780

