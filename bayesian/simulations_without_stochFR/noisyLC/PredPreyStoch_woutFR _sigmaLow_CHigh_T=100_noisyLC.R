### F. Barraquand 21/04/2015 - Code for analyzing noisy predator-prey system (incl. noise on the functional response) 
### Predator prey-only, but no fit of functional response (28/04/2015)

# Looks like the FR parameters are the most difficult to estimate indirectly, which is interesting for a 
# joint model fit. 

# small previous problem with Beverton-Holt formulation. Now Getz formulation, better. 

rm(list=ls())
graphics.off()

library("R2jags")      # Load R2jags package


### Parameters for simulation of Hassell model

n.years<-100  	# Number of years
N1<-1			# Initial pop size
P1<-0.1
K<-1			# threshold dd 
beta<-1			# density-dependence exponent
rmax_V<-2			# Max AVERAGE growth rate (thus not a true max...)
rmax_P<-0.5
sigma2.proc<-0.05		# Process sigma on the log-scale

C<-15
D<-0.25
Q<-10

### Simulation of data
set.seed(41) 
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
plot(1:n.years,N,type="b")
lines(1:n.years,P,type="b",col="blue")

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
    C~dgamma(.01,.01) # 
    D~dgamma(.01,.01)
    # check model Leslie to see how she specified priors...

    # Likelihood
    # state process

    for (t in 1:(T-1)){    
    
    
    logN[t+1] ~ dnorm(logNupdate[t],tau_V)
    logNupdate[t] <- logN[t] + r_V -log(1+N[t]/K_V) -C*P[t]/(D+N[t])
    N[t]<-exp(logN[t])
    P[t]<-exp(logP[t])
    # for some reason, log(1+(exp(r_V)-1)*N[t]/K_V) was not working
    logP[t+1]~ dnorm(logPupdate[t],tau_P)
    logPupdate[t] <- logP[t] + r_P - log(1+P[t]*Q/N[t])  
    
    }
    
    }
    ",fill=TRUE)
sink()

# 
# # Initial values
# inits <- function () {
#   list(sigma_V=runif(1,0.1,2), sigma_P=runif(1,0.1,2), r_V=runif(1,0.1,2),r_P=runif(1,0.1,2), K_V=runif(1,0.2,10), Q=runif(1,0,5),C=runif(1,10,100),D=runif(1,0.01,0.1))}


# Initial values
inits <- function () {
  list(sigma_V=runif(1,0.01,1), sigma_P=runif(1,0.01,1), r_V=runif(1,0.1,3),r_P=runif(1,0.1,1), K_V=runif(1,0.2,10), Q=runif(1,5,15),C=runif(1,1,10),D=runif(1,0.5,10))}


# Parameters monitored
parameters<-c("r_V","K_V","r_P","Q","sigma2_V","sigma2_P","C","D")

# MCMC settings
nc <- 3 #number of chains
nb <- 14000 # “burn in”
#ni <- 14000# “number of iterations” # that's for a symmetric distrib...
ni<-34000
nt <- 10 # “thinning”

# run model
out <- jags(jags.data, inits, parameters, "ssm.predprey1.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, working.directory = getwd())
print(out, dig = 2)

### Very reasonable estimates
# Inference for Bugs model at "ssm.predprey1.txt", fit using jags,
# 3 chains, each with 34000 iterations (first 14000 discarded), n.thin = 10
# n.sims = 6000 iterations saved
# mu.vect sd.vect   2.5%    25%    50%    75%  97.5% Rhat n.eff
# C          13.81    0.67  12.53  13.36  13.79  14.25  15.19    1  6000
# D           0.21    0.02   0.17   0.19   0.21   0.22   0.26    1  3000
# K_V         0.88    0.14   0.64   0.78   0.87   0.97   1.19    1  2300
# Q          10.64    1.18   8.50   9.82  10.60  11.41  13.09    1  4300
# r_P         0.52    0.04   0.44   0.49   0.52   0.55   0.61    1  6000
# r_V         2.09    0.12   1.86   2.01   2.09   2.18   2.34    1  1700
# sigma2_P    0.05    0.01   0.03   0.04   0.05   0.05   0.06    1  6000
# sigma2_V    0.05    0.01   0.04   0.05   0.05   0.06   0.07    1  4900
# deviance  -23.80    4.15 -29.77 -26.83 -24.44 -21.56 -14.10    1  2100

# check some dynamical aspects of the data first
par(mfrow=c(2,2))
plot(N,type="b")
plot(P,type="b")
plot(log(N),log(P))
plot(N,C*N/(D+N))

### Diagnostics
library(mcmcplots)
### Trace plots
png(file = "TracePlot_T=100_noisyLC.png")
traplot(out,c("r_V","K_V","r_P","Q","sigma2_V","sigma2_P","C","D"))
dev.off()
### plot densities
denplot(out,c("r_V","K_V","r_P","Q","sigma2_V","sigma2_P","C","D"))
