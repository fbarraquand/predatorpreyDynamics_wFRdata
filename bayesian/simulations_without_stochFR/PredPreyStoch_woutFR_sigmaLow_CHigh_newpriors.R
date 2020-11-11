### F. Barraquand 21/04/2015 - Code for analyzing noisy predator-prey system (incl. noise on the functional response) 
### Predator prey-only, but no fit of functional response (28/04/2015)
### Correction predator density multiplication -- 26/12/2018
### New priors for C and D -- 26/12/2018

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
    logC ~ dnorm(0,0.5) #C~dgamma(.01,.01)
    logD ~ dnorm(0,0.5) #D~gamma(.01,.01)
    C<-exp(logC)
    D<-exp(logD)

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
parameters<-c("r_V","K_V","r_P","Q","sigma2_V","sigma2_P","logC","logD") #"logN","logP"

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
traplot(out,c("r_V","K_V","r_P","Q","sigma2_V","sigma2_P","logC","logD"))
### plot densities
denplot(out,c("r_V","K_V","r_P","Q","sigma2_V","sigma2_P","logC","logD"))

# Inference for Bugs model at "ssm.predprey1.txt", fit using jags,
# 3 chains, each with 20000 iterations (first 10000 discarded), n.thin = 10
# n.sims = 3000 iterations saved
# mu.vect sd.vect   2.5%    25%    50%    75%  97.5% Rhat n.eff
# K_V          1.37    0.68   0.33   0.89   1.26   1.74   2.99 1.06    55
# Q           12.37    2.89   7.63  10.29  12.09  14.07  19.00 1.01   530
# logC        -0.21    0.81  -2.04  -0.72  -0.11   0.40   1.12 1.00   690
# logD         0.04    0.97  -1.91  -0.60   0.04   0.70   1.93 1.00   710
# r_P          0.56    0.10   0.38   0.49   0.55   0.62   0.77 1.01   600
# r_V          1.80    0.46   1.10   1.48   1.75   2.03   2.95 1.05    56
# sigma2_P     0.04    0.01   0.03   0.04   0.04   0.04   0.05 1.00  3000
# sigma2_V     0.06    0.01   0.04   0.05   0.06   0.06   0.08 1.00  2700
# deviance   -29.59    3.79 -34.79 -32.35 -30.36 -27.60 -20.38 1.00   620
out$BUGSoutput$sims.list$logC
out$BUGSoutput$sims.list$logD

# comparing mean values of C and D
C
D
exp(out$BUGSoutput$mean$logC) ## Not the same thing - but is this because I output the median?
exp(out$BUGSoutput$mean$logD)


### --- WITH dnorm(1,0.1) priors ---
# Inference for Bugs model at "ssm.predprey1.txt", fit using jags,
# 3 chains, each with 20000 iterations (first 10000 discarded), n.thin = 10
# n.sims = 3000 iterations saved
# mu.vect sd.vect   2.5%    25%    50%    75%  97.5% Rhat n.eff
# K_V          1.49    0.65   0.52   1.01   1.43   1.87   3.00 1.02   210
# Q           12.47    2.89   7.57  10.37  12.28  14.23  18.73 1.00  2500
# logC        -1.45    2.35  -6.71  -2.85  -1.11   0.20   2.50 1.00  3000
# logD         0.45    3.15  -6.00  -1.62   0.50   2.56   6.65 1.00  3000
# > C
# [1] 2.5
# > D
# [1] 1
# > exp(out$BUGSoutput$mean$logC) ## Not the same thing
# [1] 0.2341679
# > exp(out$BUGSoutput$mean$logD)
# [1] 1.57264

### --- with dnorm(1,0.5)
### C estimated to 0.65 + 0.81^2/2 =0.97805 <2.5 => NO
