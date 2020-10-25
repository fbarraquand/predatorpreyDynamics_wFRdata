### F. Barraquand 21/04/2015 - Code for analyzing noisy predator-prey system (incl. noise on the functional response) 
### Predator prey-only, but no fit of functional response (28/04/2015)

# Looks like the FR parameters are the most difficult to estimate indirectly, which is interesting for a 
# joint model fit. 

## Partially fixed by the priors here. 


rm(list=ls())
graphics.off()

library("R2jags")      # Load R2jags package
#library("modeest")

### Parameters for simulation of Hassell model

n.years<-100  	# Number of years - 25 first, perhaps use 50 or 100 / worked very well with almost no process error on the real scale
N1<-1			# Initial pop size
P1<-0.1
K<-1			# threshold dd 
beta<-1			# density-dependence exponent
rmax_V<-2			# Max AVERAGE growth rate (thus not a true max...)
rmax_P<-0.5
sigma2.proc<-0.05		# Process sigma on the log-scale

C<-2.5
D<-1
Q<-10

### Simulation of data
#set.seed(42) 
set.seed(43)

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
lines(1:n.years,P,type="b")

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
    # C~dgamma(.1,.1) # the lower the vaguer // E(C) = a/b, V(C) = a/b^2
    C~dgamma(2,1)
    D~dgamma(2,1)

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


# Initial values
inits <- function () {
  list(sigma_V=runif(1,0.1,2), sigma_P=runif(1,0.1,2), r_V=runif(1,0.1,2),r_P=runif(1,0.1,2), K_V=runif(1,0.2,10), Q=runif(1,0,5),C=runif(1,1,10),D=runif(1,0.01,0.1))}


# Parameters monitored
parameters<-c("r_V","K_V","r_P","Q","sigma2_V","sigma2_P","C","D") #,"logN","logP"

# MCMC settings
nc <- 3 #number of chains
nb <- 14000 # “burn in”
#ni <- 14000# “number of iterations” # that's for a symmetric distrib...
ni<-34000
nt <- 10 # “thinning”

# run model
out <- jags(jags.data, inits, parameters, "ssm.predprey1.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, working.directory = getwd())
print(out, dig = 2)

### plot densities
library(mcmcplots)
denplot(out,c("r_V","K_V","r_P","Q","sigma2_V","sigma2_P","C","D"))
### Trace plots
traplot(out,c("r_V","K_V","r_P","Q","sigma2_V","sigma2_P","C","D"))

### Old results 

# Inference for Bugs model at "ssm.predprey1.txt", fit using jags,
# 3 chains, each with 34000 iterations (first 14000 discarded), n.thin = 10
# n.sims = 6000 iterations saved
# mu.vect sd.vect   2.5%    25%    50%    75%  97.5% Rhat n.eff
# C            2.18    4.11   0.00   0.00   0.14   2.34  15.08 1.47     9
# D            1.24    1.95   0.00   0.00   0.00   2.27   6.47 1.19    15
# K_V          0.99    1.71   0.20   0.21   0.25   0.72   6.75 2.01     5
# Q           10.17    2.29   6.22   8.57   9.94  11.57  15.15 1.00  6000

### New results (with dgamma(.01,.01))

# Inference for Bugs model at "ssm.predprey1.txt", fit using jags,
# 3 chains, each with 34000 iterations (first 14000 discarded), n.thin = 10
# n.sims = 6000 iterations saved
# mu.vect sd.vect   2.5%    25%    50%    75%  97.5% Rhat n.eff
# C            0.00    0.01   0.00   0.00   0.00   0.00   0.01 1.31    10
# D            0.01    0.08   0.00   0.00   0.00   0.00   0.05 1.56     8
# K_V          1.50    0.70   0.52   0.99   1.39   1.88   3.19 1.01   190
# Q           12.27    2.81   7.53  10.33  12.00  14.01  18.38 1.00  4500

### New results (with dgamma(.1,.1))
# Inference for Bugs model at "ssm.predprey1.txt", fit using jags,
# 3 chains, each with 34000 iterations (first 14000 discarded), n.thin = 10
# n.sims = 6000 iterations saved
# mu.vect sd.vect   2.5%    25%    50%    75%  97.5% Rhat n.eff
# C            0.16    0.45   0.00   0.00   0.00   0.03   1.65 1.05    49
# D            0.60    1.69   0.00   0.00   0.00   0.26   6.03 1.21    20

### New results (with  C~dgamma(0.8,0.4), D~dgamma(.1,1))
# Inference for Bugs model at "ssm.predprey1.txt", fit using jags,
# 3 chains, each with 34000 iterations (first 14000 discarded), n.thin = 10
# n.sims = 6000 iterations saved
# mu.vect sd.vect   2.5%    25%    50%    75%  97.5% Rhat n.eff
# C            0.86    0.72   0.01   0.29   0.68   1.28   2.63 1.00  1300
# D            0.10    0.25   0.00   0.00   0.00   0.04   1.02 1.03    81
# K_V          1.32    0.63   0.38   0.87   1.24   1.65   2.84 1.01   420
# Q           12.31    2.87   7.49  10.34  12.03  13.99  18.75 1.00  3400

# With closer initial conditions but same priors
# Inference for Bugs model at "ssm.predprey1.txt", fit using jags,
# 3 chains, each with 34000 iterations (first 14000 discarded), n.thin = 10
# n.sims = 6000 iterations saved
# mu.vect sd.vect   2.5%    25%    50%    75%  97.5% Rhat n.eff
# C            1.01    0.81   0.02   0.35   0.84   1.50   2.91 1.01   460
# D            0.10    0.30   0.00   0.00   0.00   0.04   1.02 1.08    71
# K_V          1.09    0.67   0.25   0.56   0.95   1.45   2.76 1.03    73
# Q           12.45    2.80   7.70  10.49  12.17  14.16  18.65 1.00  5700

# New priors C~dgamma(1.25,0.5), D~dgamma(1,1)
# Sensible rule dgamma( ((mu^2)/sigma), (mu/sigma))
# Corresponding to a well-defined (though a bit flat) mode in the prior
# These are centered around the real values 
# But we can arbitrarily center them on 1 later on, which means dgamma(1,1) for both
# Inference for Bugs model at "ssm.predprey1.txt", fit using jags,
# 3 chains, each with 34000 iterations (first 14000 discarded), n.thin = 10
# n.sims = 6000 iterations saved
# mu.vect sd.vect   2.5%    25%    50%    75%  97.5% Rhat n.eff
# C            1.44    0.96   0.12   0.69   1.26   2.00   3.70 1.00   720
# D            1.13    1.06   0.03   0.35   0.80   1.60   3.87 1.00  1400
# K_V          1.16    0.63   0.30   0.67   1.06   1.55   2.64 1.04    56
# Q           12.26    2.80   7.51  10.23  12.06  13.96  18.53 1.00   820
### Better though not oerfect given centered on true values... 

### Tried to increase the shape parameter to have a well-defined peak dgamma(10,1)
# Inference for Bugs model at "ssm.predprey1.txt", fit using jags,
# 3 chains, each with 34000 iterations (first 14000 discarded), n.thin = 10
# n.sims = 6000 iterations saved
# mu.vect sd.vect   2.5%    25%    50%    75%  97.5% Rhat n.eff
# C            7.70    2.20   4.02   6.10   7.51   9.10  12.49 1.00  6000
# D           11.45    3.28   5.93   9.11  11.13  13.50  18.62 1.00  2200
# K_V          1.45    0.62   0.52   1.01   1.35   1.80   2.92 1.07    34
# Q           12.28    2.75   7.55  10.35  12.05  14.03  18.23 1.00  3200
### Not fully satisfying either

### With dgamma(2,1)
# Inference for Bugs model at "ssm.predprey1.txt", fit using jags,
# 3 chains, each with 34000 iterations (first 14000 discarded), n.thin = 10
# n.sims = 6000 iterations saved
# mu.vect sd.vect   2.5%    25%    50%    75%  97.5% Rhat n.eff
# C           1.61    0.95   0.26   0.90   1.46   2.15   3.89 1.00  6000
# D           2.07    1.40   0.27   1.03   1.78   2.76   5.51 1.00  6000
# K_V         1.23    0.58   0.33   0.82   1.15   1.55   2.66 1.02   170
# Q          12.40    2.88   7.51  10.42  12.17  14.18  18.60 1.00  3200
# r_P         0.56    0.10   0.37   0.49   0.55   0.62   0.75 1.00  3200
# r_V         1.89    0.43   1.20   1.58   1.83   2.13   2.96 1.01   160
# sigma2_P    0.04    0.01   0.03   0.04   0.04   0.04   0.05 1.00  6000
# sigma2_V    0.06    0.01   0.04   0.05   0.06   0.06   0.07 1.00  6000
# deviance  -29.98    3.50 -34.84 -32.53 -30.61 -28.16 -21.49 1.00  2900
# 
