### F. Barraquand 21/04/2015 - Code for analyzing noisy predator-prey system (incl. noise on the functional response) 
### Predator prey-only, but no fit of functional response (28/04/2015)

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
sigma2.proc<-0.05		# Process sigma on the log-scale

C<-2.5
D<-1
Q<-10

### Simulation of data
#set.seed(41) 
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
plot(1:n.years,N,type="b",ylim=c(0,max(N)), col="blue")
par(new = TRUE)
plot(1:n.years,P, axes=F, col="red",xlab=NA, ylab=NA,type="b")
axis(side = 4)
#mtext(side = 4, line = 3, 'P')

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


# Initial values
inits <- function () {
  list(sigma_V=runif(1,0.1,2), sigma_P=runif(1,0.1,2), r_V=runif(1,0.1,2),r_P=runif(1,0.1,2), K_V=runif(1,0.2,10), Q=runif(1,0,5),C=runif(1,10,100),D=runif(1,0.01,0.1))}


# Parameters monitored
parameters<-c("r_V","K_V","r_P","Q","sigma2_V","sigma2_P","C","D","logN","logP")

# MCMC settings
nc <- 3 #number of chains
nb <- 14000 # “burn in”
#ni <- 14000# “number of iterations” # that's for a symmetric distrib...
ni<-34000
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

## Serious identifiability problems on C and D here. 

### Old results 

# Inference for Bugs model at "ssm.predprey1.txt", fit using jags,
# 3 chains, each with 34000 iterations (first 14000 discarded), n.thin = 10
# n.sims = 6000 iterations saved
# mu.vect sd.vect   2.5%    25%    50%    75%  97.5% Rhat n.eff
# C            2.18    4.11   0.00   0.00   0.14   2.34  15.08 1.47     9
# D            1.24    1.95   0.00   0.00   0.00   2.27   6.47 1.19    15
# K_V          0.99    1.71   0.20   0.21   0.25   0.72   6.75 2.01     5
# Q           10.17    2.29   6.22   8.57   9.94  11.57  15.15 1.00  6000

### New results

# Inference for Bugs model at "ssm.predprey1.txt", fit using jags,
# 3 chains, each with 34000 iterations (first 14000 discarded), n.thin = 10
# n.sims = 6000 iterations saved
# mu.vect sd.vect   2.5%    25%    50%    75%  97.5% Rhat n.eff
# C            0.00    0.01   0.00   0.00   0.00   0.00   0.01 1.31    10
# D            0.01    0.08   0.00   0.00   0.00   0.00   0.05 1.56     8
# K_V          1.50    0.70   0.52   0.99   1.39   1.88   3.19 1.01   190
# Q           12.27    2.81   7.53  10.33  12.00  14.01  18.38 1.00  4500

### Estimation problems here 
