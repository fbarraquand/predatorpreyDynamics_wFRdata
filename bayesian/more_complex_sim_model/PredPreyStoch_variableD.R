### F. Barraquand 21/04/2015 - Code for analyzing noisy predator-prey system (incl. noise on the functional response) 
### 17/03/19 -- Modified to make it clear we consider a temporally variable D (half-saturation constant)

rm(list=ls())
graphics.off()
library("R2jags")      # Load R2jags package


### Parameters for simulation
n.years<-100  	# Number of years - 25 first, perhaps use 50 or 100 / worked very well with almost no process error on the real scale
N1<-1			# Initial pop size
P1<-0.1
K<-1			# threshold dd 
beta<-1			# density-dependence exponent
rmax_V<-2			# Max AVERAGE growth rate (thus not a true max...)
rmax_P<-0.5
sigma2.proc<-0.5 #1		# Process sigma on the log-scale


### FR and predator parameters
C<-2.5
D1<-2
Q<-10

## Temporal variability of D is modelled with a gamma distribution
## D~dgamma(sh,ra) ## shape, rate parameters
mD <-2 # mean or mode of D, which is the half-saturation constant
sD <-4 #reasonable value for SD
sh <- mD^2 / sD^2
ra <- mD/ sD^2
# http://doingbayesiandataanalysis.blogspot.se/2012/08/gamma-likelihood-parameterized-by-mode.html

### Simulation of data
set.seed(42) 
y<-N<-P<-D<-FR<-numeric(n.years)
N[1]<-N1
P[1]<-P1
D[1]<-D1
FR[1]<-C*N[1]/(D[1]+N[1])

rV<-rnorm(n.years-1,rmax_V,sqrt(sigma2.proc))
rP<-rnorm(n.years-1,rmax_P,sqrt(sigma2.proc))

for (t in 1:(n.years-1)){

  N[t+1]<-N[t]*(exp(rV[t])/(1+(N[t]/K)^beta))*exp(-FR[t]*P[t]/N[t])
  P[t+1]<-P[t]*exp(rP[t])/(1+P[t]*Q/N[t])
  D[t+1] <- rgamma(1,sh,ra) 
  FR[t+1]<-C*N[t+1]/(D[t]+N[t+1]) 
}
## Plotting time series of abundances and FR
par(mfrow=c(2,2))
plot(1:n.years,N,type="o")
plot(1:n.years,P,type="o")
curve(dgamma(x,sh,ra),from=0, to=10)
plot(N,FR)


#seq(0,1,0.01)
# Bundle data
jags.data <- list(T=n.years,logN=log(N),logP=log(P),FR=FR)
## Actually both C and D are constant... 
sink("predprey_constantC.txt")
cat("
    model {
    
    # Priors and constraints
    logN[1] ~ dnorm(0,0.01) # Prior on initial pop size on the log scale
    logP[1] ~ dnorm(0,0.01) # Prior on initial pop size on the log scale
    FR[1] ~ dnorm(C*exp(logN[1])/(D+exp(logN[1])),100)
    
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
    C~dgamma(.01,.01) # uninformative priors OK 
    D~dgamma(0.01,0.01)
    
    # Likelihood
    # state process
    
    for (t in 1:(T-1)){        
    
    FRUpdate[t] <- C*N[t]/(D+N[t]) #functional response equation, including noise
    FR[t+1] ~  dnorm(FRUpdate[t],tau_FR) # put Gaussian noise, probably not the best but works
    
    logNupdate[t] <- logN[t] + r_V -log(1+N[t]/K_V) -FR[t+1]*exp(logP[t])/N[t]
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
  list(sigma_V=runif(1,0.1,2), sigma_P=runif(1,0.1,2), r_V=runif(1,0.1,2),r_P=runif(1,0.1,2), K_V=runif(1,0.2,8), Q=runif(1,0,5),tau_FR=runif(1,1,10),C=runif(1,10,100),D=runif(1,0.01,30))}

# Parameters monitored
#parameters<-c("r_V","K_V","r_P","Q","sigma2_V","sigma2_P","a","b","C","D","logN","logP","FR")
parameters<-c("r_V","K_V","r_P","Q","sigma2_V","sigma2_P","tau_FR","C","D")


# MCMC settings
nc <- 3 #number of chains
nb <- 14000 # “burn in”
#ni <- 14000# “number of iterations” # that's for a symmetric distrib...
ni<-34000
nt <- 10 # “thinning”

# run model
out <- jags(jags.data, inits, parameters, "predprey_constantC.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, working.directory = getwd())
print(out, dig = 2)

# Inference for Bugs model at "predprey_constantC.txt", fit using jags,
# 3 chains, each with 34000 iterations (first 14000 discarded), n.thin = 10
# n.sims = 6000 iterations saved
# mu.vect sd.vect   2.5%    25%    50%    75%  97.5% Rhat n.eff
# C           2.31    0.10   2.13   2.24   2.31   2.37   2.51    1  6000
# D           0.52    0.13   0.30   0.43   0.51   0.60   0.79    1  6000
# K_V         1.59    0.89   0.34   0.96   1.43   2.04   3.78    1  2000
# Q          10.64    3.06   5.65   8.48  10.33  12.46  17.53    1  6000
# r_P         0.44    0.11   0.24   0.37   0.44   0.51   0.67    1  6000
# r_V         1.81    0.47   1.07   1.48   1.74   2.06   3.00    1  1800
# sigma2_P    0.41    0.06   0.31   0.37   0.41   0.45   0.55    1  6000
# sigma2_V    0.59    0.09   0.44   0.53   0.58   0.64   0.78    1  6000
# tau_FR      3.10    0.47   2.26   2.77   3.07   3.39   4.09    1  6000
# deviance  607.17    4.45 600.55 603.93 606.50 609.71 617.54    1  3200
# 

### Fairly good!!
fr_fit<-nls(FR~CE*N/(DE+N),start=list(CE=1,DE=1))
CE<-coef(fr_fit)[1]
DE<-coef(fr_fit)[2]
plot(N,FR,ylim=c(0,3))
curve(CE*x/(DE+x),col="blue",add=TRUE)
CE #2.326468 
DE #0.3638511 

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
    
    C~dgamma(.01,.01) # 
    D~dgamma(.01,.01)
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
    
    ### Need to add a functional response model in there - otherwise there's just not the same amount of data
    
    }
    
    }
    ",fill=TRUE)
sink()


# Initial values
inits <- function () {
  list(sigma_V=runif(1,0.1,2), sigma_P=runif(1,0.1,2), r_V=runif(1,0.1,2),r_P=runif(1,0.1,2), K_V=runif(1,0.2,8), Q=runif(1,0,5),C=runif(1,10,100),D=runif(1,0.01,0.1))}

# Parameters monitored
#parameters<-c("r_V","K_V","r_P","Q","sigma2_V","sigma2_P","a","b","C","D","logN","logP","FR")
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


#### Fit the exact simulated model 
sink("predprey_variableD.txt")
cat("
    model {
    
    # Priors and constraints
    logN[1] ~ dnorm(0,0.01) # Prior on initial pop size on the log scale
    logP[1] ~ dnorm(0,0.01) # Prior on initial pop size on the log scale
    FR[1] ~ dnorm(C*exp(logN[1])/(mD+exp(logN[1])),100)
    
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
    C~dgamma(.01,.01) # uninformative priors OK 
    mD~dgamma(0.1,0.1)
    sD~dgamma(0.1,0.1)
    sh <- pow(mD,2) / pow(sD,2)
    ra <-     mD    / pow(sD,2)
    
    # Likelihood
    # state process
    
    for (t in 1:(T-1)){        
    
    D[t] ~ dgamma(sh,ra)
    FRUpdate[t] <- C*N[t]/(D[t]+N[t]) #functional response equation, including noise
    FR[t+1] ~  dnorm(FRUpdate[t],tau_FR) # technical trick
    
    logNupdate[t] <- logN[t] + r_V -log(1+N[t]/K_V) -FR[t+1]*exp(logP[t])/N[t]
    logN[t+1] ~ dnorm(logNupdate[t],tau_V)
    N[t]<-exp(logN[t])
    
    logP[t+1]~ dnorm(logPupdate[t],tau_P)
    logPupdate[t] <- logP[t] + r_P - log(1+exp(logP[t])*Q/exp(logN[t]) )  
    
    }
    
    }
    ",fill=TRUE)
sink()


# Initial values
inits <- function () {
  list(sigma_V=runif(1,0.1,2), sigma_P=runif(1,0.1,2), r_V=runif(1,0.1,2),r_P=runif(1,0.1,2), K_V=runif(1,0.2,8), Q=runif(1,0,5),tau_FR=runif(1,1,10),C=runif(1,10,100),mD=runif(1,0.5,3),sD=runif(1,0.5,6))}

# Parameters monitored
parameters<-c("r_V","K_V","r_P","Q","sigma2_V","sigma2_P","mD","sD","C","tau_FR")#"logN","logP","FR"


# MCMC settings
nc <- 3 #number of chains
nb <- 14000 # “burn in”
#ni <- 14000# “number of iterations” # that's for a symmetric distrib...
ni<-34000
nt <- 10 # “thinning”

# run model
out3 <- jags(jags.data, inits, parameters, "predprey_variableD.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, working.directory = getwd())
print(out3, dig = 2)
# 
# Inference for Bugs model at "predprey_variableD.txt", fit using jags,
# 3 chains, each with 34000 iterations (first 14000 discarded), n.thin = 10
# n.sims = 6000 iterations saved
# mu.vect sd.vect   2.5%     25%     50%     75%   97.5% Rhat n.eff
# C           2.49    0.00   2.48    2.49    2.49    2.49    2.50 1.00   710
# K_V         1.60    0.87   0.40    0.97    1.44    2.04    3.70 1.01   250
# Q          10.60    3.00   5.62    8.48   10.34   12.31   17.39 1.00  6000
# mD          1.28    0.17   0.98    1.16    1.27    1.39    1.64 1.00   630
# r_P         0.44    0.11   0.24    0.37    0.44    0.51    0.66 1.00  6000
# r_V         1.79    0.44   1.09    1.48    1.74    2.06    2.84 1.01   240
# sD          3.37    0.48   2.56    3.02    3.32    3.66    4.45 1.00  6000
# sigma2_P    0.41    0.06   0.31    0.37    0.41    0.45    0.55 1.00  6000
# sigma2_V    0.59    0.09   0.45    0.53    0.58    0.64    0.79 1.00  6000
# tau_FR   1372.32  337.77 796.43 1136.43 1340.03 1579.24 2108.21 1.01   430
# deviance  -27.21   26.50 -76.55  -45.72  -28.19   -9.60   26.95 1.01   400
mD #2
sD #4
### In that case, it is a better option to let tau_FR be estimated to a large value,
### rather than force it to be small (which impedes convergence)


