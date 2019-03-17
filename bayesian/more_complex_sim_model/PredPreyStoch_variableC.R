### F. Barraquand 21/04/2015 - Code for analyzing noisy predator-prey system (incl. noise on the functional response) 
### 17/03/19 -- Modified to make it clear we consider a variable C quantity

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
### Sigma2 = 0.05 should work as well. 


### FR and predator parameters
a<-2 #0.5 10 
b<-5   #3 70
C<-2.5
D<-1
Q<-10

### Simulation of data
set.seed(42) 
y<-N<-P<-FR<-numeric(n.years)
N[1]<-N1
P[1]<-P1
FR[1]<-C*N[1]/(D+N[1])

rV<-rnorm(n.years-1,rmax_V,sqrt(sigma2.proc))
rP<-rnorm(n.years-1,rmax_P,sqrt(sigma2.proc))

for (t in 1:(n.years-1)){
  B=rbeta(1,a,b) ### Here I simulate the exact model that we'll fit later
  N[t+1]<-N[t]*(exp(rV[t])/(1+(N[t]/K)^beta))*exp(-FR[t]*P[t]/N[t])
  P[t+1]<-P[t]*exp(rP[t])/(1+P[t]*Q/N[t])
  FR[t+1]<-C*(1-B)*N[t+1]/(D+N[t+1]) 
}
## Plotting time series of abundances and FR
par(mfrow=c(2,2))
plot(1:n.years,N,type="b")
plot(1:n.years,P,type="b")
curve(dbeta(x,a,b),from=0, to=1)
plot(N,FR)

#seq(0,1,0.01)
# Bundle data
jags.data <- list(T=n.years,logN=log(N),logP=log(P),FR=FR)

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


# Output summary statistics
jags.sum<-out$BUGSoutput$summary
write.table(x=jags.sum,file="JAGSsummary_PredPreyStoch_constantC.txt")
# MCMC Output
pdf("Output_MCMC__PredPreyStoch_constantC.pdf")
out.mcmc<-as.mcmc(out)
plot(out.mcmc) 
dev.off()


#Isolated model for the functional response only

sink("func.resp1.txt")
cat("
    model {
    
    # Priors and constraints
    FR[1] ~ dnorm(C*exp(logN[1])/(D+exp(logN[1])),100)
    
    #Priors predation parameters 
    tau_FR ~ dgamma(.01,.01)
    C~dgamma(.01,.01) # uninformative priors OK for that one
    
    # Setting up informative priors for D
    D~dgamma(sh,ra)
    sh <- pow(mD,2) / pow(sd,2)
    ra <-     mD    / pow(sd,2)
    mD <-1.5# mean or mode of D, which is the half-saturation constant
    sd <-2 #reasonable values
    # http://doingbayesiandataanalysis.blogspot.se/2012/08/gamma-likelihood-parameterized-by-mode.html
    
    # Likelihood
    # state process
    
    for (t in 1:(T-1)){        
    N[t]<-exp(logN[t])
    FRUpdate[t] <- C*N[t]/(D+N[t]) #functional response equation, including noise
    FR[t+1] ~  dnorm(FRUpdate[t],tau_FR) # put Gaussian noise, probably not the best but works
    
    }
    
    }
    ",fill=TRUE)
sink()


# Initial values
inits <- function () {  list(tau_FR=runif(1,1,10),C=runif(1,10,100),D=runif(1,0.01,30))}

# Parameters monitored
#parameters<-c("r_V","K_V","r_P","Q","sigma2_V","sigma2_P","a","b","C","D","logN","logP","FR")
parameters<-c("tau_FR","C","D")


# MCMC settings
nc <- 3 #number of chains
nb <- 14000 # “burn in”
#ni <- 14000# “number of iterations” # that's for a symmetric distrib...
ni<-34000
nt <- 10 # “thinning”

# run model
out.f1 <- jags(jags.data, inits, parameters, "func.resp1.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, working.directory = getwd())
print(out.f1, dig = 2)
CEb<-out.f1$BUGSoutput$mean$C
DEb<-out.f1$BUGSoutput$mean$D
### So we have a problem with that model in particular.
# Try classic nonlinear fit

fr_fit<-nls(FR~CE*N/(DE+N),start=list(CE=1,DE=1))
CE<-coef(fr_fit)[1]
DE<-coef(fr_fit)[2]
plot(N,FR,ylim=c(0,3))
curve(CE*x/(DE+x),col="blue",add=TRUE)
curve(CEb*x/(DEb+x),col="red",add=TRUE)
CE
DE

### Added 30/06/2016
### Now we should add the predator-prey model without the FR


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

print(out,dig=2)
# 
plot(as.mcmc(out2)) 
plot(as.mcmc(out)) 


#### Fit the exact simulated model 
sink("predprey_variableC.txt")
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
    a~dgamma(0.1,0.1)
    b~dgamma(0.1,0.1)
    
    # Likelihood
    # state process
    
    for (t in 1:(T-1)){        
    
    B[t] ~ dbeta(a,b)
    FRUpdate[t] <- C*(1-B[t])*N[t]/(D+N[t]) #functional response equation, including noise
    FR[t+1] ~  dnorm(FRUpdate[t],tau_FR) # put Gaussian noise, probably not the best but works
    
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
  list(sigma_V=runif(1,0.1,2), sigma_P=runif(1,0.1,2), r_V=runif(1,0.1,2),r_P=runif(1,0.1,2), K_V=runif(1,0.2,8), Q=runif(1,0,5),tau_FR=runif(1,1,10),C=runif(1,10,100),D=runif(1,0.01,30),a=runif(1,0.5,3),b=runif(1,0.5,3))}

# Parameters monitored
parameters<-c("r_V","K_V","r_P","Q","sigma2_V","sigma2_P","a","b","C","D")#"logN","logP","FR"


# MCMC settings
nc <- 3 #number of chains
nb <- 14000 # “burn in”
#ni <- 14000# “number of iterations” # that's for a symmetric distrib...
ni<-34000
nt <- 10 # “thinning”

# run model
out3 <- jags(jags.data, inits, parameters, "predprey_variableC.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, working.directory = getwd())
print(out3, dig = 2)
# 
# Inference for Bugs model at "predprey_variableC.txt", fit using jags,
# 3 chains, each with 34000 iterations (first 14000 discarded), n.thin = 10
# n.sims = 6000 iterations saved
# mu.vect sd.vect   2.5%    25%    50%    75%  97.5% Rhat n.eff
# C           2.29    0.18   1.92   2.18   2.30   2.40   2.62 1.01   230
# D           0.96    0.16   0.65   0.86   0.96   1.07   1.27 1.00  1200
# K_V         2.13    1.06   0.62   1.38   1.93   2.68   4.74 1.00  1400
# Q          10.75    3.07   5.69   8.61  10.45  12.53  17.80 1.00  6000
# a           1.29    0.63   0.24   0.87   1.23   1.62   2.79 1.03   150
# b           4.15    1.22   2.04   3.34   4.04   4.81   6.95 1.01   210
# r_P         0.45    0.11   0.24   0.37   0.44   0.51   0.66 1.00  6000
# r_V         1.54    0.39   0.91   1.27   1.50   1.75   2.42 1.00  1300
# sigma2_P    0.41    0.06   0.31   0.37   0.41   0.45   0.55 1.00  2300
# sigma2_V    0.56    0.08   0.42   0.50   0.55   0.61   0.74 1.00  6000
# deviance  309.08   71.21 166.29 258.87 312.86 363.31 433.03 1.00   690

### Fairly good!
