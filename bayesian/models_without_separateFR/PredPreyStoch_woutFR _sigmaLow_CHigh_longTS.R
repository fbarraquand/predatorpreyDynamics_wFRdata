### F. Barraquand 21/04/2015 - Code for analyzing noisy predator-prey system (incl. noise on the functional response) 
### Predator prey-only, but no fit of functional response (28/04/2015)

# Looks like the FR parameters are the most difficult to estimate indirectly, which is interesting for a 
# joint model fit. 

# small previous problem with Beverton-Holt formulation. Now Getz formulation, better. 

rm(list=ls())
graphics.off()

library("R2jags")      # Load R2jags package


### Parameters for simulation of Hassell model

n.years<-1000  	# Number of years
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
png(file = "TracePlot_T=1000.png")
traplot(out,c("r_V","K_V","r_P","Q","sigma2_V","sigma2_P","C","D"))
dev.off()
### plot densities
denplot(out,c("r_V","K_V","r_P","Q","sigma2_V","sigma2_P","C","D"))


## Still identifiability problems on C and D here, though less severe on C
out$BUGSoutput$mean$C
out$BUGSoutput$mean$D
# > out$BUGSoutput$mean$C
# [1] 2.033202
# > out$BUGSoutput$mean$D
# [1] 0.01934797

## check chain by chain
seq1=seq(1,by=3,6000)
seq2=seq(2,by=3,6000)
seq3=seq(3,by=3,6000)

plot(out$BUGSoutput$sims.list$D)

mean(out$BUGSoutput$sims.list$D[seq1])
plot(out$BUGSoutput$sims.list$D[seq1])
mean(out$BUGSoutput$sims.list$D[seq2])
lines(out$BUGSoutput$sims.list$D[seq2],col="red")
mean(out$BUGSoutput$sims.list$D[seq3])
lines(out$BUGSoutput$sims.list$D[seq3],col="blue")
### Not same results as traplot, weird

### Check it works for C
mean(out$BUGSoutput$sims.list$C[seq1])
plot(out$BUGSoutput$sims.list$C[seq1],ylim=c(0,7))
mean(out$BUGSoutput$sims.list$C[seq2])
lines(out$BUGSoutput$sims.list$C[seq2],col="red")
mean(out$BUGSoutput$sims.list$C[seq3])
lines(out$BUGSoutput$sims.list$C[seq3],col="blue")
### Not same results as traplot, weird

### different order? 
seq1=seq(1,by=1,2000)
seq2=seq(2001,by=1,4000)
seq3=seq(4001,by=1,6000)

mean(out$BUGSoutput$sims.list$C[seq1])
plot(out$BUGSoutput$sims.list$C[seq1],ylim=c(0,6))
mean(out$BUGSoutput$sims.list$C[seq2])
lines(out$BUGSoutput$sims.list$C[seq2],col="red")
mean(out$BUGSoutput$sims.list$C[seq3])
lines(out$BUGSoutput$sims.list$C[seq3],col="blue")
### still does not work

# check again with another function
traceplot(out) #does the same thing. 

### Well this is very annoying!!!
# answer here https://sourceforge.net/p/mcmc-jags/discussion/610037/thread/cc61b820/
  
# And using sims.array is visibly the only thing correct way to get the trace...
plot(out$BUGSoutput$sims.array[,1,'C'],type='l',col="black")
lines(out$BUGSoutput$sims.array[,2,'C'],type='l',col="red")
lines(out$BUGSoutput$sims.array[,3,'C'],type='l',col="blue")

mean(out$BUGSoutput$sims.array[,1,'C'])
mean(out$BUGSoutput$sims.array[,2,'C'])
mean(out$BUGSoutput$sims.array[,3,'C'])

plot(out$BUGSoutput$sims.array[,1,'D'],type='l',col="black")
lines(out$BUGSoutput$sims.array[,2,'D'],type='l',col="red")
lines(out$BUGSoutput$sims.array[,3,'D'],type='l',col="blue")

mean(out$BUGSoutput$sims.array[,1,'D'])
mean(out$BUGSoutput$sims.array[,2,'D'])
mean(out$BUGSoutput$sims.array[,3,'D'])
# all chains are crap for D. 

### Does that change the correlation between parameters? use r_V and K_V as we know they are correlated
cor(as.vector(out$BUGSoutput$sims.array[,,'r_V']),as.vector(out$BUGSoutput$sims.array[,,'K_V']))
cor(as.vector(out$BUGSoutput$sims.list$r_V),as.vector(out$BUGSoutput$sims.list$K_V))
### Thankfully sims.list preserves the correlation between parameters!
