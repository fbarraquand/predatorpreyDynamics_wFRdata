### F. Barraquand 21/04/2015 - Code for analyzing noisy predator-prey system (incl. noise on the functional response) 
### 17/03/19 -- Modified to make it clear we consider a variable C quantity

rm(list=ls())
graphics.off()

library("R2jags")      # Load R2jags package

### Parameters for simulation
n.years<-100  	# Number of years - 25 first, perhaps use 50 or 100 / worked very well with almost no process error on the real scale
N1<-1			# Initial pop size
P1<-0.1
gamma<-1			# threshold dd coefficient
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
  N[t+1]<-N[t]*(exp(rV[t])/(1+gamma*N[t]))*exp(-FR[t]*P[t]/N[t])
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

sink("predprey_gaussianFR.txt")
cat("
    model {
    
    # Priors and constraints
    logN[1] ~ dnorm(0,0.01) # Prior on initial pop size on the log scale
    logP[1] ~ dnorm(0,0.01) # Prior on initial pop size on the log scale
    FR[1] ~ dnorm(C*exp(logN[1])/(D+exp(logN[1])),100)

    # Priors prey population dynamics
    r_V ~ dnorm(1,0.001) # below the truth, rather flat prior
    gamma ~ dunif(0.1,5)
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
    
    logNupdate[t] <- logN[t] + r_V -log(1+gamma*N[t]) -FR[t+1]*exp(logP[t])/N[t]
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
  list(sigma_V=runif(1,0.1,2), sigma_P=runif(1,0.1,2), r_V=runif(1,0.1,2),r_P=runif(1,0.1,2), gamma=runif(1,0.3,4), Q=runif(1,0,5),tau_FR=runif(1,1,10),C=runif(1,10,100),D=runif(1,0.01,30))}

# Parameters monitored
#parameters<-c("r_V","gamma","r_P","Q","sigma2_V","sigma2_P","a","b","C","D","logN","logP","FR")
parameters<-c("r_V","gamma","r_P","Q","sigma2_V","sigma2_P","tau_FR","C","D")


# MCMC settings
nc <- 3 #number of chains
nb <- 24000 # “burn in”
#ni <- 14000# “number of iterations” # that's for a symmetric distrib...
ni<-44000
nt <- 10 # “thinning”

# run model
out <- jags(jags.data, inits, parameters, "predprey_gaussianFR.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, working.directory = getwd())
print(out, dig = 2)

# MCMC Output
#pdf("Output_MCMC__PredPreyStoch_gaussianFR.pdf")
#out.mcmc<-as.mcmc(out)
#plot(out.mcmc) 
#dev.off()

### We don't know what are the "true" values of the functional response in a Gaussian setting
###, so we'll take those of a best fit least square
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
parameters<-c("tau_FR","C","D")


# MCMC settings
nc <- 3 #number of chains
nb <- 24000 # “burn in”
#ni <- 14000# “number of iterations” # that's for a symmetric distrib...
ni<-44000
nt <- 10 # “thinning”

# run model
out.f1 <- jags(jags.data, inits, parameters, "func.resp1.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, working.directory = getwd())
print(out.f1, dig = 2)
CEb<-out.f1$BUGSoutput$mean$C
DEb<-out.f1$BUGSoutput$mean$D
taub<-out.f1$BUGSoutput$mean$tau_FR

pdf(file="GaussianFuncResp_fittedonVariableC.pdf",width=8,height=5)
fr_fit<-nls(FR~CE*N/(DE+N),start=list(CE=1,DE=1))
CE<-coef(fr_fit)[1]
DE<-coef(fr_fit)[2]
plot(N,FR,ylim=c(0,3))
curve(CE*x/(DE+x),col="blue",add=TRUE)
curve(CEb*x/(DEb+x),col="red",add=TRUE)
abline(h=C,col="black")
CE
DE
dev.off()

### Now back to the full model with Gaussian FR

# Output summary statistics
jags.sum<-out$BUGSoutput$summary
write.table(x=jags.sum,file="JAGSsummary_PredPreyStoch_gaussianFR.txt")

library(xtable)
xtable(jags.sum)

library(mcmcplots)
?caterplot
par(mfrow=c(1,1))
parameters<-c("r_V","gamma","r_P","Q","sigma2_V","sigma2_P","tau_FR","C","D")
caterplot(out,parms = parameters)  ## more fuzzy stuff denstrip = TRUE
true_values=c(rmax_V,gamma,rmax_P,Q,sigma2.proc,sigma2.proc,taub,CEb,DEb)
true_values_neworder = c(DEb,sigma2.proc,rmax_P,sigma2.proc,CEb,gamma,rmax_V,taub,Q)
caterpoints(rev(true_values_neworder),pch="x",col="red")


### Now try to fit a model without the FR data. 
sink("predprey_without_sepFR.txt")
cat("
    model {
    
    # Priors and constraints
    logN[1] ~ dnorm(0,0.01) # Prior on initial pop size on the log scale
    logP[1] ~ dnorm(0,0.01) # Prior on initial pop size on the log scale
    
    # Priors prey population dynamics
    r_V ~ dnorm(1,0.001) # below the truth, rather flat prior
    gamma ~ dunif(0.1,5)
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
    
    
    logNupdate[t] <- logN[t] + r_V -log(1+gamma*N[t]) -C*exp(logP[t])/(D+N[t])
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
  list(sigma_V=runif(1,0.1,2), sigma_P=runif(1,0.1,2), r_V=runif(1,0.1,2),r_P=runif(1,0.1,2), gamma=runif(1,0.3,4), Q=runif(1,0,5),C=runif(1,10,100),D=runif(1,0.01,0.1))}

# Parameters monitored
parameters<-c("r_V","gamma","r_P","Q","sigma2_V","sigma2_P","C","D")


# MCMC settings
nc <- 3 #number of chains
nb <- 24000 # “burn in”
#ni <- 14000# “number of iterations” # that's for a symmetric distrib...
ni<-44000
nt <- 10 # “thinning”

# run model
out2 <- jags(jags.data, inits, parameters, "predprey_without_sepFR.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, working.directory = getwd())
print(out2, dig = 2)

print(out,dig=2)
# 
plot(as.mcmc(out2)) 
plot(as.mcmc(out)) 


# Output summary statistics
jags.sum<-out2$BUGSoutput$summary
xtable(jags.sum)
write.table(x=jags.sum,file="JAGSsummary_PredPreyStoch_constantFR.txt")

par(mfrow=c(1,1))
parameters<-c("r_V","gamma","r_P","Q","sigma2_V","sigma2_P","C","D")
caterplot(out2,parms = parameters)  ## more fuzzy stuff denstrip = TRUE
true_values_neworder = c(DEb,sigma2.proc,rmax_P,sigma2.proc,CEb,gamma,rmax_V,Q)
caterpoints(rev(true_values_neworder),pch="x",col="red")


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
    gamma ~ dunif(0.1,5) 
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
    sigma_FR ~ dunif(0,0.1) ## very very small
    tau_FR <- pow(sigma_FR,-2)
    # tau_FR ~ dgamma(0.001,0.001) ## I didn't manage to find good results with either flat or uninformative
    C~dgamma(.01,.01) # uninformative priors OK 
    D~dgamma(0.01,0.01)
    a~dgamma(0.1,0.1)
    b~dgamma(0.1,0.1)
    
    # Likelihood
    # state process
    
    for (t in 1:(T-1)){        
    
    B[t] ~ dbeta(a,b)
    FRUpdate[t] <- C*(1-B[t])*N[t]/(D+N[t]) #functional response equation, including noise
    FR[t+1] ~  dnorm(FRUpdate[t],tau_FR) # added layer of Gaussian noise, not the best but works
    #FR[t+1] ~  dnorm(FRUpdate[t],100000) # does not work -- modify the prior on tau_FR instead?
    
    logNupdate[t] <- logN[t] + r_V -log(1+N[t]*gamma) -FR[t+1]*exp(logP[t])/N[t]
    logN[t+1] ~ dnorm(logNupdate[t],tau_V)
    N[t]<-exp(logN[t])
    
    logP[t+1]~ dnorm(logPupdate[t],tau_P)
    logPupdate[t] <- logP[t] + r_P - log(1+exp(logP[t])*Q/exp(logN[t]) )  
    
    }
    
    }
    ",fill=TRUE)
sink()


# Parameters monitored
parameters<-c("r_V","gamma","r_P","Q","sigma2_V","sigma2_P","a","b","C","D","sigma_FR")#"logN","logP","FR"

# MCMC settings
nc <- 3 #number of chains
nb <- 24000 # “burn in”
ni<-44000
nt <- 10 # “thinning”

# Initial values
inits <- function () {
  list(sigma_V=runif(1,0.1,2), sigma_P=runif(1,0.1,2), r_V=runif(1,0.1,2),r_P=runif(1,0.1,2), gamma=runif(1,0.5,4), Q=runif(1,0,5),sigma_FR=runif(1,0,0.00005),C=runif(1,10,100),D=runif(1,0.01,30),a=runif(1,0.5,3),b=runif(1,0.5,3))}

# run model
out3 <- jags(jags.data, inits, parameters, "predprey_variableC.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, working.directory = getwd())
print(out3, dig = 2)

plot(as.mcmc(out3))

jags.sum<-out3$BUGSoutput$summary
xtable(jags.sum)
