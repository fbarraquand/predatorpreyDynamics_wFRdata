### F. Barraquand 21/04/2015 - Code for analyzing noisy predator-prey system (incl. noise on the functional response) 
### Edited 25/05/2015 Add of more diagnostic plots. Problems with the functional response model, I should probably
# analyze it separarely to see if all components can be identified - this looks quite unclear. 
# There was by the way an error in the simulation, FR[t+1]<-C*(1-B)*N[t]/(D+N[t]) was FR[t+1]<-C*(1-B)*N[1]/(D+N[1]) 
# So the FR was rightfully more constant!!

### Predator prey-only, fit of functional response (28/04/2015) but direct estimation (same model simulated as fitted). 
# small previous problem with Beverton-Holt formulation. Now Getz formulation, better. 

### For the functional response, the beta is perhaps good for simulations, but I feel like we need something else
# to fit the model. Perhaps a logNormal. -> A normal distrib is used here, note it could be truncated. 

rm(list=ls())
graphics.off()

library("R2jags")      # Load R2jags package
library("modeest")

### Parameters for simulation of Hassell model

n.years<-100  	# Number of years - 25 first, perhaps use 50 or 100 / worked very well with almost no process error on the real scale
N1<-1			# Initial pop size
P1<-0.1
K<-1			# threshold dd 
beta<-1			# density-dependence exponent
rmax_V<-2			# Max AVERAGE growth rate (thus not a true max...)
rmax_P<-0.5
sigma2.proc<-0.5 #1		# Process sigma on the log-scale


### FR and predator parameters
a<-2 #0.5 10 
b<-5   #3 70
C<-0.5
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
  FR[t+1]<-C*(1-B)*N[t]/(D+N[t]) 
  N[t+1]<-N[t]*(exp(rV[t])/(1+(N[t]/K)^beta))*exp(-FR[t]*P[t]/N[t])
  P[t+1]<-P[t]*exp(rP[t])/(1+P[t]*Q/N[t])
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

sink("ssm.predprey3.txt")
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
    C~dgamma(.01,.01) # uninformative priors OK for that one

    # Setting up informative priors for D, otherwise the estimation does not work. 
    D~dgamma(sh,ra)
    sh <- pow(mD,2) / pow(sd,2)
    ra <-     mD    / pow(sd,2)
    #sh <- 1 + mD * ra #shape param as a function of mode
    #ra <- ( mD + sqrt( mD^2 + 4*sd^2 ) ) / ( 2 * sd^2 )
    mD <-1.5# mean or mode of D, which is the half-saturation constant
    sd <-2 #reasonable values
    # http://doingbayesiandataanalysis.blogspot.se/2012/08/gamma-likelihood-parameterized-by-mode.html
    # check model Leslie to see how she specified priors...

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
out <- jags(jags.data, inits, parameters, "ssm.predprey3.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, working.directory = getwd())
print(out, dig = 2)

# Good title for a future paper (if I compute more quantities of Trophic Strength from this...) would be
# An integrated assessment of trophic interaction strength. 

# Output summary statistics
jags.sum<-out$BUGSoutput$summary
write.table(x=jags.sum,file="JAGSsummary_PredPreyStoch_analysisTer.txt")
# MCMC Output
pdf("Output_MCMC__PredPreyStoch_analysisTer.pdf")
out.mcmc<-as.mcmc(out)
plot(out.mcmc) 
dev.off()

### Still a problem of no convergence of D even when setting (a,b) to large values... Ask Leslie?

#Isolated model for the functional response only

sink("func.resp1.txt")
cat("
    model {
    
    # Priors and constraints
    FR[1] ~ dnorm(C*exp(logN[1])/(D+exp(logN[1])),100)
    
    #Priors predation parameters 
    tau_FR ~ dgamma(.01,.01)
    C~dgamma(.01,.01) # uninformative priors OK for that one
    
    # Setting up informative priors for D, otherwise the estimation does not work. 
    D~dgamma(sh,ra)
    sh <- pow(mD,2) / pow(sd,2)
    ra <-     mD    / pow(sd,2)
    mD <-1.5# mean or mode of D, which is the half-saturation constant
    sd <-2 #reasonable values
    # http://doingbayesiandataanalysis.blogspot.se/2012/08/gamma-likelihood-parameterized-by-mode.html
    # check model Leslie to see how she specified priors...
    
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
plot(N,FR,ylim=c(0,0.3))
lines(N,CE*N/(DE+N))
lines(N,CEb*N/(DEb+N))
CE
DE
# so actually the Bayesian estimation is more correct... difficult to estimate a half-saturation constant in that case anyway...
# most likely, in the gyrfalcon case, even with an informative prior is will be difficult. 

### Added 30/06/2016
### Now we should add the predator-prey model without the FR


### Now try to fit a model without the FR data. 
sink("ssm.predprey_without_sepFR.txt")
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
out2 <- jags(jags.data, inits, parameters, "ssm.predprey_without_sepFR.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, working.directory = getwd())
print(out2, dig = 2)
# http://jeromyanglim.tumblr.com/post/37362047458/how-to-get-dic-in-jags
print(out,dig=2) # to compare deviance, DIC - check also parameter values.
# 
plot(as.mcmc(out2)) 
plot(as.mcmc(out)) 
