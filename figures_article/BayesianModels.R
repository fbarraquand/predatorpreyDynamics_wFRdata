
#Create first predator-prey model with the stochastic functional response. 

sink("predprey.txt")
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
    tau_FR ~ dgamma(.01,.01)
    C~dgamma(.01,.01) # uninformative priors OK for that one
    D~dgamma(.01,.01)
    # check model Leslie to see how she specified priors...
    
    # Likelihood
    # state process
    
    for (t in 1:(T-1)){        
    
    FRUpdate[t] <- C*N[t]/(D+N[t]) #functional response equation, including noise
    FR[t] ~  dnorm(FRUpdate[t],tau_FR) #small trick to use KR data
    
    logNupdate[t] <- logN[t] + r_V -log(1+gamma*N[t]) -FR[t]*exp(logP[t])/N[t]
    logN[t+1] ~ dnorm(logNupdate[t],tau_V)
    N[t]<-exp(logN[t])
    
    logP[t+1]~ dnorm(logPupdate[t],tau_P)
    logPupdate[t] <- logP[t] + r_P - log(1+exp(logP[t])*Q/exp(logN[t]) )  
    
    }
    
    }
    ",fill=TRUE)
sink()

### Now model without the FR data. 
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

# Parameters monitored
parameters<-c("r_V","gamma","r_P","Q","sigma2_V","sigma2_P","C","D") 


# Initial values
inits <- function () {
  list(sigma_V=runif(1,0.1,2), sigma_P=runif(1,0.1,2), r_V=runif(1,0.1,2),r_P=runif(1,0.1,2), gamma=runif(1,0.3,4), Q=runif(1,0,5),tau_FR=runif(1,1,10),C=runif(1,10,100),D=runif(1,0.01,0.1))}


# MCMC settings
nc <- 3 #number of chains
nb <- 24000 # “burn in”
#ni <- 14000# “number of iterations” # that's for a symmetric distrib...
ni<-44000
nt <- 10 # “thinning”
