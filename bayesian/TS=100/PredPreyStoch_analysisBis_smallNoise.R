### F. Barraquand 21/04/2015 - Code for analyzing noisy predator-prey system (incl. noise on the functional response)
### Updated 01/06/2017. 
### Case with "small noise" on the functional response
### FB 03/01/2019 Added plots for the paired posterior distributions
### FB 14/01/2019 Corrected FR to nondelayed version

### Log of previous edits to other versions ######################################################################################

### Edited 25/05/2015 Add of more diagnostic plots. Problems with the functional response model, I should probably
# analyze it separarely to see if all components can be identified - this looks quite unclear. 
###### Edit again 25/05. actually no pb with the FR model, there's just that with the amount of noise, the Beta model for the FR
# is probably not identifiable and the FR half-sat constant is estimated at zero, which means a constant amount. 
# There was by the way initially an error in the simulation, FR[t+1]<-C*(1-B)*N[t]/(D+N[t]) was FR[t+1]<-C*(1-B)*N[1]/(D+N[1]) 
# So the FR was rightfully more constant!! Corrected now. 
# we can probably force a handling-time or half-sat value but even with an informative prior it will be estimated low. 
# Keep in mind for the gyr and ptarmigan

### Predator prey-only, fit of functional response (28/04/2015) but direct estimation (same model simulated as fitted). 
# small previous problem with Beverton-Holt formulation. Now Getz formulation, better. 

### For the functional response, the beta model is perhaps good for simulations, but I feel like we need something else
# to fit the model. Perhaps a logNormal. 

##################################################################################################################################

rm(list=ls())
graphics.off()

library("R2jags")      # Load R2jags package

### Parameters for simulation 

n.years<-100  	# Number of years 
N1<-1			# Initial pop size
P1<-0.1
K<-1			# threshold dd 
beta<-1			# density-dependence exponent
rmax_V<-2			# Max AVERAGE growth rate (thus not a true max...)
rmax_P<-0.5
sigma2.proc<-0.05		# worked well with 0.005 as well
# Process sigma on the log-scale, use the Peretti et al. value. 0.005


### FR and predator parameters
C<-2.5
D<-1
Q<-10

### Simulation of data
#set.seed(42) 
set.seed(41)
#set.seed(40)

y<-N<-P<-FR<-numeric(n.years)
N[1]<-N1
P[1]<-P1
FR[1]<-C*N[1]/(D+N[1])

rV<-rnorm(n.years-1,rmax_V,sqrt(sigma2.proc))
rP<-rnorm(n.years-1,rmax_P,sqrt(sigma2.proc))
FRnoise<-rnorm(n.years-1,0,sqrt(sigma2.proc))

for (t in 1:(n.years-1)){

  N[t+1]<-N[t]*(exp(rV[t])/(1+(N[t]/K)^beta))*exp(-FR[t]*P[t]/N[t])
  P[t+1]<-P[t]*exp(rP[t])/(1+P[t]*Q/N[t])
  FR[t+1]<-(C*N[t+1]/(D+N[t+1])) + FRnoise[t+1] #after for update
}
## Plotting time series of abundances and FR
par(mfrow=c(2,2))
plot(1:n.years,N,type="b")
plot(1:n.years,P,type="b")
#curve(dbeta(x,a,b),from=0, to=1)
plot(N,FR)
plot(log(N),log(P))
    
#seq(0,1,0.01)
# Bundle data
jags.data <- list(T=n.years,logN=log(N),logP=log(P),FR=FR)

sink("predprey.txt")
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
    tau_FR ~ dgamma(.01,.01)
    C~dgamma(.01,.01) # uninformative priors OK for that one
    D~dgamma(.01,.01)
    # check model Leslie to see how she specified priors...

    # Likelihood
    # state process

    for (t in 1:(T-1)){        

    FRUpdate[t] <- C*N[t]/(D+N[t]) #functional response equation, including noise
    FR[t] ~  dnorm(FRUpdate[t],tau_FR) #small trick to use FR data

    logNupdate[t] <- logN[t] + r_V -log(1+N[t]/K_V) -FR[t]*exp(logP[t])/N[t]
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
  list(sigma_V=runif(1,0.1,2), sigma_P=runif(1,0.1,2), r_V=runif(1,0.1,2),r_P=runif(1,0.1,2), K_V=runif(1,0.2,8), Q=runif(1,0,5),tau_FR=runif(1,1,10),C=runif(1,10,100),D=runif(1,0.01,0.1))}

# Parameters monitored
#parameters<-c("r_V","K_V","r_P","Q","sigma2_V","sigma2_P","a","b","C","D","logNupdate","logPupdate","FR")
parameters<-c("r_V","K_V","r_P","Q","sigma2_V","sigma2_P","tau_FR","C","D")


# MCMC settings
nc <- 3 #number of chains
nb <- 14000 # “burn in”
#ni <- 14000# “number of iterations” # that's for a symmetric distrib...
ni<-34000
nt <- 10 # “thinning”

# run model
out <- jags(jags.data, inits, parameters, "predprey.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, working.directory = getwd())
print(out, dig = 2)


# Output summary statistics
jags.sum<-out$BUGSoutput$summary
write.table(x=jags.sum,file="JAGSsummary_PredPreyStoch_analysisBis.txt")
# MCMC Output
#pdf("Output_MCMC__PredPreyStoch_analysisBis.pdf")
out.mcmc<-as.mcmc(out)
plot(out.mcmc) 
#dev.off()

CEb<-out$BUGSoutput$mean$C
DEb<-out$BUGSoutput$mean$D

# logN1<-out$BUGSoutput$mean$logNupdate
# logP1<-out$BUGSoutput$mean$logPupdate
# 
# par(mfrow=c(2,2))
# plot(1:(n.years-1),logN1,type="o")
# plot(1:(n.years-1),logP1,type="o")

### Fit functional response
par(mfrow=c(1,1))

fr_fit<-nls(FR~CE*N/(DE+N),start=list(CE=1,DE=1))
CE<-coef(fr_fit)[1]
DE<-coef(fr_fit)[2]
plot(N,FR,ylim=c(0,max(FR,na.rm=T)))
lines(N,CE*N/(DE+N))
lines(N,CEb*N/(DEb+N),col="blue")

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
    
    }
    
    }
    ",fill=TRUE)
sink()


# Initial values
inits <- function () {
  list(sigma_V=runif(1,0.1,2), sigma_P=runif(1,0.1,2), r_V=runif(1,0.1,2),r_P=runif(1,0.1,2), K_V=runif(1,0.2,8), Q=runif(1,0,5),C=runif(1,10,100),D=runif(1,0.01,0.1))}

# Parameters monitored
#parameters<-c("r_V","K_V","r_P","Q","sigma2_V","sigma2_P","a","b","C","D","logNupdate","logPupdate","FR")
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
print(out,dig=2) # 
# 01/06/2017 // Olivier says DIC is crap in this (and other?) context, avoid this...  
#plot(as.mcmc(out2)) 
#plot(as.mcmc(out)) 

#save.image('/media/frederic/DATA/Simuls_wOlivier/predatorpreyFRdata/simTS100.RData')

### plot densities
library(mcmcplots)
denplot(out,c("r_V","K_V","r_P","Q","sigma2_V","sigma2_P","C","D"))
denplot(out2,c("r_V","K_V","r_P","Q","sigma2_V","sigma2_P","C","D"))

### Trace plots
traplot(out,c("r_V","K_V","r_P","Q","sigma2_V","sigma2_P","C","D"))
traplot(out2,c("r_V","K_V","r_P","Q","sigma2_V","sigma2_P","C","D"))

## Prior-posterior overlap
library(MCMCvis) ## complicated 
par(mfrow=c(2,1),lwd=2)
C1 = out$BUGSoutput$sims.array[,1,'C']
C2 = out$BUGSoutput$sims.array[,2,'C']
C3 = out$BUGSoutput$sims.array[,3,'C']
curve(dgamma(x, shape = 0.01, rate = 0.01, log = FALSE),from=0.000001,to=5,ylim=c(0,3),ylab="Probability density",xlab = "C",main = "With FR data")
lines(density(C1,from=0.000001,to=5),col="green")
lines(density(C2,from=0.000001,to=5),col="blue")
lines(density(C3,from=0.000001,to=5),col="pink")
legend("topright",legend = c("prior","chain 1","chain 2","chain 3"),col=c("black","green","blue","pink"), lty = 1)

C1 = out2$BUGSoutput$sims.array[,1,'C']
C2 = out2$BUGSoutput$sims.array[,2,'C']
C3 = out2$BUGSoutput$sims.array[,3,'C']
curve(dgamma(x, shape = 0.01, rate = 0.01, log = FALSE),from=0.000001,to=5,ylim=c(0,3),ylab="Probability density",xlab = "C",main = "Without FR data")
lines(density(C1,bw=0.01,from=0.000001,to=5),col="green")
lines(density(C2,bw=0.01,from=0.000001,to=5),col="blue")
lines(density(C3,bw=0.01,from=0.000001,to=5),col="pink")
legend("topright",legend = c("prior","chain 1","chain 2","chain 3"),col=c("black","green","blue","pink"), lty = 1)
#PriorPosteriorOverlap_vagueGammaPriors.pdf

### Zoom 
par(mfrow=c(2,1),lwd=2)
C1 = out$BUGSoutput$sims.array[,1,'C']
C2 = out$BUGSoutput$sims.array[,2,'C']
C3 = out$BUGSoutput$sims.array[,3,'C']
curve(dgamma(x, shape = 0.01, rate = 0.01, log = FALSE),from=0.000001,to=5,ylim=c(0,1),ylab="Probability density",xlab = "C",main = "With FR data")
lines(density(C1,from=0.000001,to=5),col="green")
lines(density(C2,from=0.000001,to=5),col="blue")
lines(density(C3,from=0.000001,to=5),col="pink")
legend("topright",legend = c("prior","chain 1","chain 2","chain 3"),col=c("black","green","blue","pink"), lty = 1)

C1 = out2$BUGSoutput$sims.array[,1,'C']
C2 = out2$BUGSoutput$sims.array[,2,'C']
C3 = out2$BUGSoutput$sims.array[,3,'C']
curve(dgamma(x, shape = 0.01, rate = 0.01, log = FALSE),from=0.000001,to=5,ylim=c(0,1),ylab="Probability density",xlab = "C",main = "Without FR data")
lines(density(C1,bw=0.01,from=0.000001,to=5),col="green")
lines(density(C2,bw=0.01,from=0.000001,to=5),col="blue")
lines(density(C3,bw=0.01,from=0.000001,to=5),col="pink")
legend("topright",legend = c("prior","chain 1","chain 2","chain 3"),col=c("black","green","blue","pink"), lty = 1)


### Plot pair posterior densities
postsamples=cbind(out$BUGSoutput$sims.list$r_V,
out$BUGSoutput$sims.list$K_V,
out$BUGSoutput$sims.list$r_P,
out$BUGSoutput$sims.list$Q,
out$BUGSoutput$sims.list$C,
out$BUGSoutput$sims.list$D)
png(file="PairPosteriorPlot_withFRdata.png", width = 1200, height = 1200,res=300)
pairs(postsamples,c("r_V","K_V","r_P","Q","C","D"))
dev.off()

postsamples2=cbind(out2$BUGSoutput$sims.list$r_V,
                  out2$BUGSoutput$sims.list$K_V,
                  out2$BUGSoutput$sims.list$r_P,
                  out2$BUGSoutput$sims.list$Q,
                  out2$BUGSoutput$sims.list$C,
                  out2$BUGSoutput$sims.list$D)
png(file="PairPosteriorPlot_withoutFRdata.png", width = 1200, height = 1200, res=300)
pairs(postsamples2,c("r_V","K_V","r_P","Q","C","D"))
dev.off()

pdf(file="PairCorrelPosteriorPlot.pdf",width = 4,height = 8)
par(mfrow=c(2,1))
parcorplot(out,parms = c("r_V","K_V","r_P","Q","sigma2_V","sigma2_P","C","D"))
parcorplot(out2,parms = c("r_V","K_V","r_P","Q","sigma2_V","sigma2_P","C","D"))
dev.off()

# Output summary statistics
jags.sum2<-out2$BUGSoutput$summary
write.table(x=jags.sum2,file="JAGSsummary_PredPreyStoch_analysisBis_noFR.txt")

### Predicted values of the model fitted without the functional response
logN2<-out2$BUGSoutput$mean$logNupdate
logP2<-out2$BUGSoutput$mean$logPupdate


par(mfrow=c(2,1))
plot(2:(n.years),log(N)[2:(n.years)],type="o")
lines(1:(n.years-1),logN1,type="o",col="red")
lines(1:(n.years-1),logN2,type="o",col="blue")
plot(2:(n.years),log(P)[2:(n.years)],type="o")
lines(1:(n.years-1),logP1,type="o",col="red")
lines(1:(n.years-1),logP2,type="o",col="blue")



#### Reproduce the FR curve without and without the correlations between parameters
## All possible curves
Clist = out$BUGSoutput$sims.array[,,'C']
Dlist = out$BUGSoutput$sims.array[,,'D']
n = length(Clist)
ndens = 100
Nprey <- seq(1,20,length=ndens) #density index
FRstoch = matrix(NA,nrow = n, ncol = ndens)

png('Estimated_FR.png',res=300,width=2000,height=1000)
par(mfrow=c(1,2))
library(scales)
for (i in 1:n){
  for (dens in 1:length(Nprey))
  {
    FRstoch[i,dens] = Clist[i]*Nprey[dens]/(Dlist[i]+Nprey[dens])
  }
  if (i == 1){plot(Nprey,FRstoch[i,],type='l',lwd=0.5,col=alpha('black',0.05),ylab='Functional response',ylim=c(1,3),xlim=c(1,10),xlab='N prey',main='With (C,D) correlations')
  }
  else {lines(Nprey,FRstoch[i,],type='l',lwd=0.5,col=alpha('black',0.05))}
}

# what if there was no correlation between parameters? 
# randomize
Clist = sample(Clist)
Dlist = sample(Dlist)
for (i in 1:n){
  for (dens in 1:length(Nprey))
  {
    FRstoch[i,dens] = Clist[i]*Nprey[dens]/(Dlist[i]+Nprey[dens])
  }
  if (i == 1){plot(Nprey,FRstoch[i,],type='l',lwd=0.5,col=alpha('black',0.05),ylab='Functional response',ylim=c(1,3),xlim=c(1,10),xlab='N prey',main='Without (C,D) correlations')
  }
  else {lines(Nprey,FRstoch[i,],type='l',lwd=0.5,col=alpha('black',0.05))}
}

dev.off()

#### Same plot without the FR data

#Clist = out2$BUGSoutput$sims.list$C
#Dlist = out2$BUGSoutput$sims.list$D
Clist = out2$BUGSoutput$sims.array[,,'C']
Dlist = out2$BUGSoutput$sims.array[,,'D']
n = length(Clist)
ndens = 100
Nprey <- seq(1,20,length=ndens) #density index
FRstoch = matrix(NA,nrow = n, ncol = ndens)

png('Estimated_FR_withoutFRdata.png',res=300,width=2000,height=1000)
par(mfrow=c(1,2))
library(scales)
for (i in 1:n){
  for (dens in 1:length(Nprey))
  {
    FRstoch[i,dens] = Clist[i]*Nprey[dens]/(Dlist[i]+Nprey[dens])
  }
  if (i == 1){plot(Nprey,FRstoch[i,],type='l',lwd=0.5,col=alpha('black',0.05),ylab='Functional response',ylim=c(1,3),xlim=c(1,10),xlab='N prey',main='With (C,D) correlations')
  }
  else {lines(Nprey,FRstoch[i,],type='l',lwd=0.5,col=alpha('black',0.05))}
}

# what if there was no correlation between parameters? 

Clist = sample(Clist)
Dlist = sample(Dlist)
for (i in 1:n){
  for (dens in 1:length(Nprey))
  {
    FRstoch[i,dens] = Clist[i]*Nprey[dens]/(Dlist[i]+Nprey[dens])
  }
  if (i == 1){plot(Nprey,FRstoch[i,],type='l',lwd=0.5,col=alpha('black',0.05),ylab='Functional response',ylim=c(1,3),xlim=c(1,10),xlab='N prey',main='Without (C,D) correlations')
  }
  else {lines(Nprey,FRstoch[i,],type='l',lwd=0.5,col=alpha('black',0.05))}
}

dev.off()


#### Reproduce the density-dependent curves without and without the correlations between parameters
# For the prey
#rlist = out$BUGSoutput$sims.list$r_V
#Klist = out$BUGSoutput$sims.list$K_V
rlist = out$BUGSoutput$sims.array[,,'r_V']
Klist = out$BUGSoutput$sims.array[,,'K_V']

n = length(rlist)
ndens = 100
Nprey <- seq(1,50,length=ndens) #density index
preyDD = matrix(NA,nrow = n, ncol = ndens)

png('Estimated_preyDD.png',res=300,width=2000,height=1000)
par(mfrow=c(1,2))
library(scales)
for (i in 1:n){
  for (dens in 1:length(Nprey))
  {
    preyDD[i,dens] = exp(rlist[i])/(1+Nprey[dens]/Klist[i])
  }
  if (i == 1){plot(Nprey,preyDD[i,],type='l',lwd=0.5,col=alpha('black',0.05),ylab='Prey growth rate',ylim=c(-1,6),xlim=c(1,50),xlab='N prey',main='With (r,K) correlations')
  }
  else {lines(Nprey,preyDD[i,],type='l',lwd=0.5,col=alpha('black',0.05))}
}

# what if there was no correlation between parameters? 
rlist = sample(rlist)
Klist = sample(Klist)
for (i in 1:n){
  for (dens in 1:length(Nprey))
  {
    preyDD[i,dens] = exp(rlist[i])/(1+Nprey[dens]/Klist[i])
  }
  if (i == 1){plot(Nprey,preyDD[i,],type='l',lwd=0.5,col=alpha('black',0.05),ylab='Prey growth rate',ylim=c(-1,6),xlim=c(1,50),xlab='N prey',main='Without (r,K) correlations')
  }
  else {lines(Nprey,preyDD[i,],type='l',lwd=0.5,col=alpha('black',0.05))}
}

dev.off()

## without FR data

# For the prey
#rlist = out2$BUGSoutput$sims.list$r_V
#Klist = out2$BUGSoutput$sims.list$K_V
rlist = out2$BUGSoutput$sims.array[,,'r_V']
Klist = out2$BUGSoutput$sims.array[,,'K_V']

n = length(rlist)
ndens = 100
Nprey <- seq(1,50,length=ndens) #density index
preyDD = matrix(NA,nrow = n, ncol = ndens)

png('Estimated_preyDD_withoutFRdata.png',res=300,width=2000,height=1000)
par(mfrow=c(1,2))
library(scales)
for (i in 1:n){
  for (dens in 1:length(Nprey))
  {
    preyDD[i,dens] = exp(rlist[i])/(1+Nprey[dens]/Klist[i])
  }
  if (i == 1){plot(Nprey,preyDD[i,],type='l',lwd=0.5,col=alpha('black',0.05),ylab='Prey growth rate',ylim=c(-1,6),xlim=c(1,50),xlab='N prey',main='With (r,K) correlations')
  }
  else {lines(Nprey,preyDD[i,],type='l',lwd=0.5,col=alpha('black',0.05))}
}

# what if there was no correlation between parameters? 
rlist = sample(rlist)
Klist = sample(Klist)
for (i in 1:n){
  for (dens in 1:length(Nprey))
  {
    preyDD[i,dens] = exp(rlist[i])/(1+Nprey[dens]/Klist[i])
  }
  if (i == 1){plot(Nprey,preyDD[i,],type='l',lwd=0.5,col=alpha('black',0.05),ylab='Prey growth rate',ylim=c(-1,6),xlim=c(1,50),xlab='N prey',main='Without (r,K) correlations')
  }
  else {lines(Nprey,preyDD[i,],type='l',lwd=0.5,col=alpha('black',0.05))}
}

dev.off()



