### F. Barraquand 21/04/2015 - Code for analyzing noisy predator-prey system (incl. noise on the functional response)
### Updated 01/06/2017. 
### Case with "small noise" on the functional response
### FB 03/01/2019 Added plots for the paired posterior distributions
### FB 05/01/2019 Reparameterized to check whether it removes the Paired PD correlations
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

n.years<-1000  	# Number of years - 25 first, perhaps use 50 or 100 / worked very well with almost no process error on the real scale
N1<-1			# Initial pop size
P1<-0.1

beta<-1			# density-dependence exponent
rmax_V<-2			# Max AVERAGE growth rate (thus not a true max...)
K<-1*(exp(rmax_V)-1)	# true carrying capacity - so that threshold for dd = 1 as in previous code

rmax_P<-0.5
sigma2.proc<-0.05		# worked well with 0.005
# Process sigma on the log-scale, use the Peretti et al. value. 0.005


### FR and predator parameters
C<-2.5
D<-1
a = C/D #reparam to classical Holling type II
h = 1/C #reparam, handling time

### Reparam, we have 
# Q*(exp(r_P)-1) should be 10 to keep the exact same parameter values as the previous model, irrespective of parameterization
Q = 10/(exp(rmax_P)-1)

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
 
  N[t+1]<-N[t]*(exp(rV[t]) / (1+ N[t]*(exp(rmax_V)-1)/K ) )*exp(-FR[t]*P[t]/N[t])
  P[t+1]<-P[t]*exp(rP[t]) /(1+P[t]*(exp(rmax_P)-1)*Q/N[t] )
  FR[t+1]<- C*N[t+1]/(D+N[t+1])  + FRnoise[t+1]
}

## Plotting time series of abundances and FR
par(mfrow=c(2,2))
plot(1:n.years,N,type="b")
plot(1:n.years,P,type="b")
#curve(dbeta(x,a,b),from=0, to=1)
plot(N,FR)

    
#seq(0,1,0.01)
# Bundle data
jags.data <- list(T=n.years,logN=log(N),logP=log(P),FR=FR)

sink("predpreymod.txt")
cat("
    model {
    
    # Priors and constraints
    logN[1] ~ dnorm(0,0.01) # Prior on initial pop size on the log scale
    logP[1] ~ dnorm(0,0.01) # Prior on initial pop size on the log scale

    # Priors prey population dynamics
    r_V ~ dnorm(1,1) #dnorm(1,0.001) # below the truth, rather flat prior
    K_V ~ dunif(0.2,10)
    sigma_V ~ dunif(0.01,5) # rather vague 
    sigma2_V<-pow(sigma_V, 2)
    tau_V<-pow(sigma_V,-2)

 
    #Priors predator population dynamics
    Q ~ dgamma(0.1,0.1) ### can I change by e.g. dunif(4,40)?
    r_P ~ dnorm(0.3,1) ## dnorm(1,0.1)
    sigma_P ~ dunif(0.01,2) # rather vague 
    sigma2_P<-pow(sigma_P, 2)
    tau_P<-pow(sigma_P,-2)
    
    #Priors predation parameters 
    tau_FR ~ dgamma(.01,.01)
    #C~dgamma(.01,.01) # uninformative priors OK for that one
    #D~dgamma(.01,.01)
    a~dgamma(.01,.01) # uninformative priors OK for that one
    h~dgamma(.01,.01)

    # Intermediate nodes
    MV <- K_V / max(exp(r_V)-1,0.1)
    #Another option is to restrict r_V and r_P to positive values 
    QP <- Q * max(exp(r_P)-1,0.1) #this damned worked when I made a typo and used r_V
   

    # Likelihood
    # state process

    for (t in 1:(T-1)){        

    FRUpdate[t] <- a*N[t]/(1+a*h*N[t]) #functional response equation, including noise
    FR[t] ~  dnorm(FRUpdate[t],tau_FR) #small trick to use FR data

    logNupdate[t] <- logN[t] + r_V -log(1+ N[t]/MV ) -FR[t]*exp(logP[t])/N[t]
    logN[t+1] ~ dnorm(logNupdate[t],tau_V)
    N[t]<-exp(logN[t])

    # for some reason, log(1+(exp(r_V)-1)*N[t]/K_V) was not working initially. 
    # Likely the problem was when r_V became negative in some of the chains

    logP[t+1]~ dnorm(logPupdate[t],tau_P)
    logPupdate[t] <- logP[t] + r_P - log( 1 + exp(logP[t])*QP/N[t] ) 
    
    }
    
    }
    ",fill=TRUE)
sink()


# Initial values
inits <- function () {
  list(sigma_V=runif(1,0.1,2), sigma_P=runif(1,0.1,2), r_V=runif(1,1,2),r_P=runif(1,0.3,0.8), K_V=runif(1,0.2,8), Q=runif(1,4,10),tau_FR=runif(1,1,10),a=runif(1,1,10),h=runif(1,0.1,1))}

# Parameters monitored
#parameters<-c("r_V","K_V","r_P","Q","sigma2_V","sigma2_P","a","b","C","D","logNupdate","logPupdate","FR")
parameters<-c("r_V","K_V","r_P","Q","sigma2_V","sigma2_P","tau_FR","a","h")


# MCMC settings
nc <- 3 #number of chains
nb <- 14000 # “burn in”
#ni <- 14000# “number of iterations” # that's for a symmetric distrib...
ni<-34000
nt <- 10 # “thinning”

# run model
out <- jags(jags.data, inits, parameters, "predpreymod.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, working.directory = getwd())
print(out, dig = 2)

# Output summary statistics
jags.sum<-out$BUGSoutput$summary
write.table(x=jags.sum,file="JAGSsummary_PredPreyStoch_analysisBis_reparam.txt")
# MCMC Output
#pdf("Output_MCMC__PredPreyStoch_analysisBis.pdf")
out.mcmc<-as.mcmc(out)
plot(out.mcmc) 
#dev.off()

aEb<-out$BUGSoutput$mean$a
hEb<-out$BUGSoutput$mean$h

logN1<-out$BUGSoutput$mean$logNupdate
logP1<-out$BUGSoutput$mean$logPupdate

par(mfrow=c(2,2))
plot(1:(n.years-1),logN1,type="o")
plot(1:(n.years-1),logP1,type="o")

### Fit functional response

fr_fit<-nls(FR~aE*N/(1+aE*hE*N),start=list(aE=1,hE=1))
aE<-coef(fr_fit)[1]
hE<-coef(fr_fit)[2]
plot(N,FR,ylim=c(0,max(FR,na.rm=T)))
lines(N,aE*N/(1+aE*hE*N))
lines(N,aEb*N/(1+aEb*hEb*N),col="blue")

### Now try to fit a model without the FR data. 
sink("predpreymod_without_sepFR.txt")
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
    
    #C~dgamma(.01,.01) # 
    # D~dgamma(.01,.01)
    a~dgamma(.01,.01) # 
    h~dgamma(.01,.01)

    # Intermediate nodes
    MV <- K_V / max(exp(r_V)-1,0.1)
    #Another option is to restrict r_V and r_P to positive values 
    QP <- Q * max(exp(r_P)-1,0.1) #this damned worked when I made a typo and used r_V


    # Likelihood
    # state process
    
    for (t in 1:(T-1)){        

    
    logNupdate[t] <- logN[t] + r_V -log(1+N[t]/MV) -a*exp(logP[t])/(1+a*h*N[t])
    logN[t+1] ~ dnorm(logNupdate[t],tau_V)
    N[t]<-exp(logN[t])
    # for some reason, log(1+(exp(r_V)-1)*N[t]/K_V) was not working

    logP[t+1]~ dnorm(logPupdate[t],tau_P)
    logPupdate[t] <- logP[t] + r_P - log(1+exp(logP[t])*QP/exp(logN[t]) )  
    
    }
    
    }
    ",fill=TRUE)
sink()


# Initial values
inits <- function () {
  list(sigma_V=runif(1,0.1,2), sigma_P=runif(1,0.1,2), r_V=runif(1,0.1,2),r_P=runif(1,0.1,2), K_V=runif(1,0.2,8), Q=runif(1,1,10),a=runif(1,1,10),h=runif(1,0.1,1))}

# Parameters monitored
parameters<-c("r_V","K_V","r_P","Q","sigma2_V","sigma2_P","tau_FR","a","h")


# MCMC settings
nc <- 3 #number of chains
nb <- 14000 # “burn in”
#ni <- 14000# “number of iterations” # that's for a symmetric distrib...
ni<-34000
nt <- 10 # “thinning”

# run model
out2 <- jags(jags.data, inits, parameters, "predpreymod_without_sepFR.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, working.directory = getwd())
print(out2, dig = 2)
# http://jeromyanglim.tumblr.com/post/37362047458/how-to-get-dic-in-jags
print(out,dig=2) # to compare deviance, DIC - check also parameter values.

save.image('/media/frederic/DATA/Simuls_wOlivier/predatorpreyFRdata/simlongTS_reparam.RData')
file_path = "/media/frederic/DATA/Simuls_wOlivier/predatorpreyFRdata/"
save(out,file = paste(file_path,'out_simlongTS_reparam.RData',sep=""))
save(out2,file = paste(file_path,'out2_simlongTS_reparam.RData',sep=""))

# 01/06/2017 // Olivier says DIC is crap in this (and other?) context, avoid this...  
plot(as.mcmc(out2)) 
plot(as.mcmc(out)) 

### plot densities
library(mcmcplots)
denplot(out,c("r_V","K_V","r_P","Q","sigma2_V","sigma2_P","a","h"))
denplot(out2,c("r_V","K_V","r_P","Q","sigma2_V","sigma2_P","a","h"))

### Trace plots
png(file="TracePlot_withFRdata_reparam.png", width = 1200, height = 1200,res=300)
traplot(out,c("r_V","K_V","r_P","Q","sigma2_V","sigma2_P","a","h"))
dev.off()
png(file="TracePlot_withoutFRdata_reparam.png", width = 1200, height = 1200,res=300)
traplot(out2,c("r_V","K_V","r_P","Q","sigma2_V","sigma2_P","a","h"))
dev.off()

### Plot pair posterior densities
postsamples=cbind(out$BUGSoutput$sims.list$r_V,
out$BUGSoutput$sims.list$K_V,
out$BUGSoutput$sims.list$r_P,
out$BUGSoutput$sims.list$Q,
out$BUGSoutput$sims.list$a,
out$BUGSoutput$sims.list$h)
png(file="PairPosteriorPlot_withFRdata_reparam.png", width = 1200, height = 1200,res=300)
pairs(postsamples,c("r_V","K_V","r_P","Q","a","h"))
dev.off()

postsamples2=cbind(out2$BUGSoutput$sims.list$r_V,
                  out2$BUGSoutput$sims.list$K_V,
                  out2$BUGSoutput$sims.list$r_P,
                  out2$BUGSoutput$sims.list$Q,
                  out2$BUGSoutput$sims.list$a,
                  out2$BUGSoutput$sims.list$h)
png(file="PairPosteriorPlot_withoutFRdata_reparam.png", width = 1200, height = 1200, res=300)
pairs(postsamples2,c("r_V","K_V","r_P","Q","a","h"))
dev.off()

pdf(file="PairCorrelPosteriorPlot_reparam.pdf",width = 4,height = 8)
par(mfrow=c(2,1))
parcorplot(out,parms = c("r_V","K_V","r_P","Q","sigma2_V","sigma2_P","a","h"))
parcorplot(out2,parms = c("r_V","K_V","r_P","Q","sigma2_V","sigma2_P","a","h"))
dev.off()

# Output summary statistics
jags.sum2<-out2$BUGSoutput$summary
write.table(x=jags.sum2,file="JAGSsummary_PredPreyStoch_analysisBis_noFR_reparam.txt")

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
alist = out$BUGSoutput$sims.list$a
hlist = out$BUGSoutput$sims.list$h
n = length(Clist)
ndens = 100
Nprey <- seq(1,20,length=ndens) #density index
FRstoch = matrix(NA,nrow = n, ncol = ndens)

png('Estimated_FR_reparam.png',res=300,width=2000,height=1000)
par(mfrow=c(1,2))
library(scales)
for (i in 1:n){
  for (dens in 1:length(Nprey))
  {
    FRstoch[i,dens] = alist[i]*Nprey[dens]/(1.0+alist[i]*hlist[i]*Nprey[dens])
  }
  if (i == 1){plot(Nprey,FRstoch[i,],type='l',lwd=0.5,col=alpha('black',0.05),ylab='Functional response',ylim=c(1,3),xlim=c(1,10),xlab='N prey',main='With (a,h) correlations')
  }
  else {lines(Nprey,FRstoch[i,],type='l',lwd=0.5,col=alpha('black',0.05))}
  lines(Nprey,a*Nprey/(1+a*h*Nprey),col=alpha('red',1.0))
}

# what if there was no correlation between parameters? 
# randomize
alist = sample(alist)
hlist = sample(hlist)
for (i in 1:n){
  for (dens in 1:length(Nprey))
  {
    FRstoch[i,dens] = alist[i]*Nprey[dens]/(1.0+alist[i]*hlist[i]*Nprey[dens])
   
  }
  if (i == 1){plot(Nprey,FRstoch[i,],type='l',lwd=0.5,col=alpha('black',0.05),ylab='Functional response',ylim=c(1,3),xlim=c(1,10),xlab='N prey',main='Without (a,h) correlations')
  }
  else {lines(Nprey,FRstoch[i,],type='l',lwd=0.5,col=alpha('black',0.05))}
  lines(Nprey,a*Nprey/(1+a*h*Nprey),col=alpha('red',1.0))
}

dev.off()

#### Same plot without the FR data

alist = out2$BUGSoutput$sims.list$a
hlist = out2$BUGSoutput$sims.list$h
n = length(Clist)
ndens = 100
Nprey <- seq(1,20,length=ndens) #density index
FRstoch = matrix(NA,nrow = n, ncol = ndens)

png('Estimated_FR_withoutFRdata_reparam.png',res=300,width=2000,height=1000)
par(mfrow=c(1,2))
library(scales)
for (i in 1:n){
  for (dens in 1:length(Nprey))
  {
    FRstoch[i,dens] = alist[i]*Nprey[dens]/(1.0+alist[i]*hlist[i]*Nprey[dens])
  }
  if (i == 1){plot(Nprey,FRstoch[i,],type='l',lwd=0.5,col=alpha('black',0.05),ylab='Functional response',ylim=c(1,3),xlim=c(1,10),xlab='N prey',main='With (a,h) correlations')
  }
  else {lines(Nprey,FRstoch[i,],type='l',lwd=0.5,col=alpha('black',0.05))}
  lines(Nprey,a*Nprey/(1+a*h*Nprey),col=alpha('red',1.0))
}

# what if there was no correlation between parameters? 

alist = sample(alist)
hlist = sample(hlist)
for (i in 1:n){
  for (dens in 1:length(Nprey))
  {
    FRstoch[i,dens] = alist[i]*Nprey[dens]/(1.0+alist[i]*hlist[i]*Nprey[dens])
  }
  if (i == 1){plot(Nprey,FRstoch[i,],type='l',lwd=0.5,col=alpha('black',0.05),ylab='Functional response',ylim=c(1,3),xlim=c(1,10),xlab='N prey',main='Without (a,h) correlations')
  }
  else {lines(Nprey,FRstoch[i,],type='l',lwd=0.5,col=alpha('black',0.05))}
  lines(Nprey,a*Nprey/(1+a*h*Nprey),col=alpha('red',1.0))
}

dev.off()


#### Reproduce the density-dependent curves without and without the correlations between parameters
# For the prey
rlist = out$BUGSoutput$sims.list$r_V
Klist = out$BUGSoutput$sims.list$K_V

n = length(rlist)
ndens = 100
Nprey <- seq(1,50,length=ndens) #density index
preyDD = matrix(NA,nrow = n, ncol = ndens)

png('Estimated_preyDD_reparam.png',res=300,width=2000,height=1000)
par(mfrow=c(1,2))
library(scales)
for (i in 1:n){
  for (dens in 1:length(Nprey))
  {
    preyDD[i,dens] = exp(rlist[i])/(1+(exp(rlist[i])-1)*Nprey[dens]/Klist[i])
  }
  if (i == 1){plot(Nprey,preyDD[i,],type='l',lwd=0.5,col=alpha('black',0.05),ylab='Prey growth rate',ylim=c(-1,6),xlim=c(1,50),xlab='N prey',main='With (r,K) correlations')
  }
  else {lines(Nprey,preyDD[i,],type='l',lwd=0.5,col=alpha('black',0.05))}
  lines(Nprey,exp(rmax_V)/(1+Nprey*(exp(rmax_V)-1)/K),col=alpha('red',1.0))
}

# what if there was no correlation between parameters? 
rlist = sample(rlist)
Klist = sample(Klist)
for (i in 1:n){
  for (dens in 1:length(Nprey))
  {
    preyDD[i,dens] = exp(rlist[i])/(1+(exp(rlist[i])-1)*Nprey[dens]/Klist[i])
  }
  if (i == 1){plot(Nprey,preyDD[i,],type='l',lwd=0.5,col=alpha('black',0.05),ylab='Prey growth rate',ylim=c(-1,6),xlim=c(1,50),xlab='N prey',main='Without (r,K) correlations')
  }
  else {lines(Nprey,preyDD[i,],type='l',lwd=0.5,col=alpha('black',0.05))}
  lines(Nprey,exp(rmax_V)/(1+Nprey*(exp(rmax_V)-1)/K),col=alpha('red',1.0))
}

dev.off()

## without FR data

# For the prey
rlist = out2$BUGSoutput$sims.list$r_V
Klist = out2$BUGSoutput$sims.list$K_V

n = length(rlist)
ndens = 100
Nprey <- seq(1,50,length=ndens) #density index
preyDD = matrix(NA,nrow = n, ncol = ndens)

png('Estimated_preyDD_withoutFRdata_reparam.png',res=300,width=2000,height=1000)
par(mfrow=c(1,2))
library(scales)
for (i in 1:n){
  for (dens in 1:length(Nprey))
  {
    preyDD[i,dens] = exp(rlist[i])/(1+(exp(rlist[i])-1)*Nprey[dens]/Klist[i])
  }
  if (i == 1){plot(Nprey,preyDD[i,],type='l',lwd=0.5,col=alpha('black',0.05),ylab='Prey growth rate',ylim=c(-1,6),xlim=c(1,50),xlab='N prey',main='With (r,K) correlations')
  }
  else {lines(Nprey,preyDD[i,],type='l',lwd=0.5,col=alpha('black',0.05))}
  lines(Nprey,exp(rmax_V)/(1+Nprey*(exp(rmax_V)-1)/K),col=alpha('red',1.0))
}

# what if there was no correlation between parameters? 
rlist = sample(rlist)
Klist = sample(Klist)
for (i in 1:n){
  for (dens in 1:length(Nprey))
  {
    preyDD[i,dens] = exp(rlist[i])/(1+(exp(rlist[i])-1)*Nprey[dens]/Klist[i])
  }
  if (i == 1){plot(Nprey,preyDD[i,],type='l',lwd=0.5,col=alpha('black',0.05),ylab='Prey growth rate',ylim=c(-1,6),xlim=c(1,50),xlab='N prey',main='Without (r,K) correlations')
  }
  else {lines(Nprey,preyDD[i,],type='l',lwd=0.5,col=alpha('black',0.05))}
  lines(Nprey,exp(rmax_V)/(1+Nprey*(exp(rmax_V)-1)/K),col=alpha('red',1.0))
}

dev.off()


