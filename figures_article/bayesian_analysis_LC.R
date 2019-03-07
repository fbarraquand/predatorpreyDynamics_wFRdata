
########## Bayesian analysis of the predator-prey system with and without kill rate data ####################

### F. Barraquand 21/04/2015 - Code for analyzing noisy predator-prey system (incl. noise on the functional response)
### Updated 01/06/2017. 
### Case with "small noise" on the functional response
### FB 03/01/2019 Added plots for the paired posterior distributions
### FB 14/01/2019 Corrected FR to nondelayed version

rm(list=ls())
graphics.off()

library("R2jags")      # Load R2jags package

### Parameters for simulation 

n.years<-100  	# Number of years 
N1<-1			# Initial pop size
P1<-0.1
gamma<-1			# threshold dd 
rmax_V<-2			# Max AVERAGE growth rate (thus not a true max...)
rmax_P<-0.5
sigma2.proc<-0.05		# worked well with 0.005 as well
# Process sigma on the log-scale, use the Peretti et al. value. 0.005

### FR and predator parameters
C<-15
D<-0.25
Q<-10

# Loading data
file_path1 = "../simulations/small_noise_on_FR/parameter_sets/noisy_cycles/T=100/"
krep=1 # sample chosen
data1 = read.csv(paste(file_path1,"predatorPrey_withGaussianFR",krep,".csv",sep=""))
names(data1)=c("Time","n","p","KR")
head(data1)

# Bundle data
jags.data <- list(T=n.years,logN=data1$n,logP=data1$p,FR=data1$KR)

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


# Initial values
inits <- function () {
  list(sigma_V=runif(1,0.1,2), sigma_P=runif(1,0.1,2), r_V=runif(1,0.1,2),r_P=runif(1,0.1,2), gamma=runif(1,0.3,4), Q=runif(1,0,5),tau_FR=runif(1,1,10),C=runif(1,10,100),D=runif(1,0.01,0.1))}

# Parameters monitored
parameters<-c("r_V","gamma","r_P","Q","sigma2_V","sigma2_P","tau_FR","C","D")


# MCMC settings
nc <- 3 #number of chains
nb <- 24000 # “burn in”
#ni <- 14000# “number of iterations” # that's for a symmetric distrib...
ni<-44000
nt <- 10 # “thinning”

# run model
out <- jags(jags.data, inits, parameters, "predprey.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, working.directory = getwd())
print(out, dig = 2)
# Output summary statistics
#jags.sum<-out$BUGSoutput$summary
#write.table(x=jags.sum,file="JAGSsummary_PredPreyStoch_analysisBis.txt")
# Convergence a little less good than with T=1000

# MCMC Output
#pdf("Output_MCMC__PredPreyStoch.pdf")
out.mcmc<-as.mcmc(out)
plot(out.mcmc) 
#dev.off()

CEb<-out$BUGSoutput$mean$C
DEb<-out$BUGSoutput$mean$D

# logN1<-out$BUGSoutput$mean$logN
# logP1<-out$BUGSoutput$mean$logP
# 
# par(mfrow=c(2,2))
# plot(1:(n.years-1),logN1,type="o")
# plot(1:(n.years-1),logP1,type="o")

### Fit functional response
par(mfrow=c(1,1))

N=exp(data1$n)
FR=data1$KR
fr_fit<-nls(FR~CE*N/(DE+N),start=list(CE=10,DE=0.05))
CE<-coef(fr_fit)[1]
DE<-coef(fr_fit)[2]
plot(N,FR,ylim=c(0,max(FR,na.rm=T)))
curve(CE*x/(DE+x),add=TRUE)
curve(CEb*x/(DEb+x),col="blue",add=TRUE)
### Equivalent to least-square fit. 

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
  list(sigma_V=runif(1,0.1,2), sigma_P=runif(1,0.1,2), r_V=runif(1,0.1,2),r_P=runif(1,0.1,2),gamma=runif(1,0.3,4), Q=runif(1,0,5),C=runif(1,10,100),D=runif(1,0.01,0.1))}

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
print(out,dig=2) # 

#save.image('/media/frederic/DATA/Simuls_wOlivier/predatorpreyFRdata/simTS100.RData')

### plot densities
library(mcmcplots)
denplot(out,c("r_V","gamma","r_P","Q","sigma2_V","sigma2_P","C","D"))
denplot(out2,c("r_V","gamma","r_P","Q","sigma2_V","sigma2_P","C","D"))

### Trace plots
traplot(out,c("r_V","gamma","r_P","Q","sigma2_V","sigma2_P","C","D"))
traplot(out2,c("r_V","gamma","r_P","Q","sigma2_V","sigma2_P","C","D"))

## Prior-posterior overlap

par(mfrow=c(2,1),lwd=2)
C1 = out$BUGSoutput$sims.array[,1,'C']
C2 = out$BUGSoutput$sims.array[,2,'C']
C3 = out$BUGSoutput$sims.array[,3,'C']
curve(dgamma(x, shape = 0.01, rate = 0.01, log = FALSE),from=0.000001,to=20,ylim=c(0,1),ylab="Probability density",xlab = "C",main = "With KR data")
lines(density(C1,from=0.000001,to=20),col="green")
lines(density(C2,from=0.000001,to=20),col="blue")
lines(density(C3,from=0.000001,to=20),col="pink")
legend("topleft",legend = c("prior","chain 1","chain 2","chain 3"),col=c("black","green","blue","pink"), lty = 1)

C1 = out2$BUGSoutput$sims.array[,1,'C']
C2 = out2$BUGSoutput$sims.array[,2,'C']
C3 = out2$BUGSoutput$sims.array[,3,'C']
curve(dgamma(x, shape = 0.01, rate = 0.01, log = FALSE),from=0.000001,to=20,ylim=c(0,1),ylab="Probability density",xlab = "C",main = "Without KR data")
lines(density(C1,bw=0.01,from=0.000001,to=20),col="green")
lines(density(C2,bw=0.01,from=0.000001,to=20),col="blue")
lines(density(C3,bw=0.01,from=0.000001,to=20),col="pink")
legend("topleft",legend = c("prior","chain 1","chain 2","chain 3"),col=c("black","green","blue","pink"), lty = 1)
#PriorPosteriorOverlap_vagueGammaPriors.pdf

### With both C and D

pdf(file="PPO_LC.pdf",width=12,height=8)
par(mfrow=c(2,2),lwd=2,cex=1.2)

C1 = out$BUGSoutput$sims.array[,1,'C']
C2 = out$BUGSoutput$sims.array[,2,'C']
C3 = out$BUGSoutput$sims.array[,3,'C']
curve(dgamma(x, shape = 0.01, rate = 0.01, log = FALSE),from=0.000001,to=20,ylim=c(0,2),ylab="Probability density",xlab = "C",main = "With KR data")
lines(density(C1,from=0.000001,to=20),col="green")
lines(density(C2,from=0.000001,to=20),col="blue")
lines(density(C3,from=0.000001,to=20),col="pink")
abline(v=C,col="red")
legend("top",legend = c("prior","chain 1","chain 2","chain 3"),col=c("black","green","blue","pink"), lty = 1)

C1 = out2$BUGSoutput$sims.array[,1,'C']
C2 = out2$BUGSoutput$sims.array[,2,'C']
C3 = out2$BUGSoutput$sims.array[,3,'C']
curve(dgamma(x, shape = 0.01, rate = 0.01, log = FALSE),from=0.000001,to=20,ylim=c(0,2),ylab="Probability density",xlab = "C",main = "Without KR data")
lines(density(C1,bw=0.01,from=0.000001,to=20),col="green")
lines(density(C2,bw=0.01,from=0.000001,to=20),col="blue")
lines(density(C3,bw=0.01,from=0.000001,to=20),col="pink")
abline(v=C,col="red")
legend("top",legend = c("prior","chain 1","chain 2","chain 3"),col=c("black","green","blue","pink"), lty = 1)

### D
C1 = out$BUGSoutput$sims.array[,1,'D']
C2 = out$BUGSoutput$sims.array[,2,'D']
C3 = out$BUGSoutput$sims.array[,3,'D']
curve(dgamma(x, shape = 0.01, rate = 0.01, log = FALSE),from=0.000001,to=1,ylim=c(0,20),ylab="Probability density",xlab = "D",main = "With KR data")
lines(density(C1,from=0.000001,to=1),col="green")
lines(density(C2,from=0.000001,to=1),col="blue")
lines(density(C3,from=0.000001,to=1),col="pink")
abline(v=D,col="red")
legend("topright",legend = c("prior","chain 1","chain 2","chain 3"),col=c("black","green","blue","pink"), lty = 1)

C1 = out2$BUGSoutput$sims.array[,1,'D']
C2 = out2$BUGSoutput$sims.array[,2,'D']
C3 = out2$BUGSoutput$sims.array[,3,'D']
curve(dgamma(x, shape = 0.01, rate = 0.01, log = FALSE),from=0.000001,to=1,ylim=c(0,20),ylab="Probability density",xlab = "D",main = "Without KR data")
lines(density(C1,bw=0.01,from=0.000001,to=1),col="green")
lines(density(C2,bw=0.01,from=0.000001,to=1),col="blue")
lines(density(C3,bw=0.01,from=0.000001,to=1),col="pink")
abline(v=D,col="red")
legend("topright",legend = c("prior","chain 1","chain 2","chain 3"),col=c("black","green","blue","pink"), lty = 1)
dev.off()

### Plot pair posterior densities
postsamples=cbind(out$BUGSoutput$sims.list$r_V,
out$BUGSoutput$sims.list$gamma,
out$BUGSoutput$sims.list$r_P,
out$BUGSoutput$sims.list$Q,
out$BUGSoutput$sims.list$C,
out$BUGSoutput$sims.list$D)
png(file="PairPosteriorPlot_withKRdata_LC.png", width = 1200, height = 1200,res=300)
pairs(postsamples,c("r","gamma","s","Q","C","D"))
dev.off()

postsamples2=cbind(out2$BUGSoutput$sims.list$r_V,
                  out2$BUGSoutput$sims.list$gamma,
                  out2$BUGSoutput$sims.list$r_P,
                  out2$BUGSoutput$sims.list$Q,
                  out2$BUGSoutput$sims.list$C,
                  out2$BUGSoutput$sims.list$D)
png(file="PairPosteriorPlot_withoutKRdata_LC.png", width = 1200, height = 1200, res=300)
pairs(postsamples2,c("r","gamma","s","Q","C","D"))
dev.off()
# 
# pdf(file="PairCorrelPosteriorPlot.pdf",width = 4,height = 8)
# par(mfrow=c(2,1))
# parcorplot(out,parms = c("r","gamma","s","Q","C","D"))
# parcorplot(out2,parms = c("r","gamma","s","Q","C","D"))
# dev.off()

# Output summary statistics
# jags.sum2<-out2$BUGSoutput$summary
# write.table(x=jags.sum2,file="JAGSsummary_PredPreyStoch_analysisBis_noFR.txt")
### Predicted values of the model fitted without the functional response
#logN2<-out2$BUGSoutput$mean$logN
#logP2<-out2$BUGSoutput$mean$logP
# par(mfrow=c(2,1))
# plot(2:(n.years),log(N)[2:(n.years)],type="o")
# lines(1:(n.years-1),logN1,type="o",col="red")
# lines(1:(n.years-1),logN2,type="o",col="blue")
# plot(2:(n.years),log(P)[2:(n.years)],type="o")
# lines(1:(n.years-1),logP1,type="o",col="red")
# lines(1:(n.years-1),logP2,type="o",col="blue")


#### Reproduce the FR curve without and without the correlations between parameters

## All curves for (C,D) in one Fig.
library(scales)
png('estimation_curves_FR_LC.png',res=300,width=4000,height=4000)
par(mfrow=c(2,2),cex=1.2) #NB cex =0.6 needed for width = 2000

Clist = out$BUGSoutput$sims.array[,,'C']
Dlist = out$BUGSoutput$sims.array[,,'D']
n = length(Clist)
ndens = 100
Nprey <- seq(1,20,length=ndens) #density index
FRstoch = matrix(NA,nrow = n, ncol = ndens)

for (i in 1:n){
  for (dens in 1:length(Nprey))
  {
    FRstoch[i,dens] = Clist[i]*Nprey[dens]/(Dlist[i]+Nprey[dens])
  }
  if (i == 1){plot(Nprey,FRstoch[i,],type='l',lwd=0.5,col=alpha('black',0.05),ylab='Functional response',ylim=c(0,20),xlim=c(1,10),xlab='N prey',main='With (C,D) correlations')
  }
  else {lines(Nprey,FRstoch[i,],type='l',lwd=0.5,col=alpha('black',0.05))}
  lines(Nprey,C*Nprey/(D+Nprey),col=alpha('red',1.0))
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
  if (i == 1){plot(Nprey,FRstoch[i,],type='l',lwd=0.5,col=alpha('black',0.05),ylab='Functional response',ylim=c(0,20),xlim=c(1,10),xlab='N prey',main='Without (C,D) correlations')
  }
  else {lines(Nprey,FRstoch[i,],type='l',lwd=0.5,col=alpha('black',0.05))}
  lines(Nprey,C*Nprey/(D+Nprey),col=alpha('red',1.0))
}


#### Same plot without the FR data

#Clist = out2$BUGSoutput$sims.list$C # works here but better to avoid because it destroys autocorrelation
#Dlist = out2$BUGSoutput$sims.list$D
Clist = out2$BUGSoutput$sims.array[,,'C']
Dlist = out2$BUGSoutput$sims.array[,,'D']
n = length(Clist)
ndens = 100
Nprey <- seq(1,20,length=ndens) #density index
FRstoch = matrix(NA,nrow = n, ncol = ndens)

for (i in 1:n){
  for (dens in 1:length(Nprey))
  {
    FRstoch[i,dens] = Clist[i]*Nprey[dens]/(Dlist[i]+Nprey[dens])
  }
  if (i == 1){plot(Nprey,FRstoch[i,],type='l',lwd=0.5,col=alpha('black',0.05),ylab='Functional response',ylim=c(0,20),xlim=c(1,10),xlab='N prey',main='With (C,D) correlations')
  }
  else {lines(Nprey,FRstoch[i,],type='l',lwd=0.5,col=alpha('black',0.05))}
  lines(Nprey,C*Nprey/(D+Nprey),col=alpha('red',1.0))
}

# what if there was no correlation between parameters? 

Clist = sample(Clist)
Dlist = sample(Dlist)
for (i in 1:n){
  for (dens in 1:length(Nprey))
  {
    FRstoch[i,dens] = Clist[i]*Nprey[dens]/(Dlist[i]+Nprey[dens])
  }
  if (i == 1){plot(Nprey,FRstoch[i,],type='l',lwd=0.5,col=alpha('black',0.05),ylab='Functional response',ylim=c(0,20),xlim=c(1,10),xlab='N prey',main='Without (C,D) correlations')
  }
  else {lines(Nprey,FRstoch[i,],type='l',lwd=0.5,col=alpha('black',0.05))}
  lines(Nprey,C*Nprey/(D+Nprey),col=alpha('red',1.0))
}

dev.off()

### Check why those values without the FR data are so low (different results for T=1000)
out2$BUGSoutput$mean$C*Nprey/(out2$BUGSoutput$mean$D+Nprey)


#### Reproduce the density-dependent curves without and without the correlations between parameters
# For the prey


png('estimation_curves_preydd_LC.png',res=300,width=4000,height=4000)
par(mfrow=c(2,2),cex=1.2) 

rlist = out$BUGSoutput$sims.array[,,'r_V']
gammalist = out$BUGSoutput$sims.array[,,'gamma']

n = length(rlist)
ndens = 100
Nprey <- seq(1,50,length=ndens) #density index
preyDD = matrix(NA,nrow = n, ncol = ndens)

for (i in 1:n){
  for (dens in 1:length(Nprey))
  {
    preyDD[i,dens] = exp(rlist[i])/(1+Nprey[dens]*gammalist[i])
  }
  if (i == 1){plot(Nprey,preyDD[i,],type='l',lwd=0.5,col=alpha('black',0.05),ylab='Prey growth rate',ylim=c(-1,6),xlim=c(1,50),xlab='N prey',main='With (r,gamma) correlations')
  }
  else {lines(Nprey,preyDD[i,],type='l',lwd=0.5,col=alpha('black',0.05))}
  lines(Nprey,exp(rmax_V)/(1+gamma*Nprey),col=alpha('red',1.0))
}

# what if there was no correlation between parameters? 
rlist = sample(rlist)
gammalist = sample(gammalist)
for (i in 1:n){
  for (dens in 1:length(Nprey))
  {
    preyDD[i,dens] = exp(rlist[i])/(1+Nprey[dens]*gammalist[i])
  }
  if (i == 1){plot(Nprey,preyDD[i,],type='l',lwd=0.5,col=alpha('black',0.05),ylab='Prey growth rate',ylim=c(-1,6),xlim=c(1,50),xlab='N prey',main='Without (r,gamma) correlations')
  }
  else {lines(Nprey,preyDD[i,],type='l',lwd=0.5,col=alpha('black',0.05))}
  lines(Nprey,exp(rmax_V)/(1+gamma*Nprey),col=alpha('red',1.0))
}


## without FR data
# For the prey
rlist = out2$BUGSoutput$sims.array[,,'r_V']
gammalist = out2$BUGSoutput$sims.array[,,'gamma']
n = length(rlist)
ndens = 100
Nprey <- seq(1,50,length=ndens) #density index
preyDD = matrix(NA,nrow = n, ncol = ndens)

for (i in 1:n){
  for (dens in 1:length(Nprey))
  {
    preyDD[i,dens] = exp(rlist[i])/(1+Nprey[dens]*gammalist[i])
  }
  if (i == 1){plot(Nprey,preyDD[i,],type='l',lwd=0.5,col=alpha('black',0.05),ylab='Prey growth rate',ylim=c(-1,6),xlim=c(1,50),xlab='N prey',main='With (r,gamma) correlations')
  }
  else {lines(Nprey,preyDD[i,],type='l',lwd=0.5,col=alpha('black',0.05))}
  lines(Nprey,exp(rmax_V)/(1+gamma*Nprey),col=alpha('red',1.0))
}

# what if there was no correlation between parameters? 
rlist = sample(rlist)
gammalist = sample(gammalist)
for (i in 1:n){
  for (dens in 1:length(Nprey))
  {
    preyDD[i,dens] = exp(rlist[i])/(1+Nprey[dens]*gammalist[i])
  }
  if (i == 1){plot(Nprey,preyDD[i,],type='l',lwd=0.5,col=alpha('black',0.05),ylab='Prey growth rate',ylim=c(-1,6),xlim=c(1,50),xlab='N prey',main='Without (r,gamma) correlations')
  }
  else {lines(Nprey,preyDD[i,],type='l',lwd=0.5,col=alpha('black',0.05))}
  lines(Nprey,exp(rmax_V)/(1+gamma*Nprey),col=alpha('red',1.0))
}

dev.off()

