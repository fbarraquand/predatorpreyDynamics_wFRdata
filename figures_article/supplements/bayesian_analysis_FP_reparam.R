
########## Bayesian analysis of the predator-prey system with and without kill rate data ####################

########## Updated code with the reparameterization 27/03/2019

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

# Loading data
file_path1 = "../../simulations/small_noise_on_FR/parameter_sets/perturbed_fixed_point/T=100/"
krep=1 # sample chosen
data1 = read.csv(paste(file_path1,"predatorPrey_withGaussianFR",krep,".csv",sep=""))
names(data1)=c("Time","n","p","KR")
head(data1)

# Bundle data
jags.data <- list(T=n.years,logN=data1$n,logP=data1$p,FR=data1$KR)


sink("predpreymod.txt")
cat("
    model {
    
    # Priors and constraints
    logN[1] ~ dnorm(0,0.01) # Prior on initial pop size on the log scale
    logP[1] ~ dnorm(0,0.01) # Prior on initial pop size on the log scale
    
    # Priors prey population dynamics
    r_V ~ dnorm(1,1) #dnorm(1,0.001) # below the truth, rather flat prior
    K_V ~ dunif(0.2,20)
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
    MV <- K_V / max(exp(r_V)-1,0.01)
    #Another option is to restrict r_V and r_P to positive values 
    QP <- Q * max(exp(r_P)-1,0.1) 
    
    
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
  list(sigma_V=runif(1,0.1,2), sigma_P=runif(1,0.1,2), r_V=runif(1,1,2),r_P=runif(1,0.3,0.8), K_V=runif(1,0.2,15), Q=runif(1,4,10),tau_FR=runif(1,1,10),a=runif(1,1,10),h=runif(1,0.1,1))}

# Parameters monitored
parameters<-c("r_V","K_V","r_P","Q","sigma2_V","sigma2_P","tau_FR","a","h")

# MCMC settings
nc <- 3 #number of chains
nb <- 24000 # “burn in”
#ni <- 14000# “number of iterations” # that's for a symmetric distrib...
ni<-44000
nt <- 10 # “thinning”

# run model
out <- jags(jags.data, inits, parameters, "predpreymod.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, working.directory = getwd())
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
    K_V ~ dunif(0.2,20)
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
  list(sigma_V=runif(1,0.1,2), sigma_P=runif(1,0.1,2), r_V=runif(1,0.1,2),r_P=runif(1,0.1,2), K_V=runif(1,0.2,15), Q=runif(1,1,10),a=runif(1,1,10),h=runif(1,0.1,1))}

# Parameters monitored
parameters<-c("r_V","K_V","r_P","Q","sigma2_V","sigma2_P","a","h")


# MCMC settings
nc <- 3 #number of chains
nb <- 24000 # “burn in”
#ni <- 14000# “number of iterations” # that's for a symmetric distrib...
ni<-44000
nt <- 10 # “thinning”

# run model
out2 <- jags(jags.data, inits, parameters, "predpreymod_without_sepFR.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, working.directory = getwd())
print(out2, dig = 2)

#save.image('/media/frederic/DATA/Simuls_wOlivier/predatorpreyFRdata/simTS100.RData')

### plot densities
library(mcmcplots)
denplot(out,c("r_V","gamma","r_P","Q","sigma2_V","sigma2_P","a","h"))
denplot(out2,c("r_V","gamma","r_P","Q","sigma2_V","sigma2_P","a","h"))

### Trace plots
traplot(out,c("r_V","gamma","r_P","Q","sigma2_V","sigma2_P","a","h"))
traplot(out2,c("r_V","gamma","r_P","Q","sigma2_V","sigma2_P","a","h"))

## Prior-posterior overlap

par(mfrow=c(2,1),lwd=2)
a1 = out$BUGSoutput$sims.array[,1,'a']
a2 = out$BUGSoutput$sims.array[,2,'a']
a3 = out$BUGSoutput$sims.array[,3,'a']
curve(dgamma(x, shape = 0.01, rate = 0.01, log = FALSE),from=0.000001,to=5,ylim=c(0,3),ylab="Probability density",xlab = "a",main = "With KR data")
lines(density(a1,from=0.000001,to=5),col="green")
lines(density(a2,from=0.000001,to=5),col="blue")
lines(density(a3,from=0.000001,to=5),col="pink")
abline(v=a,col="red")
legend("topright",legend = c("prior","chain 1","chain 2","chain 3"),col=c("black","green","blue","pink"), lty = 1)

a1 = out2$BUGSoutput$sims.array[,1,'a']
a2 = out2$BUGSoutput$sims.array[,2,'a']
a3 = out2$BUGSoutput$sims.array[,3,'a']
curve(dgamma(x, shape = 0.01, rate = 0.01, log = FALSE),from=0.000001,to=5,ylim=c(0,3),ylab="Probability density",xlab = "a",main = "Without KR data")
lines(density(a1,bw=0.01,from=0.000001,to=5),col="green")
lines(density(a2,bw=0.01,from=0.000001,to=5),col="blue")
lines(density(a3,bw=0.01,from=0.000001,to=5),col="pink")
abline(v=C,col="red")
legend("topright",legend = c("prior","chain 1","chain 2","chain 3"),col=c("black","green","blue","pink"), lty = 1)


### With both a and h

pdf(file="PPO_FP_reparam.pdf",width=12,height=8)
par(mfrow=c(2,2),lwd=2,cex=1.2)

a1 = out$BUGSoutput$sims.array[,1,'a']
a2 = out$BUGSoutput$sims.array[,2,'a']
a3 = out$BUGSoutput$sims.array[,3,'a']
curve(dgamma(x, shape = 0.01, rate = 0.01, log = FALSE),from=0.000001,to=5,ylim=c(0,3),ylab="Probability density",xlab = "a",main = "With KR data")
lines(density(a1,from=0.000001,to=5),col="green")
lines(density(a2,from=0.000001,to=5),col="blue")
lines(density(a3,from=0.000001,to=5),col="pink")
abline(v=a,col="red")
legend("topright",legend = c("prior","chain 1","chain 2","chain 3"),col=c("black","green","blue","pink"), lty = 1)

a1 = out2$BUGSoutput$sims.array[,1,'a']
a2 = out2$BUGSoutput$sims.array[,2,'a']
a3 = out2$BUGSoutput$sims.array[,3,'a']
curve(dgamma(x, shape = 0.01, rate = 0.01, log = FALSE),from=0.000001,to=5,ylim=c(0,3),ylab="Probability density",xlab = "a",main = "Without KR data")
lines(density(a1,bw=0.01,from=0.000001,to=5),col="green")
lines(density(a2,bw=0.01,from=0.000001,to=5),col="blue")
lines(density(a3,bw=0.01,from=0.000001,to=5),col="pink")
abline(v=a,col="red")
legend("topright",legend = c("prior","chain 1","chain 2","chain 3"),col=c("black","green","blue","pink"), lty = 1)

### h
C1 = out$BUGSoutput$sims.array[,1,'h']
C2 = out$BUGSoutput$sims.array[,2,'h']
C3 = out$BUGSoutput$sims.array[,3,'h']
curve(dgamma(x, shape = 0.01, rate = 0.01, log = FALSE),from=0.000001,to=2,ylim=c(0,15),ylab="Probability density",xlab = "h",main = "With KR data")
lines(density(C1,from=0.000001,to=2),col="green")
lines(density(C2,from=0.000001,to=2),col="blue")
lines(density(C3,from=0.000001,to=2),col="pink")
abline(v=h,col="red")
legend("topright",legend = c("prior","chain 1","chain 2","chain 3"),col=c("black","green","blue","pink"), lty = 1)

C1 = out2$BUGSoutput$sims.array[,1,'h']
C2 = out2$BUGSoutput$sims.array[,2,'h']
C3 = out2$BUGSoutput$sims.array[,3,'h']
curve(dgamma(x, shape = 0.01, rate = 0.01, log = FALSE),from=0.000001,to=2,ylim=c(0,15),ylab="Probability density",xlab = "h",main = "Without KR data")
lines(density(C1,bw=0.01,from=0.000001,to=2),col="green")
lines(density(C2,bw=0.01,from=0.000001,to=2),col="blue")
lines(density(C3,bw=0.01,from=0.000001,to=2),col="pink")
abline(v=h,col="red")
legend("topright",legend = c("prior","chain 1","chain 2","chain 3"),col=c("black","green","blue","pink"), lty = 1)
dev.off()

### Plot pair posterior densities
postsamples=cbind(out$BUGSoutput$sims.list$r_V,
out$BUGSoutput$sims.list$gamma,
out$BUGSoutput$sims.list$r_P,
out$BUGSoutput$sims.list$Q,
out$BUGSoutput$sims.list$C,
out$BUGSoutput$sims.list$D)
png(file="PairPosteriorPlot_withKRdata_FP.png", width = 1200, height = 1200,res=300)
pairs(postsamples,c("r","K","s","Q","a","h"))
dev.off()

postsamples2=cbind(out2$BUGSoutput$sims.list$r_V,
                  out2$BUGSoutput$sims.list$gamma,
                  out2$BUGSoutput$sims.list$r_P,
                  out2$BUGSoutput$sims.list$Q,
                  out2$BUGSoutput$sims.list$C,
                  out2$BUGSoutput$sims.list$D)
png(file="PairPosteriorPlot_withoutKRdata_FP.png", width = 1200, height = 1200, res=300)
pairs(postsamples2,c("r","gamma","s","Q","C","D"))
dev.off()


#### Reproduce the FR curve without and without the correlations between parameters

## All curves for (C,D) in one Fig.
library(scales)
png('estimation_curves_FR_FP_reparam.png',res=300,width=2000,height=2000)
par(mfrow=c(2,2),cex=0.8) #NB cex =0.6 needed for width = 2000, 1.2 for width =4000

alist = out$BUGSoutput$sims.array[,,'a']
hlist = out$BUGSoutput$sims.array[,,'h']
n = length(alist)
ndens = 100
Nprey <- seq(1,20,length=ndens) #density index
FRstoch = matrix(NA,nrow = n, ncol = ndens)

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


#### Same plot without the FR data
alist = out2$BUGSoutput$sims.array[,,'a']
hlist = out2$BUGSoutput$sims.array[,,'h']
n = length(alist)
ndens = 100
Nprey <- seq(1,20,length=ndens) #density index
FRstoch = matrix(NA,nrow = n, ncol = ndens)

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

#### Reproduce the density-dependent curves without and without the correlations between parameters
# For the prey

png('estimation_curves_preydd_FP_reparam.png',res=300,width=2000,height=2000)
par(mfrow=c(2,2),cex=0.8) 

rlist = out$BUGSoutput$sims.array[,,'r_V']
Klist = out$BUGSoutput$sims.array[,,'K_V']

n = length(rlist)
ndens = 100
Nprey <- seq(1,50,length=ndens) #density index
preyDD = matrix(NA,nrow = n, ncol = ndens)


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


## without FR data
# For the prey
rlist = out2$BUGSoutput$sims.array[,,'r_V']
Klist = out2$BUGSoutput$sims.array[,,'K_V']
n = length(rlist)
ndens = 100
Nprey <- seq(1,50,length=ndens) #density index
preyDD = matrix(NA,nrow = n, ncol = ndens)

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

