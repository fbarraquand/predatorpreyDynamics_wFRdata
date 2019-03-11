### F. Barraquand 21/04/2015 - Code for analyzing noisy predator-prey system (incl. noise on the functional response)
### Changed to Rosenzweig-MacArthur version 14/03/2018
### Updated 01/06/2017. 
### Case with "small noise" on the functional response
##################################################################################################################################

rm(list=ls())
graphics.off()

library("R2jags")      # Load R2jags package

### Parameters for simulation of Hassell model

n.years<-200  	# Number of years - 25 first, perhaps use 50 or 100 / worked very well with almost no process error on the real scale
#N1<-1			# Initial pop size
#P1<-0.1 ### Too long transients

N1<-4
P1<-0.7

K<-1			# threshold dd 
beta<-1			# density-dependence exponent
rmax_V<-2			# Max AVERAGE growth rate (thus not a true max...)
rmax_P<-(-0.20) #
sigma2.proc<-0.05		# worked well with 0.005
# Process sigma on the log-scale, use the Peretti et al. value. 0.005


### FR and predator parameters
C<-2.5
D<-1.0 #0.6 = quasi-cycles
epsilon<-0.1

### Simulation of data
set.seed(42) 
#set.seed(41)
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
  P[t+1]<-P[t]*exp(rP[t]+epsilon*FR[t])
  FR[t+1]<-(C*N[t+1]/(D+N[t+1])) + FRnoise[t+1]
}
##### Other question: should I use an interval a little less than one to stabilize this?

## Plotting time series of abundances and FR
par(mfrow=c(2,2))
plot(1:n.years,N,type="b")
plot(1:n.years,P,type="b")
#curve(dbeta(x,a,b),from=0, to=1)
plot(N,FR)

    
#seq(0,1,0.01)
# Bundle data

after_transients=101:200
jags.data <- list(T=length(after_transients),logN=log(N[after_transients]),logP=log(P[after_transients]),FR=FR[after_transients])

sink("predprey.rma.txt")
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
    epsilon ~ dunif(0.01,0.25)
    r_P ~ dnorm(0,0.1)
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
    logNupdate[t] <- logN[t] + r_V -log(1+N[t]/K_V) - FR[t]*exp(logP[t])/N[t]
    logN[t+1] ~ dnorm(logNupdate[t],tau_V)
    N[t]<-exp(logN[t])

    # for some reason, log(1+(exp(r_V)-1)*N[t]/K_V) was not working
    logP[t+1]~ dnorm(logPupdate[t],tau_P)
    logPupdate[t] <- logP[t] + r_P + epsilon*FR[t]
    
    }
    
    }
    ",fill=TRUE)
sink()


# Initial values
inits <- function () {
  list(sigma_V=runif(1,0.1,2), sigma_P=runif(1,0.1,2), r_V=runif(1,0.1,2),r_P=runif(1,-0.3,0), K_V=runif(1,0.2,8), epsilon=runif(1,0.01,0.25),tau_FR=runif(1,1,10),C=runif(1,10,100),D=runif(1,0.01,0.1))}

# Parameters monitored
#parameters<-c("r_V","K_V","r_P","epsilon","sigma2_V","sigma2_P","C","D","logN","logP","FR")
parameters<-c("r_V","K_V","r_P","epsilon","sigma2_V","sigma2_P","tau_FR","C","D")


# MCMC settings
nc <- 3 #number of chains
nb <- 24000 # “burn in”
#ni <- 14000# “number of iterations” # that's for a symmetric distrib...
ni<-44000
nt <- 10 # “thinning”

# run model
out <- jags(jags.data, inits, parameters, "predprey.rma.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, working.directory = getwd())
print(out, dig = 2)

# Good title for a future paper (if I compute more quantities of Trophic Strength from this...) would be
# An integrated assessment of trophic interaction strength. 

# Output summary statistics
jags.sum<-out$BUGSoutput$summary
#write.table(x=jags.sum,file="JAGSsummary_PredPreyStoch_analysisBis.txt")
# MCMC Output
#pdf("Output_MCMC__PredPreyStoch_analysisBis.pdf")
out.mcmc<-as.mcmc(out)
plot(out.mcmc) 
#dev.off()

CEb<-out$BUGSoutput$mean$C
DEb<-out$BUGSoutput$mean$D

logN1<-out$BUGSoutput$mean$logN
logP1<-out$BUGSoutput$mean$logP

par(mfrow=c(2,2))
plot(1:(n.years-1),logN1,type="o")
plot(1:(n.years-1),logP1,type="o")

### Fit functional response

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
    epsilon ~ dunif(0.01,0.25)
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
    logPupdate[t] <- logP[t] + r_P + epsilon*C*exp(logP[t])/(D+N[t])
    
    }
    
    }
    ",fill=TRUE)
sink()


# Initial values
inits <- function () {
  list(sigma_V=runif(1,0.1,2), sigma_P=runif(1,0.1,2), r_V=runif(1,0.1,2),r_P=runif(1,-0.3,0), K_V=runif(1,0.2,8), epsilon=runif(1,0.01,0.25),C=runif(1,10,100),D=runif(1,0.01,0.1))}

# Parameters monitored
#parameters<-c("r_V","K_V","r_P","epsilon","sigma2_V","sigma2_P","C","D","logN","logP","FR")
parameters<-c("r_V","K_V","r_P","epsilon","sigma2_V","sigma2_P","C","D")


# MCMC settings
nc <- 3 #number of chains
nb <- 24000 # “burn in”
#ni <- 14000# “number of iterations” # that's for a symmetric distrib...
ni<-44000
nt <- 10 # “thinning”

# run model
out2 <- jags(jags.data, inits, parameters, "predprey_without_sepFR.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, working.directory = getwd())
print(out2, dig = 2)

# http://jeromyanglim.tumblr.com/post/37362047458/how-to-get-dic-in-jags
print(out,dig=2) # to compare deviance, DIC - check also parameter values.
# 01/06/2017 // Olivier says DIC is crap in this (and other?) context, avoid this...  
plot(as.mcmc(out2)) 
plot(as.mcmc(out)) 

### Predicted values of the model fitted without the functional response
logN2<-out2$BUGSoutput$mean$logN
logP2<-out2$BUGSoutput$mean$logP

par(mfrow=c(2,1))
plot(2:(n.years),log(N)[2:(n.years)],type="o")
lines(1:(n.years-1),logN1,type="o",col="red")
lines(1:(n.years-1),logN2,type="o",col="blue")
plot(2:(n.years),log(P)[2:(n.years)],type="o")
lines(1:(n.years-1),logP1,type="o",col="red")
lines(1:(n.years-1),logP2,type="o",col="blue")

## NB Perhaps it is much easier to identify all the parameters when the model has close to zero environmental noise. 
## Using sigma2 = 0.005 At least C is identified approximately correctly, though not D. 
## Still much better fit with the FR it seems though. 


### plot densities
library(mcmcplots)
denplot(out,c("r_V","K_V","r_P","epsilon","sigma2_V","sigma2_P","C","D"))
denplot(out2,c("r_V","K_V","r_P","epsilon","sigma2_V","sigma2_P","C","D"))

### Trace plots
traplot(out,c("r_V","K_V","r_P","epsilon","sigma2_V","sigma2_P","C","D"))
traplot(out2,c("r_V","K_V","r_P","epsilon","sigma2_V","sigma2_P","C","D"))


### Plot pair posterior densities
postsamples=cbind(out$BUGSoutput$sims.list$r_V,
                  out$BUGSoutput$sims.list$K_V,
                  out$BUGSoutput$sims.list$r_P,
                  out$BUGSoutput$sims.list$epsilon,
                  out$BUGSoutput$sims.list$C,
                  out$BUGSoutput$sims.list$D)
png(file="PairPosteriorPlot_withFRdata.png", width = 1200, height = 1200,res=300)
pairs(postsamples,c("r_V","K_V","r_P","epsilon","C","D"))
dev.off()

postsamples2=cbind(out2$BUGSoutput$sims.list$r_V,
                   out2$BUGSoutput$sims.list$K_V,
                   out2$BUGSoutput$sims.list$r_P,
                   out2$BUGSoutput$sims.list$epsilon,
                   out2$BUGSoutput$sims.list$C,
                   out2$BUGSoutput$sims.list$D)
png(file="PairPosteriorPlot_withoutFRdata.png", width = 1200, height = 1200, res=300)
pairs(postsamples2,c("r_V","K_V","r_P","epsilon","C","D"))
dev.off()

pdf(file="PairCorrelPosteriorPlot.pdf",width = 4,height = 8)
par(mfrow=c(2,1))
parcorplot(out,parms = c("r_V","K_V","r_P","epsilon","sigma2_V","sigma2_P","C","D"))
parcorplot(out2,parms = c("r_V","K_V","r_P","epsilon","sigma2_V","sigma2_P","C","D"))
dev.off()

#### Prior posterior overlap
### With both C and D
pdf(file="PPO_FP_aftertransients.pdf",width=12,height=8)
par(mfrow=c(2,2),lwd=2,cex=1.2)

C1 = out$BUGSoutput$sims.array[,1,'C']
C2 = out$BUGSoutput$sims.array[,2,'C']
C3 = out$BUGSoutput$sims.array[,3,'C']
curve(dgamma(x, shape = 0.01, rate = 0.01, log = FALSE),from=0.000001,to=5,ylim=c(0,3),ylab="Probability density",xlab = "C",main = "With KR data")
lines(density(C1,from=0.000001,to=5),col="green")
lines(density(C2,from=0.000001,to=5),col="blue")
lines(density(C3,from=0.000001,to=5),col="pink")
abline(v=C,col="red")
legend("topright",legend = c("prior","chain 1","chain 2","chain 3"),col=c("black","green","blue","pink"), lty = 1)

C1 = out2$BUGSoutput$sims.array[,1,'C']
C2 = out2$BUGSoutput$sims.array[,2,'C']
C3 = out2$BUGSoutput$sims.array[,3,'C']
curve(dgamma(x, shape = 0.01, rate = 0.01, log = FALSE),from=0.000001,to=5,ylim=c(0,3),ylab="Probability density",xlab = "C",main = "Without KR data")
lines(density(C1,bw=0.01,from=0.000001,to=5),col="green")
lines(density(C2,bw=0.01,from=0.000001,to=5),col="blue")
lines(density(C3,bw=0.01,from=0.000001,to=5),col="pink")
abline(v=C,col="red")
legend("topright",legend = c("prior","chain 1","chain 2","chain 3"),col=c("black","green","blue","pink"), lty = 1)

### D
C1 = out$BUGSoutput$sims.array[,1,'D']
C2 = out$BUGSoutput$sims.array[,2,'D']
C3 = out$BUGSoutput$sims.array[,3,'D']
curve(dgamma(x, shape = 0.01, rate = 0.01, log = FALSE),from=0.000001,to=5,ylim=c(0,2),ylab="Probability density",xlab = "D",main = "With KR data")
lines(density(C1,from=0.000001,to=5),col="green")
lines(density(C2,from=0.000001,to=5),col="blue")
lines(density(C3,from=0.000001,to=5),col="pink")
abline(v=D,col="red")
legend("topright",legend = c("prior","chain 1","chain 2","chain 3"),col=c("black","green","blue","pink"), lty = 1)

C1 = out2$BUGSoutput$sims.array[,1,'D']
C2 = out2$BUGSoutput$sims.array[,2,'D']
C3 = out2$BUGSoutput$sims.array[,3,'D']
curve(dgamma(x, shape = 0.01, rate = 0.01, log = FALSE),from=0.000001,to=5,ylim=c(0,2),ylab="Probability density",xlab = "D",main = "Without KR data")
lines(density(C1,bw=0.01,from=0.000001,to=5),col="green")
lines(density(C2,bw=0.01,from=0.000001,to=5),col="blue")
lines(density(C3,bw=0.01,from=0.000001,to=5),col="pink")
abline(v=D,col="red")
legend("topright",legend = c("prior","chain 1","chain 2","chain 3"),col=c("black","green","blue","pink"), lty = 1)
dev.off()




