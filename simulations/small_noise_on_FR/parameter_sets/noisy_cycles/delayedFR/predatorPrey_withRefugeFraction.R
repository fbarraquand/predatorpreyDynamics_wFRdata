#################################################################################
### FB 15/03/2018 - finding parameter set with noisy limit cycle (i.e. limit cycle present when sigma ->0)
### FB 13/11/2017 - predator-prey model with noisy functional response data
### From much earlier code -- see Bayesian folder

rm(list=ls())
graphics.off()

################# Simulating the model without any noise ######################################

############# New parameter set --- stable limit cycles without stochasticity #############
############# Formally an invariant loop ##################################################

### Model used (use stored simulated data later on, just tryouts for now)
n.years<-300  	# Number of years // large so that we have no bias due to data length
N1<-1			# Initial pop size
P1<-0.1
K<-1			# threshold dd 
beta<-1			# density-dependence exponent
rmax_V<-2			# Max AVERAGE growth rate (thus not a true max...)
rmax_P<-0.5
sigma2.proc<-0.0		# worked well with 0.005
# Process sigma on the log-scale, use the Peretti et al. value. 0.005


### FR and predator parameters
C<-15 #instead 2.5
D<-0.25 #instead 1  (0.1 provoke extinction with C=15 and Q=10)
#NB 0.25 can still provoke extinction, 0.5 does not provide noise-free invariant loop
Q<-10 #possible 100 or 1000 otherwise 

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
FRnoise<-rnorm(n.years,0,sqrt(sigma2.proc))

for (t in 1:(n.years-1)){
  
  FR[t+1]<-(C*N[t]/(D+N[t])) + FRnoise[t+1]
  N[t+1]<-0.05+N[t]*(exp(rV[t])/(1+(N[t]/K)^beta))*(0.01+(1-0.01)*exp(-FR[t]*P[t]/N[t]))
  P[t+1]<-P[t]*exp(rP[t])/(1+P[t]*Q/N[t])
}
## Plotting time series of abundances and FR
par(mfrow=c(2,2))
plot(1:n.years,N,type="b")
plot(1:n.years,P,type="b")

#curve(dbeta(x,a,b),from=0, to=1)
#par(mfrow=c(2,2))
spectrum(log(N),method="ar")
spectrum(log(P),method="ar")

plot(N,FR)
plot(N,P)

### Produce dataset to fit
data = cbind(log(N),log(P),FR) 

###################### Noisy limit cycle #############################################

rm(list=ls())
graphics.off()

### Model used (use stored simulated data later on, just tryouts for now)
n.years<-300  	# Number of years // large so that we have no bias due to data length
N1<-1			# Initial pop size
P1<-0.1
K<-1			# threshold dd 
beta<-1			# density-dependence exponent
rmax_V<-2			# Max AVERAGE growth rate (thus not a true max...)
rmax_P<-0.5
sigma2.proc<-0.005		# 0.05 usual but worked well with 0.005
# Process sigma on the log-scale, use the Peretti et al. value. 0.005


### FR and predator parameters
C<-15 #2.5
D<-0.25 #0.5 #1
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
FRnoise<-rnorm(n.years,0,sqrt(sigma2.proc))

for (t in 1:(n.years-1)){
 
  FR[t+1]<-max((C*N[t]/(D+N[t])) + FRnoise[t+1],0.0001)
  N[t+1]<-N[t]*(exp(rV[t])/(1+(N[t]/K)^beta))*(0.01+(1-0.01)*exp(-FR[t]*P[t]/N[t]))
  P[t+1]<-P[t]*exp(rP[t])/(1+P[t]*Q/N[t])
}
## Plotting time series of abundances and FR
par(mfrow=c(2,2))
plot(1:n.years,N,type="b")
plot(1:n.years,P,type="b")
#curve(dbeta(x,a,b),from=0, to=1)

#par(mfrow=c(2,2))
spectrum(log(N),method="ar")
spectrum(log(P),method="ar")

par(mfrow=c(1,2))
plot(N,P)
plot(N,FR)
