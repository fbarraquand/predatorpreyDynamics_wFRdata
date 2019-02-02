#### 15/03/2018 -- FBarraquand and OGimenez
#### Computation of the Hessian for the discrete-time Leslie-May model
### FB 16/01/2019 -- edited so that functional response timing match math model (no delays)
### FB 02/02/2019 -- Here we compute first the Hessian for 10 x T=1000 case. Then we average to see the FIM 

rm(list=ls())
graphics.off()
library("numDeriv")
source('LikelihoodFunctions.R')


### Parameters
n.repeats<-10
n.years<-1000  	# Long time series to look at structural identifiability
N1<-1			# Initial pop size
P1<-0.1
K<-1			# threshold dd 
beta<-1			# density-dependence exponent
rmax_V<-2			# Max AVERAGE growth rate (thus not a true max...)
rmax_P<-0.5
sigma2.proc<-0.05		# worked well with 0.005
### FR and predator parameters
C<-2.5
D<-1
Q<-10

### Simulation of data and computation of Hessian
#set.seed(42) 
set.seed(41)
#set.seed(40)

### Initialize structures
H1 = array(0,c(9,9,n.repeats))
H2 = array(0,c(8,8,n.repeats))


for (krep in 1:n.repeats){

  print(krep)
  
y<-N<-P<-FR<-numeric(n.years)
N[1]<-N1
P[1]<-P1
FR[1]<-C*N[1]/(D+N[1])

rV<-rnorm(n.years-1,rmax_V,sqrt(sigma2.proc))
rP<-rnorm(n.years-1,rmax_P,sqrt(sigma2.proc))
FRnoise<-rnorm(n.years,0,sqrt(sigma2.proc))

for (t in 1:(n.years-1)){
  
  N[t+1]<-N[t]*(exp(rV[t])/(1+(N[t]/K)^beta))*exp(-FR[t]*P[t]/N[t])
  P[t+1]<-P[t]*exp(rP[t])/(1+P[t]*Q/N[t])
  FR[t+1]<-(C*N[t+1]/(D+N[t+1])) + FRnoise[t+1] #updated after so that timing matches
}
## Plotting time series of abundances and FR
par(mfrow=c(2,2))
plot(1:n.years,N,type="b")
plot(1:n.years,P,type="b")
#curve(dbeta(x,a,b),from=0, to=1)
plot(N,FR)

### Produce dataset to fit
data = cbind(log(N),log(P),FR) 

#################### Computation of Hessian at the true parameter values ######################

theta_true  = c(rmax_V,1/K,sqrt(0.05),rmax_P,Q,sqrt(0.05),C,D,sqrt(0.05))

LL = function(theta){
  return(logLik(theta,data))
}
LL(theta_true)
H1[,,krep] = hessian(LL,theta_true)

### Removing last parameter for test with LL_FRwoutNoise
theta_true  = c(rmax_V,1/K,sqrt(0.05),rmax_P,Q,sqrt(0.05),C,D)

LL_FRwoutNoise = function(theta){
  return(logLik_FRwoutNoise(theta,data))
}

hessian_thetaTrue_FRwoutNoise=
H2[,,krep] = hessian(LL_FRwoutNoise,theta_true)

}


### Using Viallefond et al.'s criteria
lambda1 = max(eigen(H1[,,1])$values)
q = 9
eps = 10^(-9)
eigen(H1[,,1])$values>q*lambda1*eps

lambda1 = max(eigen(H2[,,1])$values)
q = 8
eigen(H2[,,1])$values>q*lambda1*eps

lambda1 = max(eigen(H1[,,2])$values)
q = 9
eigen(H1[,,2])$values>q*lambda1*eps

lambda1 = max(eigen(H2[,,2])$values)
q = 8
eigen(H2[,,2])$values>q*lambda1*eps
##

## Average array over last dimension?
H1mean = apply(H1,c(1, 2),mean)
H2mean = apply(H2,c(1, 2),mean)

### Add to array
require(abind)
H1 = abind(H1,H1mean)
H2 = abind(H2,H2mean)

### Do a number of dotchart
par(mfrow=c(4,3))
for (krep in 1:(n.repeats+1)){
  eigentable = cbind((eigen(H1[,,krep])$values),c((eigen(H2[,,krep])$values),NA))
  dotchart(eigentable,gcolor=c("blue","red"),pch=19,xlab="eigenvalue",ylab="Rank eigenvalue",xlim=c(-100,50000))
}


eps=1
par(mfrow=c(4,3))
for (krep in 1:(n.repeats+1)){
eigentable = cbind(log10(eigen(H1[,,krep])$values+eps),c(log10(eigen(H2[,,krep])$values+eps),NA))
dotchart(eigentable,gcolor=c("blue","red"),pch=19,xlab="log10(eigenvalue+1)",ylab="Rank eigenvalue")
}

### Looking at the FIM 
eigen(H1mean)$values
lambda1 = max(eigen(H1mean)$values)
q = 9
eps = 10^(-9)
eigen(H1mean)$values>q*lambda1*eps

eigen(H2mean)$values
lambda1 = max(eigen(H2mean)$values)
q = 8
eigen(H2mean)$values>q*lambda1*eps
### Looks like the last eigenvalue is one order of magnitude below, but what can we deduce from that? 
### We still have 10^-5 difference between the largest and lowest eigenvalue instead of 10^(-4)