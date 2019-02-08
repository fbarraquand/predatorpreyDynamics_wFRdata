### FB 13/11/2017 - predator-prey model with noisy functional response data
### From much earlier code
### 15/03/2018 Added comparison to deterministic version of the code
### 02/02/2019 Removed comparison, loop over 100 repeats for T=100

rm(list=ls())
graphics.off()

### Model used (use stored simulated data later on, just tryouts for now)
n.repeats<-100
n.years<-100  	# Number of years //
N1<-1			# Initial pop size
P1<-0.1   # O.1 before but difficult to see if there are transient oscillations
K<-1			# threshold dd 
beta<-1			# density-dependence exponent
rmax_V<-2			# Max AVERAGE growth rate (thus not a true max...)
rmax_P<-0.5
sigma2.proc<-0.05		# worked well with 0.005 too

### FR and predator parameters
C<-15
D<-0.25
Q<-10

### Simulation of data
#set.seed(42) 
set.seed(41)
#set.seed(40)

for (krep in 1:n.repeats){

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
  FR[t+1]<-(C*N[t+1]/(D+N[t+1])) + FRnoise[t+1] #must be after for updating. 
}

### Produce dataset to fit
data = cbind(log(N),log(P),FR) 

## Writing to data file
if (sum(sapply(data,is.finite))<300){
  message(paste0("Pb with rep:",krep))
 write.csv(data,file=paste("predatorPrey_withGaussianFR",krep,".csv",sep=""))
}
else
{
  
  ## Plotting time series of abundances and FR
  par(mfrow=c(2,2))
  plot(1:n.years,N,type="b")
  plot(1:n.years,P,type="b")
  spectrum(log(N),method="ar")
  spectrum(log(P),method="ar")
  
  
  par(mfrow=c(2,2))
  plot(N,P)
  plot(N,FR)
  plot(P/N)
  plot(log(N),log(P))
  
write.csv(format(data,digits=3),file=paste("predatorPrey_withGaussianFR",krep,".csv",sep=""))
}

}
