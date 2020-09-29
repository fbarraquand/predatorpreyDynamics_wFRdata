### Optimizing the likelihood in a frequentist setting
### Doing this for the 100 datasets so as to investigate frequentist estimators
### Random initial condition for BFGS
### FB 29/09/2020

rm(list=ls())
graphics.off()

# taking simuls from 
file_path1 = "../simulations/small_noise_on_FR/parameter_sets/perturbed_fixed_point/T=100/"
file_path2 = "../simulations/small_noise_on_FR/parameter_sets/noisy_cycles/T=100/"

### Datasets of the form
### data = cbind(log(N),log(P),FR) 

#krep=1 # sample chosen (if one)
C_est = D_est = C_est_woutFR = D_est_woutFR = rep(NA,100) #Initializing estimators of C and D

for (krep in 1:100){

# loading dataset, perturbed FP
data1 = read.csv(paste(file_path1,"predatorPrey_withGaussianFR",krep,".csv",sep=""))
names(data1)=c("Time","n","p","KR")
data1[,1]<-NULL
head(data1)

# switch to noisy LC dataset
data2 = read.csv(paste(file_path2,"predatorPrey_withGaussianFR",krep,".csv",sep=""))
names(data2)=c("Time","n","p","KR")
data2[,1]<-NULL
head(data2)

## Sourcing the functions to compute the likelihood
source('../frequentist/LikelihoodFunctions.R')
source('../figures_article/figlabel.R')
cex_labels = 1.3

### Parameters of the model - fixed point parameter set
K<-1			# threshold dd 
rmax_V<-2			# Max AVERAGE growth rate (thus not a true max...)
rmax_P<-0.5
sigma2.proc<-0.05	
### FR and predator parameters
C<-2.5
D<-1
Q<-10

theta_true1  = c(rmax_V,1/K,sqrt(0.05),rmax_P,Q,sqrt(0.05),C,D,sqrt(0.05))
theta_true2  = c(rmax_V,1/K,sqrt(0.05),rmax_P,Q,sqrt(0.05),C,D)

## Now optimizing from somewhere random

theta_start = c(runif(1,0.5,5),runif(1,0.01,10),runif(1,0.05,1),runif(1,0.1,1),runif(1,1,15),runif(1,0.05,1),runif(1,0.01,5),runif(1,1,5),runif(1,0.05,1))

p_opt<-try(optim(theta_start, logLik, y=data1,method="BFGS",hessian=T), silent=TRUE)
if(is(p_opt, 'try-error'))
{  
  C_est[krep] = NA
  D_est[krep] = NA 
  break
} else {
  C_est[krep] = p_opt$par[7]
  D_est[krep] = p_opt$par[8]
}

## p_opt<-optim(theta_start, logLik, y=data,hessian=T) # Nelder-Mead

### --------------- Without the KR data ------------------------
theta_start = theta_start[-9]#one parameter less

p_opt<-try(optim(theta_start, logLik_FRwoutNoise, y=data1,method="BFGS",hessian=T), silent=TRUE)
if(is(p_opt, 'try-error'))
{  
  C_est_woutFR[krep] = NA
  D_est_woutFR[krep] = NA 
  break
} else {
  C_est_woutFR[krep] = p_opt$par[7]
  D_est_woutFR[krep] = p_opt$par[8]
}

## p_opt<-optim(theta_start, logLik_FRwoutNoise, y=data,hessian=T) # Nelder-Mead 

} # end of loop on repeats

pdf(file = "Estimators_BFGS_perturbedFP.pdf")
par(mfrow=c(2,2))
hist(C_est,breaks=30,xlim=c(0,50))
abline(v=C,col="red",lwd=2)
hist(D_est,breaks=30,xlim=c(0,60))
abline(v=D,col="red",lwd=2)

hist(C_est_woutFR,breaks=30)
abline(v=C,col="red",lwd=2)
hist(D_est_woutFR,breaks=30)
abline(v=D,col="red",lwd=2)
dev.off()

pdf(file = "Estimators_BFGS_perturbedFP_zoom.pdf")
par(mfrow=c(2,2))
hist(C_est,breaks=30,xlim=c(0,10))
abline(v=C,col="red",lwd=2)
hist(D_est,breaks=30,xlim=c(0,10))
abline(v=D,col="red",lwd=2)

hist(C_est_woutFR,breaks=30,xlim=c(0,10))
abline(v=C,col="red",lwd=2)
hist(D_est_woutFR,breaks=30,xlim=c(0,10))
abline(v=D,col="red",lwd=2)
dev.off()


#######################################################################################
#### ------------------- Now dataset 2 - noisy limit cycles ---------------------- ####
#######################################################################################

### Parameters of the model - fixed point parameter set
K<-1			# threshold dd 
rmax_V<-2			# Max AVERAGE growth rate (thus not a true max...)
rmax_P<-0.5
sigma2.proc<-0.05	
### FR and predator parameters
C<-15
D<-0.25
Q<-10

theta_true1  = c(rmax_V,1/K,sqrt(0.05),rmax_P,Q,sqrt(0.05),C,D,sqrt(0.05))
theta_true2  = c(rmax_V,1/K,sqrt(0.05),rmax_P,Q,sqrt(0.05),C,D)

C_est = D_est = C_est_woutFR = D_est_woutFR = rep(NA,100) #Initializing estimators of C and D

for (krep in 1:100){
  
  # switch to noisy LC dataset
  data2 = read.csv(paste(file_path2,"predatorPrey_withGaussianFR",krep,".csv",sep=""))
  names(data2)=c("Time","n","p","KR")
  data2[,1]<-NULL
  head(data2)
  
  ## Now optimizing from somewhere random
  
  theta_start = c(runif(1,0.5,5),runif(1,0.01,10),runif(1,0.05,1),runif(1,0.1,1),runif(1,1,15),runif(1,0.05,1),runif(1,0.01,25),runif(1,0.01,5),runif(1,0.05,1))
  
  p_opt<-try(optim(theta_start, logLik, y=data2,method="BFGS",hessian=T), silent=TRUE)
  if(is(p_opt, 'try-error'))
  {  
    C_est[krep] = NA
    D_est[krep] = NA 
    break
  } else {
    C_est[krep] = p_opt$par[7]
    D_est[krep] = p_opt$par[8]
  }
  
  ## p_opt<-optim(theta_start, logLik, y=data,hessian=T) # Nelder-Mead
  
  ### --------------- Without the KR data ------------------------
  theta_start = theta_start[-9]#one parameter less
  
  p_opt<-try(optim(theta_start, logLik_FRwoutNoise, y=data2,method="BFGS",hessian=T), silent=TRUE)
  if(is(p_opt, 'try-error'))
  {  
    C_est_woutFR[krep] = NA
    D_est_woutFR[krep] = NA 
    break
  } else {
    C_est_woutFR[krep] = p_opt$par[7]
    D_est_woutFR[krep] = p_opt$par[8]
  }
  
  ## p_opt<-optim(theta_start, logLik_FRwoutNoise, y=data,hessian=T) # Nelder-Mead 
  
} # end of loop on repeats

par(mfrow=c(2,2))

hist(C_est)
hist(D_est)

hist(C_est_woutFR)
hist(D_est_woutFR)



pdf(file = "Estimators_BFGS_noisyLC.pdf")
par(mfrow=c(2,2))
hist(C_est,breaks=30,xlim=c(0,50))
abline(v=C,col="red",lwd=2)
hist(D_est,breaks=30,xlim=c(0,60))
abline(v=D,col="red",lwd=2)

hist(C_est_woutFR,breaks=30,xlim=c(0,50))
abline(v=C,col="red",lwd=2)
hist(D_est_woutFR,breaks=30)
abline(v=D,col="red",lwd=2,xlim=c(0,60))
dev.off()

pdf(file = "Estimators_BFGS_noisyLC_zoom.pdf")
par(mfrow=c(2,2))
hist(C_est,breaks=30,xlim=c(0,10))
abline(v=C,col="red",lwd=2)
hist(D_est,breaks=30,xlim=c(0,10))
abline(v=D,col="red",lwd=2)

hist(C_est_woutFR,breaks=30,xlim=c(0,10))
abline(v=C,col="red",lwd=2)
hist(D_est_woutFR,breaks=30,xlim=c(0,10))
abline(v=D,col="red",lwd=2)
dev.off()
