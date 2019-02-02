#### 15/03/2018 -- FBarraquand and OGimenez
#### Computation of the Hessian for the discrete-time Leslie-May model
### FB 16/01/2019 -- edited so that functional response timing match math model (no delays)

# FB 02/02/2019 -- computing Fisher Information Matrix

rm(list=ls())
graphics.off()

# taking simuls from 
file_path = "../simulations/small_noise_on_FR/parameter_sets/perturbed_fixed_point/T=100/"
source('LikelihoodFunctions.R')

### Dataset to fit of the form
### data = cbind(log(N),log(P),FR) 

library("numDeriv")

### Parameters
K<-1			# threshold dd 
rmax_V<-2			# Max AVERAGE growth rate (thus not a true max...)
rmax_P<-0.5
sigma2.proc<-0.05	
### FR and predator parameters
C<-2.5
D<-1
Q<-10

theta_true1  = c(rmax_V,1/K,sqrt(0.05),rmax_P,Q,sqrt(0.05),C,D,sqrt(0.05))
theta_true2  = c(rmax_V,1/K,sqrt(0.05),rmax_P,Q,sqrt(0.05),C,D)   ### Removing last parameter for test with LL_FRwoutNoise

n.repeats <- 100
FIM1 = matrix(0,9,9) #initialize FIM
FIM2 = matrix(0,8,8)

for (krep in 1:100){

  data = read.csv(paste(file_path,"predatorPrey_withGaussianFR",krep,".csv",sep=""))
  data[1]<-NULL
  print(krep) 
  
  if (sum(sapply(data,is.finite))<300){message(paste0("Pb with rep:",krep))}
  
  LL = function(theta){
  return(logLik(theta,data))
  }
  LL1 = LL(theta_true1)
  print(LL1)
  
  if (is.finite(LL1)){ # removes automatically problematic time series and bad estimates
  hessian_thetaTrue=hessian(LL,theta_true1)
  FIM1 = FIM1 + hessian_thetaTrue
  }
  
  LL_FRwoutNoise = function(theta){
  return(logLik_FRwoutNoise(theta,data))
  }
  LL2 = LL_FRwoutNoise(theta_true2)
  print(LL2)
  
  if (is.finite(LL2)){
  hessian_thetaTrue_FRwoutNoise=hessian(LL_FRwoutNoise,theta_true2)
  FIM2 = FIM2 + hessian_thetaTrue_FRwoutNoise
  }  
}

FIM1 = FIM1/n.repeats
FIM2 = FIM2/n.repeats

eigen(FIM1)
eigen(FIM2)

### NB we average the matrices not the eigenvalues 
### the eigenvalues themselves seem noisy for T=100

### Using Viallefond et al.'s criteria
lambda1 = max(eigen(FIM1)$values)
q = 9
eps = 10^(-9)
eigen(FIM1)$values>q*lambda1*eps
# [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
eigen(FIM1)$values
# [1] 3935.3588217 3897.2837594 3881.6778990 3367.7063613 1982.7696500 1618.3102018    7.7358066
# [8]    4.2169881    0.1838534

lambda1 = max(eigen(FIM2)$values)
q = 8
eps = 10^(-9)
eigen(FIM2)$values>q*lambda1*eps
#[1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
eigen(FIM2)$values
# [1] 3.958224e+03 3.897284e+03 3.375020e+03 1.982770e+03 4.355443e+00 9.727962e-01 1.838534e-01
# [8] 4.096381e-02

# plot both sets of eigenvalues
eigentable = cbind(eigen(FIM1)$values,c(eigen(FIM2)$values,NA))
dotchart(eigentable,gcolor=c("blue","red"),pch=19)

pdf("FIMeigenvalues_perturbedfixedpoint_T100.pdf",width=10,height=5)
eigentable = cbind(log10(eigen(FIM1)$values),c(log10(eigen(FIM2)$values),NA))
dotchart(eigentable,gcolor=c("blue","red"),pch=19,xlab="log10(eigenvalue)",ylab="Rank eigenvalue",xlim=c(-2,4))
dev.off()
# Not very different -- I should do the same for the T=1000 computation