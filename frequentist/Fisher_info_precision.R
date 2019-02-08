#### 15/03/2018 -- FBarraquand and OGimenez
#### Computation of the Hessian for the discrete-time Leslie-May model
### FB 16/01/2019 -- edited so that functional response timing match math model (no delays)

# FB 02/02/2019 -- computing Fisher Information Matrix
# FB 08/02/2019 -- computed the FIM for a time series perturbed by a (very small) Obs Error
# Objective was to see if this can be used to compute a precision for the FIM, but no (FIM changes too much)

rm(list=ls())
graphics.off()

# taking simuls from 
file_path = "../simulations/small_noise_on_FR/parameter_sets/perturbed_fixed_point/T=100/"
file_path2 = "../simulations/small_noise_on_FR/parameter_sets/perturbed_fixed_point/T=100_wObsError/"

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

FIM1_ObsError = matrix(0,9,9) #initialize FIM
FIM2_ObsError = matrix(0,8,8)

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
  
  ### Second dataset
  data = read.csv(paste(file_path2,"predatorPrey_withGaussianFR",krep,".csv",sep=""))
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
    FIM1_ObsError = FIM1_ObsError + hessian_thetaTrue
  }
  
  LL_FRwoutNoise = function(theta){
    return(logLik_FRwoutNoise(theta,data))
  }
  LL2 = LL_FRwoutNoise(theta_true2)
  print(LL2)
  
  if (is.finite(LL2)){
    hessian_thetaTrue_FRwoutNoise=hessian(LL_FRwoutNoise,theta_true2)
    FIM2_ObsError = FIM2_ObsError + hessian_thetaTrue_FRwoutNoise
  }  
}

FIM1 = FIM1/n.repeats
FIM2 = FIM2/n.repeats

FIM1_ObsError = FIM1_ObsError/n.repeats
FIM2_ObsError = FIM2_ObsError/n.repeats

Deviation_FIM1 = FIM1_ObsError - FIM1
Deviation_FIM2 = FIM2_ObsError - FIM2

rFIM1 = round(FIM1_ObsError,digits=6)
rFIM2 = round(FIM1_ObsError,digits=6)
write.csv(rFIM1,file="FIM1_ObsError.csv")
write.csv(rFIM2,file="FIM2_ObsError.csv")

# Condition number
1 / (norm(FIM1) * norm(solve(FIM1)))
# [1] 4.440851e-05
1 / (norm(FIM2) * norm(solve(FIM2)))
# [1] 7.847715e-06

# Condition number
1 / (norm(FIM1_ObsError) * norm(solve(FIM1_ObsError)))
# [1] 4.34616e-05
1 / (norm(FIM2_ObsError) * norm(solve(FIM2_ObsError)))
# [1] 4.280249e-06

# Compare to the deviation which is the acceptable error on the FIM matrix
round(Deviation_FIM1,digits=3)
round(FIM1,digits=3)
round(FIM1_ObsError,digits=3)
max(abs(Deviation_FIM1))
#[1] 59.55259
mean(abs(Deviation_FIM1))
#[1] 4.210814

# Using the deviation to the original FIM does not provide a meaningful bound to compare to the relative condition number

############### Alternative lead to understand the FIM eigenvalues -- 08/02/201 ##################
############### Tentative null model for eigenvalues of similarly distributed matrices ###########

### Let's think about a null model for these matrices where everything goes smoothly
Sigma = 10000*diag(9) 
X=rWishart(1, 18, Sigma)
X=X[,,1]
eigen(X)$values
#429022.84 296022.00 242272.82 164730.08 104347.50  85050.62  68206.50  56516.72  38392.17

#Let's introduce some correlation
sc=100
block1 = matrix(rnorm(9,sc,15),3,3)
block2 = matrix(rnorm(9,sc,15),3,3)
block3 = matrix(rnorm(9,sc,15),3,3)
null = matrix(0,3,3)
a=cbind(block1,null,null)
b=cbind(null,block2,null)
c=cbind(null,null,block3)
A = rbind(a,b,c)
eigen(A)$values
Sigma = t(A)%*%A #gets it positive definite
eigen(Sigma)$values

X=rWishart(1, 18, Sigma) 
X=X[,,1]
eigen(X)$values
### Looks like this modular variance-covariance matrix generates skew in the eigenvalues

### Here A would be equivalent to values of the partial LL / partial theta_i? 
### Hum, not really. Some argument below
library(mvtnorm)
sc=1
block1 = matrix(rnorm(9,sc,1),3,3)
block2 = matrix(rnorm(9,sc,1),3,3)
block3 = matrix(rnorm(9,sc,1),3,3)
null = matrix(0,3,3)
a=cbind(block1,null,null)
b=cbind(null,block2,null)
c=cbind(null,null,block3)
A = rbind(a,b,c)
A[lower.tri(A)]=t(A)[lower.tri(A)]
A
B=0.1*t(A)%*%A # to render it positive definite
X=rmvnorm(1,rep(10,9),B)
False_FIM = t(X) %*% X ## does not work, not modular like the real FIM
## Perhaps I should get back to the differentiation of the log-likelihood