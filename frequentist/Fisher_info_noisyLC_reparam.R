#### 15/03/2018 -- FBarraquand and OGimenez
#### Computation of the Hessian for the discrete-time Leslie-May model
### FB 16/01/2019 -- edited so that functional response timing match math model (no delays)

# FB 02/02/2019 -- computing Fisher Information Matrix
# FB 17/02/2019 -- with reparameterized model to avoid correlations

rm(list=ls())
graphics.off()

# taking simuls from 
file_path = "../simulations/small_noise_on_FR/parameter_sets/noisy_cycles/T=100/"
source('LikelihoodFunctions.R')

### Dataset to fit of the form
### data = cbind(log(N),log(P),FR) 

library("numDeriv")

### Parameters 
# Have been modified to take into account reparam
rmax_V<-2			# Max AVERAGE growth rate (thus not a true max...)
K<-1*(exp(rmax_V)-1)	# true carrying capacity - so that threshold for dd = 1 as in previous code
rmax_P<-0.5
sigma2.proc<-0.05	
### FR and predator parameters
C<-15
D<-0.25
a = C/D #reparam to classical Holling type II
h = 1/C #reparam, handling time
Q = 10/(exp(rmax_P)-1)

theta_true1  = c(rmax_V,K,sqrt(0.05),rmax_P,Q,sqrt(0.05),a,h,sqrt(0.05))
theta_true2  = c(rmax_V,K,sqrt(0.05),rmax_P,Q,sqrt(0.05),a,h)   ### Removing last parameter for test with LL_FRwoutNoise

n.repeats <- 100
FIM1 = matrix(0,9,9) #initialize FIM
FIM2 = matrix(0,8,8)

for (krep in 1:100){

  data = read.csv(paste(file_path,"predatorPrey_withGaussianFR",krep,".csv",sep=""))
  data[1]<-NULL
  print(krep) 
  
  if (sum(sapply(data,is.finite))<300){message(paste0("Pb with rep:",krep))}
  
  LL = function(theta){
  return(logLik_reparam(theta,data))
  }
  LL1 = LL(theta_true1)
  print(LL1)
  
  if (is.finite(LL1)){ # removes automatically problematic time series and bad estimates
  hessian_thetaTrue=hessian(LL,theta_true1)
  FIM1 = FIM1 + hessian_thetaTrue
  }
  
  LL_FRwoutNoise = function(theta){
  return(logLik_FRwoutNoise_reparam(theta,data))
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
# [1] 5.894246e+07 3.904568e+03 3.846748e+03 3.841592e+03 7.602579e+02 2.894537e+02 1.840307e+01 1.188629e+00
# [9] 9.319847e-01

lambda1 = max(eigen(FIM2)$values)
q = 8
eps = 10^(-9)
eigen(FIM2)$values>q*lambda1*eps
#[1]  TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
eigen(FIM2)$values
# [1] 2.966232e+05 4.620387e+04 3.841592e+03 7.602579e+02 1.763334e+02 1.360253e+01 1.188629e+00
# [8] 6.845582e-02

# plot both sets of eigenvalues
eigentable = cbind(eigen(FIM1)$values,c(eigen(FIM2)$values,NA))
dotchart(eigentable,gcolor=c("blue","red"),pch=19)

pdf("FIMeigenvalues_noisyLC_T100_reparam.pdf",width=10,height=5)
eigentable = cbind(log10(eigen(FIM1)$values),c(log10(eigen(FIM2)$values),NA))
dotchart(eigentable,gcolor=c("blue","red"),pch=19,xlab="log10(eigenvalue)",ylab="Rank eigenvalue",xlim=c(-2,4))
dev.off()
# Not very different -- I should do the same for the T=1000 computation

det(FIM1)
#[1]  1.525779e+25
det(FIM2)
#[1] 7.812083e+18

rFIM1 = round(FIM1,digits=6)
rFIM2 = round(FIM2,digits=6)
write.csv(rFIM1,file="FIM1_noisyLC_reparam.csv")
write.csv(rFIM2,file="FIM2_noisy_LC_reparam.csv")

# Condition number
1 / (norm(FIM1) * norm(solve(FIM1)))
# [1] 1.580722e-08
1 / (norm(FIM2) * norm(solve(FIM2)))
# [1]  2.173759e-07

### One says that a matrix is ill-conditioned 
### when the condition number is larger than the precision of the matrix entries
# Is this the same thing? Should be. Check again SVD btw. 
max(eigen(FIM1)$values)/min(eigen(FIM1)$values)
min(eigen(FIM1)$values)/max(eigen(FIM1)$values) # yep - to the errors. 
min(abs(eigen(FIM2)$values))/max(abs(eigen(FIM2)$values)) # yep - to the errors. 

### Could I use a Hilbert matrix as a gold standard of ill-conditioning? 
hilbert = function(k){
  B = matrix(0,k,k)
  for (i in 1:k){
    for (j in 1:k){
      B[i,j] = 1 / (i+j-1)
    }
  }
  return(B)
}

A = hilbert(9)
lambda1 = max(eigen(A)$values)
q = 9
eps = 10^(-9)
eigen(A)$values>q*lambda1*eps
#  Criteria works for these.
det(A)
#[1] 9.720265e-43
#but beware
# det(1000*A)
# [1] 9.720274e-16
1 / (norm(A * norm(solve(A))))
# [1] 9.093794e-13
1 / (norm(1000* A * norm(solve(1000*A))))
# [1] 9.093802e-13

### Construct correlation matrices
library(corrplot)

Sigma = solve(FIM1)
Rho = cov2cor(Sigma) #correlation matrix
par(mfrow=c(1,1))
rownames(Rho) = c("r_V","K","sigma1","r_P","Q","sigma2","a","h","sigma3")
colnames(Rho) = c("r_V","K","sigma1","r_P","Q","sigma2","a","h","sigma3")
corrplot(Rho, method="circle")
pdf("InverseFIM1_noisyLC_corrplot_reparam.pdf",width=6,height=6)
corrplot(Rho, method = "number")
dev.off()

Sigma = solve(FIM2) #plenty of negative values there
Rho = cov2cor(Sigma) #correlation matrix
par(mfrow=c(1,1))
rownames(Rho) = c("r_V","K","sigma1","r_P","Q","sigma2","a","h")
colnames(Rho) = c("r_V","K","sigma1","r_P","Q","sigma2","a","h")
corrplot(Rho, method="circle")
pdf("InverseFIM2_noisyLC_corrplot_reparam.pdf",width=6,height=6)
corrplot(Rho, method = "number")
dev.off()
