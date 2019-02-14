#### 15/03/2018 -- FBarraquand and OGimenez
#### Computation of the Hessian for the discrete-time Leslie-May model
### FB 16/01/2019 -- edited so that functional response timing match math model (no delays)

# FB 02/02/2019 -- computing Fisher Information Matrix

rm(list=ls())
graphics.off()

# taking simuls from 
file_path = "../simulations/small_noise_on_FR/parameter_sets/noisy_cycles/T=100/"
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
C<-15
D<-0.25
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
# [1] 97319.954453  3906.347199  3846.735087  3841.058487  2786.549700  1962.214666   642.783411
# [8]    80.434582     1.094164

lambda1 = max(eigen(FIM2)$values)
q = 8
eps = 10^(-9)
eigen(FIM2)$values>q*lambda1*eps
#[1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
eigen(FIM2)$values
# [1] 46327.223462 10685.933317  3841.058487  2467.297875  1962.214666    55.254972     2.523326
# [8]     1.094164

# plot both sets of eigenvalues
eigentable = cbind(eigen(FIM1)$values,c(eigen(FIM2)$values,NA))
dotchart(eigentable,gcolor=c("blue","red"),pch=19)

pdf("FIMeigenvalues_noisyLC_T100.pdf",width=10,height=5)
eigentable = cbind(log10(eigen(FIM1)$values),c(log10(eigen(FIM2)$values),NA))
dotchart(eigentable,gcolor=c("blue","red"),pch=19,xlab="log10(eigenvalue)",ylab="Rank eigenvalue",xlim=c(-2,4))
dev.off()
# Not very different -- I should do the same for the T=1000 computation

det(FIM1)
#[1] 1.737478e+27
det(FIM2)
#[1] 1.404414e+21

rFIM1 = round(FIM1,digits=6)
rFIM2 = round(FIM2,digits=6)
write.csv(rFIM1,file="FIM1_noisyLC.csv")
write.csv(rFIM2,file="FIM2_noisyLC.csv")

# Condition number
1 / (norm(FIM1) * norm(solve(FIM1)))
#[1] 1.007347e-05
1 / (norm(FIM2) * norm(solve(FIM2)))
# [1] 2.183349e-05

1 / (norm(rFIM1) * norm(solve(rFIM1)))
# [1] 1.007348e-05
1 / (norm(rFIM2) * norm(solve(rFIM2)))
# [1] 2.18335e-05
### Good, robust to rounding

### One says that a matrix is ill-conditioned 
### when the condition number is larger than the precision of the matrix entries
# Is this the same thing? Should be. Check again SVD btw. 
max(eigen(FIM1)$values)/min(eigen(FIM1)$values)
min(eigen(FIM1)$values)/max(eigen(FIM1)$values) # yep - to the errors. 
min(eigen(FIM2)$values)/max(eigen(FIM2)$values) # yep - to the errors. 

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
rownames(Rho) = c("r_V","1/K_V","sigma1","r_P","Q","sigma2","C","D","sigma3")
colnames(Rho) = c("r_V","1/K_V","sigma1","r_P","Q","sigma2","C","D","sigma3")
corrplot(Rho, method="circle")
pdf("InverseFIM1_noisyLC_corrplot.pdf",width=6,height=6)
corrplot(Rho, method = "number")
dev.off()

Sigma = solve(FIM2)
Rho = cov2cor(Sigma) #correlation matrix
par(mfrow=c(1,1))
rownames(Rho) = c("r_V","1/K_V","sigma1","r_P","Q","sigma2","C","D")
colnames(Rho) = c("r_V","1/K_V","sigma1","r_P","Q","sigma2","C","D")
corrplot(Rho, method="circle")
pdf("InverseFIM2_noisyLC_corrplot.pdf",width=6,height=6)
corrplot(Rho, method = "number")
dev.off()
