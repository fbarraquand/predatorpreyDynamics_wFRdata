#### 15/03/2018 -- FBarraquand and OGimenez
#### Computation of the Hessian for the discrete-time Leslie-May model


### --- Questions ------------------------------------------------------------------#### 
### Better to do this numerically or symbolically? (e.g. later in MAPLE?) 
### If numerically, how far into parameter space, i.e. what precision for the grid? 
### --------------------------------------------------------------------------------####

rm(list=ls())
graphics.off()

##################### Helper function #####################################################
### By OG 
logprot <- function(v){
eps <- 2.2204e-016
u <- log(eps) * (1+vector(length=length(v)))
index <- (v>eps)
u[index] <- log(v[index])
u
}

##################### Simulation used #####################################################

### Model used (use stored simulated data later on, just tryouts for now)
n.years<-1000  	# Long time series to look at structural identifiability
N1<-1			# Initial pop size
P1<-0.1
K<-1			# threshold dd 
beta<-1			# density-dependence exponent
rmax_V<-2			# Max AVERAGE growth rate (thus not a true max...)
rmax_P<-0.5
sigma2.proc<-0.05		# worked well with 0.005
# Process sigma on the log-scale, use the Peretti et al. value. 0.005


### FR and predator parameters
C<-2.5
D<-1
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
 
  FR[t+1]<-(C*N[t]/(D+N[t])) + FRnoise[t+1]
  N[t+1]<-N[t]*(exp(rV[t])/(1+(N[t]/K)^beta))*exp(-FR[t]*P[t]/N[t])
  P[t+1]<-P[t]*exp(rP[t])/(1+P[t]*Q/N[t])
}
## Plotting time series of abundances and FR
par(mfrow=c(2,2))
plot(1:n.years,N,type="b")
plot(1:n.years,P,type="b")
#curve(dbeta(x,a,b),from=0, to=1)
plot(N,FR)

### Produce dataset to fit
data = cbind(log(N),log(P),FR) 

##################### Define likelihood and RSS for both models with and without noisy FRs ######################################

################# Define the likelihood ########################################
logLik=function(theta,y){
#### y is the data with n rows and 3 columns // log-abundance data + FR data

#### Parameters
# theta1 = r 
# theta2 = gamma
# theta3 = sigma1
# theta4 = s
# theta5 = q 
# theta6 = sigma2
# theta7 = C
# theta8 = D
# theta9 = sigma3

n=nrow(y)
ll = 0.0 ### Or p_1(a_1|x_1) p(x_1)
for (t in 2:n){
N_t = exp(y[t-1,1])
P_t = exp(y[t-1,2]) 

######## Error that was previously there!! N_t/P_t instead P_t/N_t ###
### mu1 = y[t-1,1] + theta[1] - log(1+theta[2]*N_t) - y[t-1,3]*N_t/P_t
######################################################################
mu1 = y[t-1,1] + theta[1] - logprot(1+theta[2]*N_t) - y[t-1,3]*P_t/N_t
mu2 = y[t-1,2] + theta[4] - logprot((1+theta[5]*P_t/N_t))
mu3 = (theta[7]*N_t)/(theta[8] + N_t)
#ll= ll + log(dnorm(y[t,1], mu1, theta[3])) +  log(dnorm(y[t,2], mu2, theta[6])) + log(dnorm(y[t,3], mean = mu3, sd = theta[9]))
# we have log(0) problem
d1=dnorm(y[t,1], mu1, theta[3],log=T) ## directly asking for the log avoids problems
d2=dnorm(y[t,2], mu2, theta[6],log=T)
d3=dnorm(y[t,3], mu3, theta[9],log=T)
ll=ll+d1+d2+d3
}
return(-ll)
}


############### Working directly with the sum of squares ######################
RSS=function(theta,y){
  #### y is the data with n rows and 3 columns // log-abundance data + FR data
  
  #### Parameters
  # theta1 = r 
  # theta2 = gamma
  # theta3 = s
  # theta4 = q 
  # theta5 = C
  # theta6 = D
  # not the same theta
  
  n=nrow(y)
  rss = 0.0 ### Or p_1(a_1|x_1) p(x_1)
  for (t in 2:n){
    N_t = exp(y[t-1,1])
    P_t = exp(y[t-1,2]) 
    ############## Correction of error ###########################################
    ## mu1 = y[t-1,1] + theta[1] - log(1+theta[2]*N_t) - y[t-1,3]*N_t/P_t # error
    mu1 = y[t-1,1] + theta[1] - logprot(1+theta[2]*N_t) - y[t-1,3]*P_t/N_t # corrected
    mu2 = y[t-1,2] + theta[3] - logprot((1+theta[4]*P_t/N_t))
    mu3 = (theta[5]*N_t)/(theta[6] + N_t)
    rss=rss+(y[t,1] - mu1)^2+(y[t,2]-mu2)^2+(y[t,3]-mu3)^2
  }
  return(rss)
}

############################### LL ##################################################

################# Define the likelihood ########################################

logLik_FRwoutNoise=function(theta,y){
#### y is the data with n rows and 3 columns // log-abundance data + FR data

#### Parameters
# theta1 = r 
# theta2 = gamma
# theta3 = sigma1
# theta4 = s
# theta5 = q 
# theta6 = sigma2
# theta7 = C
# theta8 = D
# theta9 = sigma3

n=nrow(y)
ll = 0.0 ### Or p_1(a_1|x_1) p(x_1)
for (t in 2:n){
N_t = exp(y[t-1,1])
P_t = exp(y[t-1,2]) 

######## Error that was previously there!! N_t/P_t instead P_t/N_t ###
### mu1 = y[t-1,1] + theta[1] - log(1+theta[2]*N_t) - y[t-1,3]*N_t/P_t
######################################################################
mu1 = y[t-1,1] + theta[1] - logprot(1+theta[2]*N_t) - (theta[7]*P_t)/(theta[8] + N_t)
mu2 = y[t-1,2] + theta[4] - logprot((1+theta[5]*P_t/N_t))
#ll= ll + log(dnorm(y[t,1], mu1, theta[3])) +  log(dnorm(y[t,2], mu2, theta[6])) + log(dnorm(y[t,3], mean = mu3, sd = theta[9]))
# we have log(0) problem
d1=dnorm(y[t,1], mu1, theta[3],log=T) ## directly asking for the log avoids problems
d2=dnorm(y[t,2], mu2, theta[6],log=T)
#d3=dnorm(y[t,3], mu3, theta[9],log=T)
ll=ll+d1+d2#+d3
}
return(-ll)
}


############### Working directly with the sum of squares ######################
RSS_FRwoutNoise=function(theta,y){
  #### y is the data with n rows and 3 columns // log-abundance data + FR data
  
  #### Parameters
  # theta1 = r 
  # theta2 = gamma
  # theta3 = s
  # theta4 = q 
  # theta5 = C
  # theta6 = D
  # not the same theta
  
  n=nrow(y)
  rss = 0.0 ### Or p_1(a_1|x_1) p(x_1)
  for (t in 2:n){
    N_t = exp(y[t-1,1])
    P_t = exp(y[t-1,2]) 
    ############## Correction of error ###########################################
    ## mu1 = y[t-1,1] + theta[1] - log(1+theta[2]*N_t) - y[t-1,3]*N_t/P_t # error
    mu1 = y[t-1,1] + theta[1] - logprot(1+theta[2]*N_t) -  (theta[5]*P_t)/(theta[6] + N_t) # corrected
    mu2 = y[t-1,2] + theta[3] - logprot((1+theta[4]*P_t/N_t))
    rss=rss+(y[t,1] - mu1)^2+(y[t,2]-mu2)^2#+(y[t,3]-mu3)^2
  }
  return(rss)
}

################################ Computation of MLE ##########################################################
################################ + computation of the Hessian at the MLE #####################################

######### Estimation starting at the MLE ###########
theta_true  = c(rmax_V,1/K,sqrt(0.05),rmax_P,Q,sqrt(0.05),C,D,sqrt(0.05))
p_opt<-optim(theta_true, logLik, y=data,method="BFGS",hessian=T)
p_opt$par
theta_true

######### Estimation starting away from MLE ########
theta_init = theta_true + rnorm(9,0,sd=0.05)
p_opt<-optim(theta_init, logLik, y=data,method="BFGS",hessian=T,control=list(fnscale=-1))
p_opt$par
theta_true
### Ask optim for the Hessian
hessian=p_opt$hessian
print(hessian, dig = 9) ## or use scipen
options(scipen = 999)
print(hessian, dig = 9)
#revert back
options(scipen = 0)

eigen(hessian) ## clearly no zeroes there? 
# $values
# [1] 41617.593143 40362.310251 37309.863813 33668.343616 18679.557058 18079.111183    93.263279    12.579849     1.507702

lambda1 = max(eigen(hessian)$values)
9*lambda1*10^(-9) #way below lowest eigenvalue

# should be reproduced for several simulations later on. 

## --- old comments on previous simus ---- #
# There is however the question of whether this matrix is positive definite. 
# (which happens for the Hessian when all eigenvalues are positive) -> I'm not sure in fact. 
# The last eigenvalue is negative. On the other hand it is very small so we can consider it zero. 
## ---------------------------------------- #

# Looking at the correlations in the matrix 

Sigma =solve(hessian) #variance-covariance matrix
Rho = cov2cor(Sigma) #correlation matrix
library(corrplot)
par(mfrow=c(1,1))
rownames(Rho) = c("r_V","1/K_V","sigma1","r_P","Q","sigma2","C","D","sigma3")
colnames(Rho) = c("r_V","1/K_V","sigma1","r_P","Q","sigma2","C","D","sigma3")
corrplot(Rho, method="circle")
corrplot(Rho, method = "number")
# There are significant correlations between pairs of parameters (same as in the Bayesian analysis)

######## Estimation based on RSS ##################

### New restricted theta
theta_true  = c(rmax_V,1/K,rmax_P,Q,C,D)
p_opt<-optim(theta_true, RSS,y=data,hessian=T)
p_opt$par
theta_true
p_opt$hessian
eigen(p_opt$hessian) ### clearly interesting and non zero
# $values
# [1] 3334.2367985 2000.8569587 1736.6593414    8.9622937    1.4761840    0.1603394
# 
# $vectors
# [,1]          [,2]          [,3]          [,4]          [,5]          [,6]
# [1,]  7.739897e-01 -1.080998e-12  4.475014e-12  2.216871e-10 -6.331982e-01  1.602045e-09
# [2,] -6.331982e-01  2.257642e-12  1.515849e-12  2.719858e-10 -7.739897e-01  1.961642e-09
# [3,]  2.183674e-12  9.992858e-01  2.420841e-11 -9.492240e-12 -9.464525e-11 -3.778865e-02
# [4,] -2.225576e-12 -3.778865e-02 -3.714140e-12 -2.153198e-10 -2.530934e-09 -9.992858e-01
# [5,]  2.562748e-12  2.312261e-11 -9.299362e-01  3.677209e-01  1.253021e-10 -7.665235e-11
# [6,] -3.279388e-13 -7.693068e-12  3.677209e-01  9.299362e-01  3.277751e-10 -2.014524e-10

# Let's try from somewhere slightly different
theta_init = theta_true + rnorm(6,0,sd=0.05)
p_opt<-optim(theta_init, RSS, y=data,hessian=T)
p_opt$par
theta_true
p_opt$hessian
eigen(p_opt$hessian) ### clearly interesting and non zero
#$values
# [1] 3352.8312867 2000.8570244 1735.5738470    8.9450653    1.2821586    0.1615576
# 
# $vectors
# [,1]          [,2]          [,3]         [,4]          [,5]          [,6]
# [1,]  7.718547e-01  2.191974e-12 -1.043386e-12 2.212856e-10 -6.357990e-01 -1.340079e-09
# [2,] -6.357990e-01 -8.543305e-14 -4.183676e-12 2.691969e-10 -7.718547e-01 -1.625599e-09
# [3,]  1.715033e-12 -9.992857e-01 -8.596346e-12 1.284783e-11  7.828529e-11 -3.778910e-02
# [4,] -8.572868e-13  3.778910e-02 -3.438708e-12 2.772205e-10  2.105292e-09 -9.992857e-01
# [5,]  1.855122e-12  8.737344e-12 -9.298957e-01 3.678232e-01  1.317965e-10  1.055713e-10
# [6,] -3.522451e-13 -9.151013e-13  3.678232e-01 9.298957e-01  3.226127e-10  2.566702e-10


###################################################################################################
################### Now for the model where the FR is specified without the noise #################
###################################################################################################
# Revert back to full theta_true with 9 elements
# Nope!! There are only 8 elements here, no sigma3!!

theta_true  = c(rmax_V,1/K,sqrt(0.05),rmax_P,Q,sqrt(0.05),C,D)
p_opt<-optim(theta_true, logLik_FRwoutNoise, method="BFGS",y=data,hessian=T)
p_opt$par
### BFGS
# [1]  1.5041455  0.5197072  0.2226370  0.5065558 10.3046050  0.2314249 33.2937812 66.2760413
### Nelder-Mead
# [1]  1.9914255  0.9680671  0.2224996  0.5107744 10.4210424  0.2310457  3.2480513  1.8051012
theta_true
# [1]  2.0000000  1.0000000  0.2236068  0.5000000 10.0000000  0.2236068  2.5000000  1.0000000
p_opt$hessian
# [,1]          [,2]          [,3]          [,4]          [,5]          [,6]          [,7]          [,8]
# [1,]  2.017933e+04 -1.744166e+04  1.334904e+02 -1.776357e-09  0.000000e+00  3.552714e-09 -9.778635e+02  4.602678e+02
# [2,] -1.744166e+04  1.508817e+04 -1.271656e+02 -3.552714e-09 -2.309264e-08 -3.552714e-09  8.390894e+02 -3.920669e+02
# [3,]  1.334904e+02 -1.271656e+02  4.062105e+04 -1.065814e-08  1.421085e-08  2.131628e-08 -9.458752e+00  5.322839e+00
# [4,] -1.776357e-09 -3.552714e-09 -1.065814e-08  1.871413e+04 -7.040348e+02  2.880262e+01  0.000000e+00 -1.598721e-08
# [5,]  0.000000e+00 -2.309264e-08  1.421085e-08 -7.040348e+02  2.795303e+01 -2.714926e+00  0.000000e+00 -1.776357e-08
# [6,]  3.552714e-09 -3.552714e-09  2.131628e-08  2.880262e+01 -2.714926e+00  3.761812e+04  3.552714e-09  3.552714e-09
# [7,] -9.778635e+02  8.390894e+02 -9.458752e+00  0.000000e+00  0.000000e+00  3.552714e-09  5.381699e+01 -2.630698e+01
# [8,]  4.602678e+02 -3.920669e+02  5.322839e+00 -1.598721e-08 -1.776357e-08  3.552714e-09 -2.630698e+01  1.318995e+01
eigen(p_opt$hessian)
# $values
# [1] 6.054577e+04 4.031416e+04 3.730971e+04 1.867950e+04 5.779300e+01 1.504246e+00 5.509579e-02
# [8] 5.308338e-05
# 
# $vectors
# [,1]          [,2]          [,3]          [,4]          [,5]          [,6]
# [1,]  5.764044e-01 -4.308947e-05  1.260089e-12 -4.883187e-13  8.171292e-01 -2.827274e-10
# [2,] -8.171589e-01  4.152968e-05  1.139785e-13 -3.081142e-13  5.764101e-01 -2.160426e-10
# [3,] -5.876738e-05 -1.000000e+00  1.064255e-11 -3.343097e-13 -1.125797e-05  7.715109e-14
# [4,]  3.740335e-14 -3.318773e-13 -5.056317e-05  9.992864e-01  1.529329e-11  3.777154e-02
# [5,]  1.881479e-13  9.401594e-14  5.322336e-06 -3.777154e-02  3.893140e-10  9.992864e-01
# [6,]  6.325077e-13 -1.064259e-11 -1.000000e+00 -5.072812e-05  1.096834e-12  3.408689e-06
# [7,] -2.766412e-03 -1.958303e-06  2.828092e-13  5.638389e-13 -6.198803e-03  2.891478e-09
# [8,]  1.282357e-03  3.971253e-07  5.726018e-13  7.743699e-13  3.588536e-03 -4.496595e-09
# [,7]          [,8]
# [1,] -7.586655e-03 -5.278368e-04
# [2,] -1.620123e-03 -3.734781e-04
# [3,]  2.203137e-06 -4.566506e-07
# [4,]  1.706311e-10  1.071605e-10
# [5,]  4.509435e-09  2.859934e-09
# [6,] -9.924139e-15  6.475621e-13
# [7,] -9.071098e-01  4.208392e-01
# [8,]  4.208225e-01  9.071350e-01

#### Nelder-Mead
# $values
# [1] 40627.4553092 37618.1668741 35311.2566499 18740.5788017    14.9681016     2.1552563     1.4647539    -0.2764503
# 
# $vectors
# [,1]          [,2]          [,3]          [,4]          [,5]          [,6]          [,7]          [,8]
# [1,] -2.609296e-02 -1.182946e-12  7.553842e-01 -1.117742e-13  4.023732e-01  4.688183e-01  9.838325e-09 -2.168389e-01
# [2,]  2.284980e-02  1.033504e-12 -6.531189e-01 -1.290197e-13  5.113473e-01  5.015502e-01  1.072476e-08 -2.447141e-01
# [3,] -9.993972e-01  7.189136e-12 -3.471459e-02  5.128058e-13  7.714872e-05  1.989578e-04  3.735806e-12 -9.823869e-05
# [4,]  4.875137e-13 -1.527377e-03  4.445520e-14  9.992918e-01 -5.565064e-11 -6.339617e-10  3.759676e-02  2.314851e-10
# [5,] -3.705045e-13  1.008307e-04  4.041092e-13 -3.759665e-02 -1.491911e-09 -1.685328e-08  9.992930e-01  6.134774e-09
# [6,] -7.131091e-12 -9.999988e-01 -1.818651e-12 -1.530088e-03 -3.962634e-14 -7.965266e-13  4.333509e-05  1.493318e-13
# [7,]  1.334819e-03 -4.242448e-14 -3.649754e-02  9.164045e-16 -6.113728e-01  7.270867e-01  9.445264e-09  3.102167e-01
# [8,] -6.481256e-04 -1.180791e-13  1.712617e-02 -8.180576e-13  4.503811e-01 -1.299384e-03 -4.829758e-09  8.926710e-01

######### Estimation starting away from MLE ########
theta_init = theta_true + rnorm(8,0,sd=0.05)
p_opt<-optim(theta_init, logLik_FRwoutNoise,method="BFGS", y=data,hessian=T)
p_opt$par
theta_true
theta_init
p_opt$hessian
eigen(p_opt$hessian)
### BFGS
#$values
#[1] 6.193873e+04 4.032005e+04 3.731044e+04 1.867967e+04 6.016344e+01 1.507778e+00 1.812701e-02
#[8] 9.703971e-06
### Looks like BFGS performs worse than Nelder-Mead without the FR data (see below as well)

### Nelder-Mead
#$values
#[1]  4.013941e+04  3.727161e+04  3.536148e+04  1.867194e+04  1.933835e+01  4.327331e+00  1.533097e+00 -4.063515e-02
### Seems better starting a little further away. May be need to be repeated

theta_init = theta_true + rnorm(8,0,sd=0.05)
p_opt<-optim(theta_init, logLik_FRwoutNoise,method="BFGS",y=data,hessian=T)
p_opt$par
theta_true
theta_init
p_opt$hessian
eigen(p_opt$hessian)
### BFGS
#$values
#[1] 6.037045e+04 4.031188e+04 3.731284e+04 1.868018e+04 5.741190e+01 1.513332e+00 5.721312e-02
#[8] 5.451893e-05
# $values
# [1] 40942.1328803 36579.3709513 33872.4695387 18511.1139219    19.6803090     2.9133301     1.2436616    -0.2908425
### Seems like quite of bit of these estimates of last eigenvalue are negative. 
diag(solve(p_opt$hessian)) ### produces a number of negative variance estimates. 
## Beware that we use -LL, so Hessian = observed FIM (Fisher Information Matrix)
## Loop over datasets to obtain the true, expected FIM?

######## Estimation from RSS ################
theta_true  = c(rmax_V,1/K,rmax_P,Q,C,D)
p_opt<-optim(theta_true,RSS_FRwoutNoise, y=data,hessian=T)
p_opt$par
theta_true
p_opt$hessian
eigen(p_opt$hessian) ### clearly interesting and non zero? (inluding the last one?)
# $values
# [1] 3344.31352989 2000.86370062    1.73515380    0.38035271    0.16121792   -0.00985483
# $vectors
# [,1]          [,2]          [,3]          [,4]          [,5]          [,6]
# [1,]  7.728650e-01 -6.748275e-12 -4.098421e-01  4.651301e-01  2.680788e-09  1.355107e-01
# [2,] -6.330323e-01  1.050155e-12 -5.502109e-01  5.191807e-01  3.818701e-09  1.642847e-01
# [3,] -5.951248e-12 -9.992841e-01 -2.332476e-10 -2.002480e-10 -3.783315e-02  7.145223e-10
# [4,]  3.866315e-13  3.783315e-02 -6.257578e-09 -5.249541e-09 -9.992841e-01  1.886142e-08
# [5,] -4.011429e-02  1.187677e-12  5.917052e-01  7.025929e-01 -1.481861e-08 -3.932418e-01
# [6,]  1.845706e-02 -1.824998e-12 -4.232994e-01 -1.430688e-01 -1.348010e-08 -8.944318e-01

## The last two parameters may be a little linked? Negative eigenvalue but this may in fact be zero? 

theta_init = theta_true + rnorm(6,0,sd=0.05)
p_opt<-optim(theta_init, RSS, y=data,hessian=T)
p_opt$par
theta_true
theta_init
#No obvious influence of theta_init on the results
p_opt$hessian
eigen(p_opt$hessian) ### clearly interesting and non zero
# $values
# [1] 3324.3545403 2000.8370768 1736.7462917    8.9679044    1.3913928    0.1586457
# $vectors
# [,1]          [,2]          [,3]          [,4]          [,5]          [,6]
# [1,]  7.751464e-01 -1.888266e-12  1.748739e-13 -7.603514e-12 -6.317817e-01  1.496517e-09
# [2,] -6.317817e-01  4.366499e-12 -2.381479e-12 -1.538362e-11 -7.751464e-01  1.837078e-09
# [3,] -4.196303e-12 -9.992907e-01  1.514892e-12  1.713220e-11  8.703747e-11  3.765710e-02
# [4,]  7.709085e-13  3.765710e-02  5.252822e-12  4.526253e-10  2.367879e-09  9.992907e-01
# [5,]  1.177811e-13 -1.195932e-12 -9.298648e-01  3.679013e-01 -4.540579e-12 -1.617067e-10
# [6,] -4.160374e-12  5.543679e-13  3.679013e-01  9.298648e-01 -1.619357e-11 -4.231339e-10

# $values
# [1] 3288.9438361 2000.8565182 1730.1249428    8.8419892    0.6927461    0.1618371
# 
# $vectors
# [,1]          [,2]          [,3]          [,4]          [,5]          [,6]
# [1,]  7.793633e-01 -2.303427e-12  3.370507e-12 -5.139319e-10 -6.265723e-01 -6.975312e-09
# [2,] -6.265723e-01 -4.417751e-12 -2.509825e-13 -6.387047e-10 -7.793633e-01 -8.678396e-09
# [3,]  9.214957e-13 -9.992859e-01 -2.733869e-11 -3.403855e-11 -4.158298e-10  3.778575e-02
# [4,] -1.376074e-12  3.778575e-02 -9.480527e-12 -9.308575e-10 -1.112639e-08  9.992859e-01
# [5,]  2.716018e-12  2.550093e-11 -9.299939e-01  3.675749e-01 -3.031194e-10  3.326169e-10
# [6,] -7.025175e-13 -8.833045e-12  3.675749e-01  9.299939e-01 -7.617032e-10  8.701316e-10

#### some percentage has negative eigenvalues but way less than before with the likelihood I think

solve(p_opt$hessian)

################################ Computation of Hessian #######################################################

### Computation at the true parameter values
### There's of course the question of whether it is better to compute the Hessian
### for the true param values or at the MLE (which might be slightly different on a finite time series). 

### In Gimenez et al. Animal Biodiversity and Conservation 27.1 (2004) They use the MLE.
### But for a theoretical work, what should we use? 

### Use other algorithm to compute the Hessian exactly at the true value
library("numDeriv")

theta_true  = c(rmax_V,1/K,sqrt(0.05),rmax_P,Q,sqrt(0.05),C,D,sqrt(0.05))
LL = function(theta){
  return(logLik(theta,data))
}
LL(theta_true)

hessian_thetaTrue=hessian(LL,theta_true)
eigen(hessian_thetaTrue)
## SD error
sqrt(diag(solve(hessian_thetaTrue)))
#[1] 0.243049234 0.288849849 0.005040121 0.028748093 0.728317339 0.004755960 0.043872686 0.110662255 0.005165725

#$values
# [1] 44294.175086 39378.550751 37848.716366 34111.484789 19976.028177 16277.744431    70.872841
# [8]     7.018634     1.882450
# 

### Removing last parameter for test with LL_FRwoutNoise
theta_true  = c(rmax_V,1/K,sqrt(0.05),rmax_P,Q,sqrt(0.05),C,D)

LL_FRwoutNoise = function(theta){
  return(logLik_FRwoutNoise(theta,data))
}

hessian_thetaTrue_FRwoutNoise=hessian(LL_FRwoutNoise,theta_true)
eigen(hessian_thetaTrue_FRwoutNoise)

# $values
# [1] 4.429418e+04 3.966980e+04 3.419146e+04 1.997603e+04 2.326595e+01 6.332812e+00 1.882450e+00 1.473978e-02
### Still a question of whether the last one is positive or not... 

### Using Viallefond et al.'s criteria
lambda1 = max(eigen(hessian_thetaTrue_FRwoutNoise)$values)
q = 8
eps = 10^(-9)
eigen(hessian_thetaTrue_FRwoutNoise)$values>q*lambda1*eps
# TRUE here

# Previously we had an error changing these $values
# [1]  4.429418e+04  3.966980e+04  3.419146e+04  1.997603e+04  2.326595e+01  6.332812e+00
# [7]  1.882450e+00  1.473978e-02 -4.179847e-23

eigen_FRwoutNoise=eigen(hessian_thetaTrue_FRwoutNoise)
#eigenvector associated with near-zero eigenvalue
eigen_FRwoutNoise$vectors[,8]
#largest link between C and D. Or is this correct? -> what's the new basis? 

################# Redo the same analysis with sum of squares #################

theta_true  = c(rmax_V,1/K,rmax_P,Q,C,D)
RSS_bis = function(theta){
  return(RSS(theta,data))
}
RSS_bis(theta_true)

hessian_RSS_thetaTrue=hessian(RSS_bis,theta_true)
eigen(hessian_RSS_thetaTrue)
# $values
# [1] 3411.9132542 2000.9252194 1631.3533267    7.1423266    0.7019349    0.1882889

# Now for the model with FR without noise

RSS_bis_FRwoutNoise = function(theta){
  return(RSS_FRwoutNoise(theta,data))
}
RSS_bis_FRwoutNoise(theta_true)

hessian_RSS_thetaTrue_FRwoutNoise=hessian(RSS_bis_FRwoutNoise,theta_true)
eigen(hessian_RSS_thetaTrue_FRwoutNoise)
# $values
# [1] 3.419685e+03 2.000925e+03 2.328096e+00 6.335312e-01 1.882889e-01 1.475940e-03
# 
# $vectors
# [,1]          [,2]          [,3]          [,4]          [,5]          [,6]
# [1,]  7.642358e-01 -5.185763e-15  5.166907e-01 -3.764703e-01 -1.339072e-11 -8.511460e-02
# [2,] -6.433480e-01  1.123952e-14  6.521352e-01 -3.883610e-01 -1.487730e-11 -9.999365e-02
# [3,]  9.900002e-15  9.992687e-01  2.487080e-13 -8.519244e-13  3.823704e-02 -7.097592e-13
# [4,] -3.011594e-15 -3.823704e-02  7.151991e-12 -2.188996e-11  9.992687e-01 -1.699767e-11
# [5,] -4.182289e-02 -1.298747e-15 -4.654567e-01 -8.059012e-01 -8.139938e-12  3.634889e-01
# [6,]  1.726023e-02  6.557340e-14  3.018285e-01  2.407676e-01  1.880481e-11  9.222975e-01

q = 6
lambda1 = max(eigen(hessian_RSS_thetaTrue_FRwoutNoise)$values)
eigen(hessian_RSS_thetaTrue_FRwoutNoise)$values>q*lambda1*eps
## Compare to singular value decomposition
eigen(hessian_RSS_thetaTrue_FRwoutNoise)$values
svd(hessian_RSS_thetaTrue_FRwoutNoise)$d

### Looks much better -- the last eigenvalue is fairly weak compared to the others though
### The last eigenvector suggests a linkage between the last two elements 
### Keep in mind that it is because these elements belong to the "null" space eigenvector that they indicate linkage. 
### These are C and D, which makes sense based on previous results. 

### Or can we use the different orders of magnitude for the eigenvalues here? 

### Clearly that's a different value from the one that we get at MLE, though there are some similarities
### in the spectra of the Hessian matrix at theta_true vs theta_hat. 

### NB there are some models in the literature that may not be identifiable 
### I'm curious about this one for instance, since fitted on only one compartment
### Ives et al. (2008). High-amplitude fluctuations and alternative dynamical states of midges in Lake Myvatn. Nature, 452(7183), 84.
#####################################################################################################################################

########################### Profiles of the likelihood near theta_true ##################################




