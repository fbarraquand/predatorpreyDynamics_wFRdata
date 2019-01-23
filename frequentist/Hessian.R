#### 15/03/2018 -- FBarraquand and OGimenez
#### Computation of the Hessian for the discrete-time Leslie-May model
### FB 16/01/2019 -- edited so that functional response timing match math model (no delays)


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
N_tplus1 = exp(y[t,1]) #useful for the functional response

######## Error that was previously there!! N_t/P_t instead P_t/N_t ###
### mu1 = y[t-1,1] + theta[1] - log(1+theta[2]*N_t) - y[t-1,3]*N_t/P_t
######################################################################
mu1 = y[t-1,1] + theta[1] - logprot(1+theta[2]*N_t) - y[t-1,3]*P_t/N_t #this is correct timing
mu2 = y[t-1,2] + theta[4] - logprot((1+theta[5]*P_t/N_t))
mu3 = (theta[7]*N_tplus1)/(theta[8] + N_tplus1) #easier to update all variables simultaneously, minimizes errors
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
    N_tplus1 = exp(y[t,1]) #useful for the functional response
    ############## Correction of error ###########################################
    ## mu1 = y[t-1,1] + theta[1] - log(1+theta[2]*N_t) - y[t-1,3]*N_t/P_t # error
    mu1 = y[t-1,1] + theta[1] - logprot(1+theta[2]*N_t) - y[t-1,3]*P_t/N_t # corrected
    mu2 = y[t-1,2] + theta[3] - logprot((1+theta[4]*P_t/N_t))
    mu3 = (theta[5]*N_tplus1)/(theta[6] + N_tplus1)
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
#p_opt<-optim(theta_init, logLik, y=data,method="BFGS",hessian=T,control=list(fnscale=-1)) ## we consider -LL, no need
p_opt<-optim(theta_init, logLik, y=data,method="BFGS",hessian=T)
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
# [1] 41478.246663 40362.227921 37312.035094 33666.368176 18680.108567 16862.134670    75.108129    12.326374     1.498224

lambda1 = max(eigen(hessian)$values)
9*lambda1*10^(-9) #way below lowest eigenvalue

# should be reproduced for several simulations later on. 

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
# [1] 3346.9992650 2000.8586776 1625.3684593    7.2470990    1.0993765    0.1614776
# 
# $vectors
# [,1]          [,2]          [,3]          [,4]          [,5]          [,6]
# [1,]  7.725414e-01 -1.637563e-13  1.812814e-12  1.123897e-10 -6.349644e-01 -3.712557e-09
# [2,] -6.349644e-01  2.852194e-16  1.203261e-12  1.383346e-10 -7.725414e-01 -4.519073e-09
# [3,]  7.567933e-14  9.992853e-01  7.809808e-12 -4.175188e-12 -2.211775e-10  3.780001e-02
# [4,] -1.350909e-12 -3.780001e-02 -5.240919e-12 -9.279350e-11 -5.844329e-09  9.992853e-01
# [5,]  9.624491e-13  7.691070e-12 -9.307156e-01  3.657437e-01  6.325099e-11  2.937278e-11
# [6,]  7.090222e-13 -2.308098e-12  3.657437e-01  9.307156e-01  1.666449e-10  8.825674e-11

# Let's try from somewhere slightly different
theta_init = theta_true + rnorm(6,0,sd=0.05)
p_opt<-optim(theta_init, RSS, y=data,hessian=T)
p_opt$par
theta_true
# > p_opt$par
# [1]  2.0612266  1.0739755  0.5063395 10.2958375  2.5142004  1.0192181
# > theta_true
# [1]  2.0  1.0  0.5 10.0  2.5  1.0
p_opt$hessian
eigen(p_opt$hessian) ### clearly interesting and non zero

# $values
# [1] 3251.8365314 2000.8568137 1625.2728372    7.2446891    1.0199015    0.1613022
# 
# $vectors
# [,1]          [,2]          [,3]          [,4]          [,5]         [,6]
# [1,]  7.837736e-01  2.290378e-12 -4.314895e-12 -6.754080e-10 -6.210466e-01 3.531534e-09
# [2,] -6.210466e-01  8.070546e-14 -1.016375e-12 -8.497219e-10 -7.837736e-01 4.456541e-09
# [3,] -1.751452e-12  9.992858e-01 -1.142597e-11 -1.450107e-11  2.163518e-10 3.778770e-02
# [4,] -1.372204e-13 -3.778770e-02  4.608813e-12 -3.555505e-10  5.682050e-09 9.992858e-01
# [5,] -1.956737e-12 -1.040251e-11 -9.307361e-01  3.656915e-01 -3.937041e-10 1.340139e-10
# [6,]  2.541700e-12  5.221157e-12  3.656915e-01  9.307361e-01 -1.011538e-09 3.296710e-10


###################################################################################################
################### Now for the model where the FR is specified without the noise #################
###################################################################################################
# Revert back to full theta_true with 9 elements
# Nope!! There are only 8 elements here, no sigma3!!

theta_true  = c(rmax_V,1/K,sqrt(0.05),rmax_P,Q,sqrt(0.05),C,D)
p_opt1<-optim(theta_true, logLik_FRwoutNoise, method="BFGS",y=data,hessian=T)
p_opt1$par
### BFGS 
# [1]  1.5515275  0.5574554  0.2226381  0.5064311 10.3011363  0.2314222 43.7378222 92.4917840
p_opt2<-optim(theta_true, logLik_FRwoutNoise,y=data,hessian=T)
p_opt2$par
### Nelder-Mead
# [1]  1.9890439  0.9734439  0.2223960  0.5013721 10.1733093  0.2319168  3.4129376  2.5978860
theta_true
# [1]  2.0000000  1.0000000  0.2236068  0.5000000 10.0000000  0.2236068  2.5000000  1.0000000

### Nelder-Mead is suprisingly better here

eigen(p_opt1$hessian) ## BFGS
# $values
# [1] 5.656082e+04 4.031409e+04 3.731116e+04 1.867995e+04 4.961284e+01 1.501037e+00 2.920489e-02 2.392860e-05
lambda1 = max(eigen(p_opt1$hessian)$values)
9*lambda1*10^(-9) #0.0005090474 #so the last eigenvalue is zero

eigen(p_opt2$hessian) ## Nelder-Mead
# $values
# $values
# [1]  4.073766e+04  3.691639e+04  3.521076e+04  1.860057e+04  1.470314e+01  3.518293e+00  1.532745e+00 -9.594615e-03
lambda1 = max(eigen(p_opt2$hessian)$values)
9*lambda1*10^(-9) ## 0.000366639 = 3e-4 so the last eigenvalue is not zero. 

######### Estimation starting away from MLE rather than exactly at MLE ########
######### Useful to investigate if/why Nelder-Mead does so good a job  ########

theta_init = theta_true + rnorm(8,0,sd=0.05) ## a little further away?
theta_init
theta_init = theta_true*(1+rnorm(8,0,sd=0.1)) ## not practical
### Crafting a new theta-init with reasonable a priori errors
theta_init = c(3,0.5,0.75,1,5,0.75,4,5)

p_opt1<-optim(theta_init, logLik_FRwoutNoise,method="BFGS", y=data,hessian=T)
p_opt1$par
theta_true
theta_init
# [1]  1.5718984  0.5730893  0.2226600  0.5067040 10.3083510  0.2314212 26.0655759 53.0488133
# > theta_true
# [1]  2.0000000  1.0000000  0.2236068  0.5000000 10.0000000  0.2236068  2.5000000  1.0000000
# > theta_init
# [1] 3.00 0.50 0.75 1.00 5.00 0.75 4.00 5.00

eigen(p_opt1$hessian) ### BFGS
# $values
# [1] 5.505953e+04 4.030562e+04 3.731199e+04 1.868010e+04 4.707488e+01 1.498589e+00 8.201987e-02 1.141953e-04
lambda1 = max(eigen(p_opt1$hessian)$values)
9*lambda1*10^(-9) ### 0.0005090474 -- close to zero. 
diag(solve(p_opt1$hessian)) ## positive values though
### Looks like BFGS performs worse than Nelder-Mead without the FR data (see below)

p_opt2<-optim(theta_init, logLik_FRwoutNoise, y=data,hessian=T)
p_opt2$par
eigen(p_opt2$hessian) ### Nelder-Mead
# $values
# [1] 42784.648389 39426.476918 30604.638451 18749.117486    46.054877     3.689130     2.639528    -0.111515
### Less good than wwhen starting at the maximum but still fairly reasonable

## Beware that we use -LL, so Hessian = observed FIM (Fisher Information Matrix)
## Should we loop over datasets to obtain the true, expected FIM?

######## Estimation from RSS ################
theta_true  = c(rmax_V,1/K,rmax_P,Q,C,D)
p_opt<-optim(theta_init,RSS_FRwoutNoise, y=data,hessian=T)
p_opt$par
theta_true
theta_init
# p_opt$par
# [1]   4.0791051   8.5632770   0.4753874   9.4701785   7.2905113   6.1694545 -12.9474382   2.6903906
# > theta_true
# [1]  2.0  1.0  0.5 10.0  2.5  1.0
# > theta_init
# [1] 3.00 0.50 0.75 1.00 5.00 0.75 4.00 5.00
### This is actually not very good. 

p_opt$hessian
eigen(p_opt$hessian) ### clearly interesting and non zero? (inluding the last one?)
# $values
# [1] 2026.57979293 2001.05550390    0.28303556    0.19744008    0.04618698    0.00000000    0.00000000   -0.01026505
## Negative eigenvalue but this may in fact be zero? Pb here. 

p_opt<-optim(theta_init, RSS, y=data,method="BFGS",hessian=T)
p_opt$par
theta_true
theta_init
# [1]  2.0298227  1.0349085  0.5066157 10.3060164  2.5146345  1.0202954  4.0000000  5.0000000
# > theta_true
# [1]  2.0  1.0  0.5 10.0  2.5  1.0
# > theta_init
# [1] 3.00 0.50 0.75 1.00 5.00 0.75 4.00 5.00
### This time it works better with BFGS it seems. Or does it? 
eigen(p_opt$hessian)
# $values
# [1] 3333.4379120 2000.8544900 1624.7208034    7.2369413    1.2214526    0.1606015    0.0000000    0.0000000
# 
# $vectors
# [,1]          [,2]          [,3]          [,4]          [,5]          [,6] [,7] [,8]
# [1,]  7.741023e-01  5.202697e-12 -1.676635e-12 -3.333724e-10 -6.330605e-01 -1.268594e-09    0    0
# [2,] -6.330605e-01 -2.100985e-12 -1.974447e-12 -4.104674e-10 -7.741023e-01 -1.549797e-09    0    0
# [3,]  5.387892e-12 -9.992864e-01  5.897838e-12  9.367340e-12 -7.731644e-11  3.777234e-02    0    0
# [4,]  7.036003e-13  3.777234e-02  2.565462e-12  3.079372e-10 -2.001306e-09  9.992864e-01    0    0
# [5,] -6.978705e-13 -6.225881e-12 -9.307386e-01  3.656853e-01 -1.909597e-10 -1.100637e-10    0    0
# [6,] -1.645111e-12  6.383782e-15  3.656853e-01  9.307386e-01 -4.931110e-10 -2.877529e-10    0    0
# [7,]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00    0    1
# [8,]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00    1    0

#################### Computation of Hessian at the true parameter values ######################

### There's of course the question of whether it is better to compute the Hessian 
### (aka observed Fisher Information Matrix, computed at the MLE)
### or the Hessian for the true parameter values? 

### If the time series is very long and the process ergodic *or* if we average over simulations, 
### we end up computing the expected Hessian = true FIM , for the true parameter values. 

### In Gimenez et al. Animal Biodiversity and Conservation 27.1 (2004) they the observed FIM. 
### But for a theoretical work, what should we use? 

### Use another algorithm to compute the Hessian exactly at the true value -- without optim() then
library("numDeriv")

theta_true  = c(rmax_V,1/K,sqrt(0.05),rmax_P,Q,sqrt(0.05),C,D,sqrt(0.05))
LL = function(theta){
  return(logLik(theta,data))
}
LL(theta_true)
# -223.51
hessian_thetaTrue=hessian(LL,theta_true)
eigen(hessian_thetaTrue)
# $values
# [1] 44294.175692 39378.566113 37850.532888 34120.067394 19976.032031 16291.064713    73.163052     6.801572     1.877072
## SD error
sqrt(diag(solve(hessian_thetaTrue)))
#[1] [1] 0.246939072 0.293385800 0.005040127 0.028788660 0.729359915 0.004756012 0.043196118 0.108919768 0.005146236

### Removing last parameter for test with LL_FRwoutNoise
theta_true  = c(rmax_V,1/K,sqrt(0.05),rmax_P,Q,sqrt(0.05),C,D)

LL_FRwoutNoise = function(theta){
  return(logLik_FRwoutNoise(theta,data))
}

hessian_thetaTrue_FRwoutNoise=hessian(LL_FRwoutNoise,theta_true)
eigen(hessian_thetaTrue_FRwoutNoise)
# $values
# [1]  4.429418e+04  3.969325e+04  3.419572e+04  1.997603e+04  1.825568e+01  3.140147e+00  1.877072e+00 -7.990887e-02
### Still a question of whether the last one is positive or not... 

### Using Viallefond et al.'s criteria
lambda1 = max(eigen(hessian_thetaTrue_FRwoutNoise)$values)
q = 8
eps = 10^(-9)
eigen(hessian_thetaTrue_FRwoutNoise)$values>q*lambda1*eps
# FALSE here

eigen_FRwoutNoise=eigen(hessian_thetaTrue_FRwoutNoise)
#eigenvector associated with near-zero eigenvalue
eigen_FRwoutNoise$vectors[,8]
# [1] -1.666351e-01 -1.963420e-01 -5.513071e-05 -1.499877e-13 -3.446900e-12 -2.555706e-14  3.464293e-01  9.020362e-01
#largest link between C and D. After that, between r and K. 

################# Redo the same analysis with sum of squares #################

theta_true  = c(rmax_V,1/K,rmax_P,Q,C,D)
RSS_bis = function(theta){
  return(RSS(theta,data))
}
RSS_bis(theta_true)
hessian_RSS_thetaTrue=hessian(RSS_bis,theta_true)
eigen(hessian_RSS_thetaTrue)
# $values
# [1] 3412.7730513 2000.9256614 1632.9214503    7.3169044    0.6802283    0.1877551

# Now for the model with FR without noise
RSS_bis_FRwoutNoise = function(theta){
  return(RSS_FRwoutNoise(theta,data))
}
RSS_bis_FRwoutNoise(theta_true)
hessian_RSS_thetaTrue_FRwoutNoise=hessian(RSS_bis_FRwoutNoise,theta_true)
eigen(hessian_RSS_thetaTrue_FRwoutNoise)
# $values
# [1]  3.419920e+03  2.000926e+03  1.825577e+00  3.140665e-01  1.877551e-01 -7.978823e-03
# 
# $vectors
# [,1]          [,2]          [,3]          [,4]          [,5]          [,6]
# [1,]  7.642811e-01 -8.053068e-14 -3.777157e-01  4.954033e-01  4.710235e-11  1.666759e-01
# [2,] -6.432938e-01 -4.846036e-14 -5.024479e-01  5.432783e-01  5.197209e-11  1.963870e-01
# [3,] -3.065831e-14 -9.992686e-01 -3.229404e-14  3.404351e-12 -3.823992e-02  2.723379e-13
# [4,]  9.353924e-16  3.823992e-02 -2.447700e-12  9.063696e-11 -9.992686e-01  7.462877e-12
# [5,] -4.183903e-02  2.311070e-15  6.490927e-01  6.759797e-01  5.713693e-11 -3.463806e-01
# [6,]  1.723317e-02 -1.167274e-14 -4.284339e-01 -4.975591e-02 -1.020057e-11 -9.020376e-01

q = 6
lambda1 = max(eigen(hessian_RSS_thetaTrue_FRwoutNoise)$values)
eigen(hessian_RSS_thetaTrue_FRwoutNoise)$values>q*lambda1*eps
## [1]  TRUE  TRUE  TRUE  TRUE  TRUE FALSE
## Compare to singular value decomposition just to check that it is similar
eigen(hessian_RSS_thetaTrue_FRwoutNoise)$values
svd(hessian_RSS_thetaTrue_FRwoutNoise)$d

### The last eigenvector suggests a linkage between the last two elements 
### Keep in mind that it is because these elements belong to the "null" space eigenvector that they indicate linkage. 
### These are C and D, which makes sense based on previous results. 

### Clearly that's a different value from the one that we get at MLE, though there are some similarities
### in the spectra of the Hessian matrix at theta_true vs theta_hat. 

### NB there are some models in the literature that may not be identifiable 
### I'm curious about this one for instance, since fitted on only one compartment
### Ives et al. (2008). High-amplitude fluctuations and alternative dynamical states of midges in Lake Myvatn. Nature, 452(7183), 84.
#####################################################################################################################################

########################### Profiles of the likelihood near theta_true ##################################

### See Likelihood_gaussianNoiseonFR.R file. 


