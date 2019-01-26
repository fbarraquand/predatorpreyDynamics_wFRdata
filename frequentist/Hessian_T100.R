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
n.years<-100  	# Long time series to look at structural identifiability
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
plot(N,FR)
plot(log(N),log(P))

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

theta_true  = c(rmax_V,1/K,sqrt(0.05),rmax_P,Q,sqrt(0.05),C,D,sqrt(0.05))
######### Estimation starting away from MLE ########
theta_init = theta_true + rnorm(9,0,sd=0.05)
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

eigen(hessian)$values
# [1] 4542.8536605 4286.2765998 3817.9241132 3484.7171714 2273.9059709 1321.4767269    5.8780348
# [8]    5.5188393    0.1139454
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

######## Estimation based on RSS - removed here ##################

###################################################################################################
################### Now for the model where the FR is specified without the noise #################
###################################################################################################
# there are only 8 elements here, no sigma3!!

theta_true  = c(rmax_V,1/K,sqrt(0.05),rmax_P,Q,sqrt(0.05),C,D)
######### Estimation starting away from MLE  ########

theta_init = theta_true + rnorm(8,0,sd=0.05) ## a little further away?
theta_init
theta_init = theta_true*(1+rnorm(8,0,sd=0.1)) ## not practical
### Crafting a new theta-init with reasonable a priori errors
theta_init = c(3,0.5,0.75,1,5,0.75,4,5)

p_opt1<-optim(theta_init, logLik_FRwoutNoise,method="BFGS", y=data,hessian=T)
p_opt1$par
theta_true
theta_init
# > p_opt1$par
# [1]  2.2760846  1.3026393  0.2274158  0.5701208 11.8620933  0.2087870  1.9523082 -0.5548188
# > theta_true
# [1]  2.0000000  1.0000000  0.2236068  0.5000000 10.0000000  0.2236068  2.5000000  1.0000000
# > theta_init
# [1] 3.00 0.50 0.75 1.00 5.00 0.75 4.00 5.00
eigen(p_opt1$hessian)$values ### BFGS
# [1] 4542.7949743 3828.9389375 2791.5840144 2273.8941818   20.1873069    1.2145190    0.2970695
# [8]    0.1139395
lambda1 = max(eigen(p_opt1$hessian)$values)
9*lambda1*10^(-9) ### 4.088515e-05 we are above this
diag(solve(p_opt1$hessian)) ## errors on parameters

p_opt2<-optim(theta_init, logLik_FRwoutNoise, y=data,hessian=T)
p_opt2$par ## much less good 
eigen(p_opt2$hessian)$values ### Nelder-Mead
# [1] 3872.75548876 3421.06107496 2219.03709015 2022.71574588    1.31418173    0.30652243
# [7]    0.20553044   -0.01583548
# OK last one is zero. 

## Beware that we use -LL, so Hessian = observed FIM (Fisher Information Matrix)
## Should we loop over datasets to obtain the true, expected FIM?
## Does that make sense to compute identifiability at a wrong parameter value????

######## Estimation from RSS ################
theta_true  = c(rmax_V,1/K,rmax_P,Q,C,D)
p_opt<-optim(theta_init,RSS_FRwoutNoise, y=data,hessian=T)
p_opt$par
theta_true
theta_init
# [1] 1.5919112 0.6291094 0.4703591 9.1895996 1.5833255 3.3719153 7.2144505 8.1993919
# > theta_true
# [1]  2.0  1.0  0.5 10.0  2.5  1.0
# > theta_init
# [1] 3.00 0.50 0.75 1.00 5.00 0.75 4.00 5.00
### not very good. 

p_opt<-optim(theta_init, RSS_FRwoutNoise, y=data,method="BFGS",hessian=T)
p_opt$par
theta_true
theta_init
# > p_opt$par
# [1]  1.2988370  0.6316623  0.5236089 10.6387390  4.6682049 -0.1580577  4.0000000  5.0000000
# > theta_true
# [1]  2.00  1.00  0.50 10.00 15.00  0.25
# > theta_init
# [1] 3.00 0.50 0.75 1.00 5.00 0.75 4.00 5.00
eigen(p_opt$hessian)$values
# [1] 2.887165e+03 3.588719e+02 1.982219e+02 6.386789e+00 1.494630e+00 6.157179e-02 0.000000e+00
# [8] 0.000000e+00 -- might a problem here, algo has not moved. 

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
# [1] 4960.0244666 4808.9592996 3223.7316917 2816.9902724 1981.8798602 1624.9940079    8.5864833
# [8]    0.1417606  -32.0562507
## aha, problem. Which is logical since we don't have the expected FIM here. 

### Removing last parameter for test with LL_FRwoutNoise
theta_true  = c(rmax_V,1/K,sqrt(0.05),rmax_P,Q,sqrt(0.05),C,D)

LL_FRwoutNoise = function(theta){
  return(logLik_FRwoutNoise(theta,data))
}

hessian_thetaTrue_FRwoutNoise=hessian(LL_FRwoutNoise,theta_true)
eigen(hessian_thetaTrue_FRwoutNoise)
# $values
# [1] 4969.0842204 3223.7316917 2818.7286510 1981.8798602    1.6168326    0.1417606   -0.3183379
# [8]  -32.1364430

### We need to average the FIM over many many repeats or compute only with very long TS!!


