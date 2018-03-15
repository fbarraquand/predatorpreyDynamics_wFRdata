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
# [1]  3.959853e+04  3.758643e+04  3.441184e+04  1.873364e+04  1.861345e+01  4.247310e+00  1.453736e+00
# [8]  0.000000e+00 -4.407635e-02

# There is however the question of whether this matrix is positive definite. I'm not sure in fact. 
# The last eigenvalue is negative. On the other hand it is very small so we can consider it zero. 

######## Estimation based on RSS ##################

### New restricted theta
theta_true  = c(rmax_V,1/K,rmax_P,Q,C,D)
p_opt<-optim(theta_true, RSS, y=data,hessian=T)
p_opt$par
theta_true
p_opt$hessian
eigen(p_opt$hessian) ### clearly interesting and non zero
$values
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
theta_true  = c(rmax_V,1/K,sqrt(0.05),rmax_P,Q,sqrt(0.05),C,D,sqrt(0.05))
p_opt<-optim(theta_true, logLik_FRwoutNoise, y=data,hessian=T)
p_opt$par
theta_true
p_opt$hessian
eigen(p_opt$hessian)

######### Estimation starting away from MLE ########
theta_init = theta_true + rnorm(9,0,sd=0.05)
p_opt<-optim(theta_init, logLik_FRwoutNoise, y=data,hessian=T)
p_opt$par
theta_true
theta_init
p_opt$hessian
eigen(p_opt$hessian)
# $values
# [1] 4.014639e+04 3.728309e+04 3.512972e+04 1.867364e+04 2.195517e+01 5.424323e+00 1.465820e+00 4.758790e-02
# [9] 0.000000e+00


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
#$values
# [1] 44294.175086 39378.550751 37848.716366 34111.484789 19976.028177 16277.744431    70.872841
# [8]     7.018634     1.882450
# 

LL_FRwoutNoise = function(theta){
  return(logLik_FRwoutNoise(theta,data))
}

hessian_thetaTrue_FRwoutNoise=hessian(LL_FRwoutNoise,theta_true)
eigen(hessian_thetaTrue_FRwoutNoise)

# $values
# [1]  4.429418e+04  3.966980e+04  3.419146e+04  1.997603e+04  2.326595e+01  6.332812e+00
# [7]  1.882450e+00  1.473978e-02 -4.179847e-23

eigen_FRwoutNoise=eigen(hessian_thetaTrue_FRwoutNoise)
#eigenvector associated with near-zero eigenvalue
eigen_FRwoutNoise$vectors[,9]
### Not like Gimenez et al. 2004 -- only the last value (sigma3) has a large eigenvector
### does this mean it canot be estimated?? Though perhaps the penultimate as well. 

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

### Looks much better -- the last eigenvalue is fairly weak compared to the others though
### The last eigenvector suggests a linkage between the last two elements 
### Keep in mind that because these elements belong to the null space eigenvector, they indicate linkage. 
### These are C and D, which makes sense based on previous results. 

### Or can we use the different orders of magnitude for the eigenvalues here? 

### Clearly that's a different value from the one that we get at MLE, though there are some similarities
### in the spectra of the Hessian matrix at theta_true vs theta_hat. 

### NB there are some models in the literature that may not be identifiable 
### I'm curious about this one for instance, since fitted on only one compartment
### Ives et al. (2008). High-amplitude fluctuations and alternative dynamical states of midges in Lake Myvatn. Nature, 452(7183), 84.
#####################################################################################################################################

########################### Profiles of the likelihood near theta_true ##################################



