#### 15/03/2018 -- FBarraquand and OGimenez
#### Computation of the Hessian for the discrete-time Leslie-May model
#### New file with T=100


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
n.years<-100	# Long time series to look at structural identifiability
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
#### Warning the last parameter theta9 is not estimated!!
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
#[1]  1.7640017  0.7142281  0.2277440  0.5727717 11.9375808  0.2087442  2.5128467  1.0491252  0.2389437
theta_true
#[1]  2.0000000  1.0000000  0.2236068  0.5000000 10.0000000  0.2236068  2.5000000  1.0000000  0.2236068

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
# [1] 4544.6756366 4279.0130556 3817.9057096 3468.3549356 2274.8086575 1394.6774084    9.0787072    5.5516334    0.1141467

######## Estimation based on RSS ##################

### New restricted theta
theta_true  = c(rmax_V,1/K,rmax_P,Q,C,D)
p_opt<-optim(theta_true, RSS, y=data,hessian=T)
p_opt$par
theta_true
p_opt$hessian
eigen(p_opt$hessian) ### clearly interesting and non zero
### Old $values for t=1000
# [1] 3334.2367985 2000.8569587 1736.6593414    8.9622937    1.4761840    0.1603394
### New values for t=100
# $values
# [1] 454.80251632 198.25709097 159.49222726   1.04456291   0.59298376   0.01113483

# Let's try from somewhere slightly different
theta_init = theta_true + rnorm(6,0,sd=0.05)
p_opt<-optim(theta_init, RSS, y=data,hessian=T)
p_opt$par
theta_true
p_opt$hessian
eigen(p_opt$hessian) ### clearly interesting and non zero
# $values
# [1] 440.46957584 198.28217321 159.69830610   1.05047544   0.55995778   0.01463831
# 
# $vectors
# [,1]          [,2]          [,3]         [,4]          [,5]          [,6]
# [1,]  6.699400e-01  2.065958e-12 -6.270210e-13 7.096422e-10  7.424153e-01  1.868657e-10
# [2,] -7.424153e-01 -2.664349e-15 -2.594159e-12 6.408067e-10  6.699400e-01  1.690894e-10
# [3,]  1.372015e-12 -9.992881e-01 -1.563891e-11 2.579492e-11  1.103813e-11 -3.772525e-02
# [4,] -3.978080e-13  3.772525e-02 -3.483252e-12 7.140937e-10  2.517745e-10 -9.992881e-01
# [5,]  1.523714e-12  1.406952e-11 -9.345979e-01 3.557061e-01 -3.421681e-10  2.579773e-10
# [6,] -2.300136e-13 -6.599388e-12  3.557061e-01 9.345979e-01 -8.928331e-10  6.663765e-10


###################################################################################################
################### Now for the model where the FR is specified without the noise #################
###################################################################################################
# Revert back to full theta_true with 9 elements 
# Nope! It's actually 8 elements!!
theta_true  = c(rmax_V,1/K,sqrt(0.05),rmax_P,Q,sqrt(0.05),C,D) ### 8 elements only
p_opt<-optim(theta_true, logLik_FRwoutNoise, y=data,hessian=T)
p_opt$par ### Obviously a problem on the last param - but perhaps this could be done with constraints? 
theta_true
p_opt$hessian
eigen(p_opt$hessian)
# $values
# [1] 4441.0238334 4032.1190290 3288.9278103 2251.8416066   12.2374011    0.9329324    0.1419784   -0.1424231

# One negative eigenvalues -- does not seem to be doing well. 

######### Estimation starting away from MLE ########
theta_init = theta_true + rnorm(8,0,sd=0.05)
p_opt<-optim(theta_init, logLik_FRwoutNoise, y=data,hessian=T)
p_opt$par
# [1]  1.8712006  0.7963704  0.2284715  0.5822778 12.2094645  0.2085856  2.7754965  0.7688282
theta_true
theta_init
p_opt$hessian
eigen(p_opt$hessian)

### Going much better regarding estimated parameters! Still one of the eigenvalues is negative

# $values
# [1] 4562.5534412 3877.9858790 3784.6128492 2278.2132032    4.6680105    0.7002127    0.1077702   -0.1287456
# 
# $vectors
# [,1]          [,2]          [,3]          [,4]          [,5]          [,6]          [,7]          [,8]
# [1,] -2.001803e-12  6.989324e-01 -4.661571e-03  8.078860e-13  6.605280e-01 -2.439218e-01 -6.196734e-10  1.252063e-01
# [2,]  2.040020e-12 -7.138802e-01  3.769635e-03  3.544896e-13  6.657225e-01 -1.821585e-01 -7.111915e-10  1.182908e-01
# [3,] -2.286138e-12  5.953634e-03  9.999820e-01  5.189003e-13  6.246352e-04 -3.745008e-04  2.402706e-13 -2.906280e-05
# [4,] -2.323761e-04 -3.551782e-13 -5.454230e-13  9.993922e-01 -7.219949e-12 -8.607576e-11 -3.485950e-02 -3.066175e-10
# [5,]  7.199107e-05 -6.662195e-13 -8.062620e-13 -3.485948e-02 -1.784306e-10 -2.467440e-09 -9.993922e-01 -8.807948e-09
# [6,] -1.000000e+00 -2.892694e-12 -2.269019e-12 -2.347444e-04 -9.228824e-14 -6.137554e-13 -6.384680e-05 -9.480846e-13
# [7,]  6.795001e-13 -3.882151e-02  2.334626e-05 -3.542763e-14 -2.279896e-01 -9.079786e-01  5.361948e-09 -3.494116e-01
# [8,]  1.504452e-13  1.805621e-02 -1.899723e-04  8.307197e-13  2.617950e-01  2.879151e-01  7.359433e-09 -9.210007e-01


######## Estimation from RSS ################
theta_true  = c(rmax_V,1/K,rmax_P,Q,C,D)
p_opt<-optim(theta_true,RSS_FRwoutNoise, y=data,hessian=T)
p_opt$par
theta_true
p_opt$hessian
eigen(p_opt$hessian) ### clearly interesting and non zero? (inluding the last one?)
# $values
# [1] 3.113905e+02 1.982414e+02 1.388143e+00 1.100563e-01 4.245931e-02 9.481991e-03
# 
# $vectors
# [,1]          [,2]          [,3]          [,4]          [,5]          [,6]
# [1,]  7.972963e-01 -2.585307e-12  1.695954e-01  9.918239e-02 -5.707178e-01 -7.066790e-10
# [2,] -6.001156e-01  1.683015e-12  2.959810e-01  2.006327e-01 -7.155439e-01 -8.977019e-10
# [3,]  2.962597e-12  9.993910e-01 -1.913536e-12 -5.645251e-12  4.126922e-11 -3.489450e-02
# [4,]  9.314381e-13 -3.489450e-02 -1.304325e-11 -1.171858e-10  1.214692e-09 -9.993910e-01
# [5,] -5.681138e-02 -1.980770e-12 -2.549246e-01 -9.128736e-01 -3.137637e-01 -2.709746e-10
# [6,]  3.086059e-02  9.863041e-13  9.047906e-01 -3.414251e-01  2.526468e-01  3.352953e-10


## The last two parameters may be a little linked? This may in fact be zero? 

theta_init = theta_true + rnorm(6,0,sd=0.05)
p_opt<-optim(theta_init, RSS, y=data,hessian=T)
p_opt$par
theta_true
theta_init
#No obvious influence of theta_init on the results
p_opt$hessian
eigen(p_opt$hessian) ### clearly interesting and non zero
# $values
# [1] 4.437050e+02 1.982452e+02 1.592054e+02 1.035154e+00 5.756634e-01 9.893308e-03
# 
# $vectors
# [,1]          [,2]          [,3]          [,4]          [,5]          [,6]
# [1,]  6.674751e-01 -5.816374e-12  4.810703e-13 -3.709480e-10 -7.446321e-01  2.511515e-10
# [2,] -7.446321e-01  1.049480e-12  3.097723e-13 -3.339804e-10 -6.674751e-01  2.248009e-10
# [3,] -4.652314e-12 -9.993813e-01  8.683790e-12  5.741907e-12 -8.227046e-12 -3.517237e-02
# [4,]  4.073223e-13  3.517237e-02  1.163913e-13  2.147508e-10 -3.369836e-10 -9.993813e-01
# [5,] -3.046072e-13 -8.752388e-12 -9.346021e-01  3.556949e-01 -1.780708e-10  7.601617e-11
# [6,] -1.054619e-12  1.389139e-12  3.556949e-01  9.346021e-01 -4.662995e-10  2.009212e-10

###### Note -- we will probably need to evaluate this for many initial conditions #########

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
# $values
# [1] 4959.8096919 4809.0994550 3223.7316080 2816.5693745 1981.8791487 1617.1298279   11.2572042    0.1440341  -32.0313725
# The last eigenvalue is negative... 

## SD error
diag(solve(hessian_thetaTrue)) ### some negative variances... 
sqrt(diag(solve(hessian_thetaTrue)))

#### Previously we had --
#$values
# [1] 44294.175086 39378.550751 37848.716366 34111.484789 19976.028177 16277.744431    70.872841
# [8]     7.018634     1.882450
# 

############## Perhaps using the TRUE value only makes sense when T=1000!!! so that we are in near-perfect conditions #####
theta_true  = c(rmax_V,1/K,sqrt(0.05),rmax_P,Q,sqrt(0.05),C,D)

LL_FRwoutNoise = function(theta){
  return(logLik_FRwoutNoise(theta,data))
}

hessian_thetaTrue_FRwoutNoise=hessian(LL_FRwoutNoise,theta_true)
eigen(hessian_thetaTrue_FRwoutNoise)

######### Previously for T=1000 -- but we had an error 
# $values
# [1]  4.429418e+04  3.966980e+04  3.419146e+04  1.997603e+04  2.326595e+01  6.332812e+00
# [7]  1.882450e+00  1.473978e-02 -4.179847e-23
#### But we had a problem there as only 8 values are evaluated... 

################# Redo the same analysis with sum of squares #################

theta_true  = c(rmax_V,1/K,rmax_P,Q,C,D)
RSS_bis = function(theta){
  return(RSS(theta,data))
}
RSS_bis(theta_true)

hessian_RSS_thetaTrue=hessian(RSS_bis,theta_true)
eigen(hessian_RSS_thetaTrue)
# $values
# [1] 336.62510855 198.28492515 161.91041244   1.12619881   0.01454556  -3.19410640
# Now for the model with FR without noise

RSS_bis_FRwoutNoise = function(theta){
  return(RSS_FRwoutNoise(theta,data))
}
RSS_bis_FRwoutNoise(theta_true)

hessian_RSS_thetaTrue_FRwoutNoise=hessian(RSS_bis_FRwoutNoise,theta_true)
eigen(hessian_RSS_thetaTrue_FRwoutNoise)
# $values
# [1] 337.33536785 198.28492515   0.15246182   0.01454556  -0.02863044  -3.12787742
# 
# $vectors
# [,1]          [,2]          [,3]          [,4]         [,5]          [,6]
# [1,]  7.685884e-01  3.261043e-13  3.883697e-02 -3.965198e-12 2.333877e-02 -6.381370e-01
# [2,] -6.382035e-01  1.253299e-13 -1.843874e-02 -2.025826e-12 4.256147e-03 -7.696351e-01
# [3,]  1.671099e-13 -9.992812e-01  1.285985e-13  3.790849e-02 4.282178e-12 -3.805211e-13
# [4,] -5.012393e-16  3.790849e-02  2.040929e-12  9.992812e-01 1.093211e-10 -2.068036e-12
# [5,] -4.115850e-02  1.154317e-13  7.420921e-01 -7.463781e-11 6.687327e-01  2.004908e-02
# [6,]  1.655498e-02  6.678507e-14 -6.689178e-01 -7.992144e-11 7.431243e-01  6.407457e-03

### Similar to Likelihood results -- I assume that looking at the Hessian at theta_true only makes sense for very long time series, 
### i.e. cases where the MLE is actually very close to theta_true.  
#####################################################################################################################################

########################### Profiles of the likelihood near theta_true ##################################




