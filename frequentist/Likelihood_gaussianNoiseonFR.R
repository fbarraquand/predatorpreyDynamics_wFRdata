### FB 13/11/2017 - predator-prey model with noisy functional response data
### Writes down the likelihood of the model under the assumption of Gaussian noise on the FR
### Designed to be the simplest form noise tested first
### FB 09/01/2019 - better likelihood profiles for pairs of potentially correlated params
### FB 17/01/2019 -- edited so that functional response timing match math model (no delays)
### FB 30/01/2019 -- added simulations without attack rate data 

rm(list=ls())
graphics.off()

### Model used (use stored simulated data later on, just tryouts for now)
n.years<-1000  	# Long time series to look at structural identifiability
N1<-1			# Initial pop size
P1<-0.1
K<-1			# threshold dd 
beta<-1			# density-dependence exponent
rmax_V<-2			# Max AVERAGE growth rate (thus not a true max...)
rmax_P<-0.5
sigma2.proc<-0.05		# worked well with 0.005 (Note: Peretti et al. value. 0.005)

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
  FR[t+1]<-(C*N[t+1]/(D+N[t+1])) + FRnoise[t+1]
}
## Plotting time series of abundances and FR
par(mfrow=c(2,2))
plot(1:n.years,N,type="b")
plot(1:n.years,P,type="b")
#curve(dbeta(x,a,b),from=0, to=1)
plot(N,FR)

### Produce dataset to fit
data = cbind(log(N),log(P),FR) 

##################### Helper function #####################################################
### By OG 
logprot <- function(v){
  eps <- 2.2204e-016
  u <- log(eps) * (1+vector(length=length(v)))
  index <- (v>eps)
  u[index] <- log(v[index])
  u
}

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


##################################### Now optimize the LL ###########################################

### Here, to contrast with the Hessian.R code, we start with more "difficult" starting values since 
### we aim at identifying "global" properties of the likelihood surface. 

# theta_start = c(runif(1,0.5,2),runif(1,0.001,5),runif(1,0.1,1),runif(1,0.1,2),runif(1,0.1,5),runif(1,0.1,2),runif(1,0.01,0.1),runif(1,10,100),runif(1,0.1,5))
# Put high values for D not C -- previously the reverse was done for bayesian est. (error?)
theta_start = c(runif(1,0.5,2),runif(1,0.001,5),runif(1,0.05,1),runif(1,0.1,1),runif(1,0.1,5),runif(1,0.05,1),runif(1,0.01,0.1),runif(1,1,5),runif(1,0.05,1))
theta_true  = c(rmax_V,1/K,sqrt(0.05),rmax_P,Q,sqrt(0.05),C,D,sqrt(0.05))


p_opt<-optim(theta_start, logLik, y=data,method="BFGS",hessian=T)
### warnings initially: log(1 + theta[2] * N_t) : NaNs produced --> use of logprot to avoid this. 
warnings() ## no big deal (see below, it can handle it )
### 1: In dnorm(y[t, 3], mu3, theta[9], log = T) : NaNs produced
### See also comments in Ben Bolker's EDMB book. 
### Likely due to some negative values of sigmas during iterations.  

### does it converge?
p_opt$convergence
p_opt$fn
p_opt$par
theta_true
theta_start
# > p_opt$convergence
# [1] 1 ### max iteration limit reached. 
# > p_opt$fn
# NULL
# > p_opt$par
# [1]  2.0821862  1.0979641  0.2223287  0.5085918 10.3523280  0.2312247  2.5154473  1.0211207  0.2195082
# > theta_true
# [1]  2.0000000  1.0000000  0.2236068  0.5000000 10.0000000  0.2236068  2.5000000  1.0000000  0.2236068
# > theta_start
# [1] 1.48343038 4.01204424 0.42927949 0.61820082 3.44240320 0.76254955 0.04287112 3.92265334 0.27997384
## Estimation quality can depend on initial values, especially sigmas (difficult if chosen too large)

p_opt<-optim(theta_start, logLik, y=data,method="BFGS",hessian=T,control=list(maxit=1000)) ## increasing iterations (x10)
# > p_opt$convergence
# [1] 0
# > p_opt$fn
# NULL
# > p_opt$par
# [1]  2.0357938  1.0422482  0.2224955  0.5065761 10.3051520  0.2314203  2.5146110  1.0202603  0.2194918
# > theta_true
# [1]  2.0000000  1.0000000  0.2236068  0.5000000 10.0000000  0.2236068  2.5000000  1.0000000  0.2236068
# > theta_start
# [1] 1.48343038 4.01204424 0.42927949 0.61820082 3.44240320 0.76254955 0.04287112 3.92265334 0.27997384

### Another theta_start
# > p_opt$convergence
# [1] 0
# > p_opt$fn
# NULL
# > p_opt$par
# [1]  2.0295409  1.0345492  0.2225045  0.5067380 10.3093475  0.2314224  2.5145841  1.0201612  0.2194936
# > theta_true
# [1]  2.0000000  1.0000000  0.2236068  0.5000000 10.0000000  0.2236068  2.5000000  1.0000000  0.2236068
# > theta_start
# [1] 1.52053301 0.35408353 0.71520718 0.94605200 4.51362848 0.27847591 0.04492732 2.55811540 0.68880964

## Still warnings() with this version of the code
# Warning messages:
# 1: In dnorm(y[t, 2], mu2, theta[6], log = T) : NaNs produced
# 2: In dnorm(y[t, 2], mu2, theta[6], log = T) : NaNs produced
## May be due to use try out of negative values for theta[3], theta[6] and theta[9] -- always these ones. 


### Let's check with L-BFGS-B -- all our parameters are non-negative
### The function fn can return NAs or Inf in optim() but not for L-BFGS-B
p_opt<-optim(theta_start, logLik, y=data,method="L-BFGS-B",lower=0.01,hessian=T)
p_opt$convergence
p_opt$fn
p_opt$par
theta_true
theta_start
# > p_opt$convergence
# [1] 1
# > p_opt$fn
# NULL
# > p_opt$par
# [1]  2.0298601  1.0349542  0.2224989  0.5068188 10.3113274  0.2314191  2.5145319  1.0200425  0.2194903
# > theta_true
# [1]  2.0000000  1.0000000  0.2236068  0.5000000 10.0000000  0.2236068  2.5000000  1.0000000  0.2236068
# > theta_start
# [1] 1.52053301 0.35408353 0.71520718 0.94605200 4.51362848 0.27847591 0.04492732 2.55811540 0.68880964

# Simpler Nelder-Mead (+ if functions are non-differentiable, otherwise BFGS may perform better)
p_opt<-optim(theta_start, logLik, y=data,hessian=T)
p_opt$convergence
p_opt$fn
p_opt$par
theta_true
theta_start
# > p_opt$convergence
# [1] 1
# > p_opt$fn
# NULL
# > p_opt$par
# [1] 0.8263124 0.1850046 0.2279861 0.1643470 3.6805673 0.1777670 1.5887547 1.9490441 1.0947071
# > theta_true
# [1]  2.0000000  1.0000000  0.2236068  0.5000000 10.0000000  0.2236068  2.5000000  1.0000000  0.2236068
# > theta_start
# [1] 1.52053301 0.35408353 0.71520718 0.94605200 4.51362848 0.27847591 0.04492732 2.55811540 0.68880964
## we should probably stick to BFGS 

library(optimx) ## check with different package
p_opt<-optimx(theta_start, logLik, y=data,hessian=T)
#slow - still same warnings -- outputs both algo simultaneously though, which is cool
p_opt
theta_true
theta_start
# p_opt
# p1        p2        p3       p4        p5        p6       p7       p8        p9     value fevals
# Nelder-Mead 0.8263124 0.1850046 0.2279861 0.164347  3.680567 0.1777670 1.588755 1.949044 1.0947071 1571.0164    502
# BFGS        2.0295409 1.0345492 0.2225045 0.506738 10.309347 0.2314224 2.514584 1.020161 0.2194936 -225.7685    290
# gevals niter convcode  kkt1  kkt2  xtime
# Nelder-Mead     NA    NA        1 FALSE FALSE 11.956
# BFGS            56    NA        0 FALSE  TRUE 37.148
# > theta_true
# [1]  2.0000000  1.0000000  0.2236068  0.5000000 10.0000000  0.2236068  2.5000000  1.0000000  0.2236068
# > theta_start
# [1] 1.52053301 0.35408353 0.71520718 0.94605200 4.51362848 0.27847591 0.04492732 2.55811540 0.68880964

### Let's try something extreme - starting with the true value
p_opt<-optim(theta_true, logLik, y=data,method="BFGS",hessian=T)
p_opt$convergence
p_opt$fn
p_opt$par
theta_true
# > p_opt$convergence
# [1] 0
# > p_opt$fn
# NULL
# > p_opt$par
# [1]  2.0297892  1.0348692  0.2225044  0.5067429 10.3093893  0.2314211  2.5146378  1.0203027  0.2194915
# > theta_true
# [1]  2.0000000  1.0000000  0.2236068  0.5000000 10.0000000  0.2236068  2.5000000  1.0000000  0.2236068

p_opt<-optim(theta_true, logLik, y=data,hessian=T)
p_opt$convergence
p_opt$fn
p_opt$par
theta_true
# > p_opt$convergence
# [1] 1
# > p_opt$fn
# NULL
# > p_opt$par
# [1]  1.9897758  0.9869458  0.2224962  0.5053926 10.2760870  0.2314583  2.5144643  1.0197673  0.2194978
# > theta_true
# [1]  2.0000000  1.0000000  0.2236068  0.5000000 10.0000000  0.2236068  2.5000000  1.0000000  0.2236068
### In this straightforward case, Nelder-Mead works similarly to BGFS, in other cases BFGS is better.  
### (Note: Nelder-Mead may work better for some different initial conditions for theta)

############ Create likelihood profiles ########################################

######### Reminder ############
### How does contour() plot?
A=matrix(c(1,2,3,4,1,2,3,4,1,2,3,4),nrow=3,byrow=T)
A
contour(c(1,2,3),c(1,2,3,4),A)
# here we have values going up -- the first vector is in rows and the second in columns. 
# this corresponds to C [row i index] being as the x-axis (therefore as column in the image)
# and to D [column j index] being as the y-axis (therefore as row in the image)
######## All is well! #########

par(mfrow=c(2,1))
# basic informal check for C and D
logLik(theta_true,data)

DeltaC = 2 # C is +/- 2
DeltaD = 0.9 #D is +/-0.9
niter = 50
C_new=D_new=rep(0,niter)
llbis=matrix(0,nrow=niter,ncol=niter)
for (i in 1:niter){
  for (j in 1:niter){
    C_new[i] =2*DeltaC*i/niter-DeltaC
    D_new[j] = 2*DeltaD*j/niter-DeltaD
    theta_new=theta_true+c(0,0,0,0,0,0,C_new[i],D_new[j] ,0)
    llbis[i,j]=logLik(theta_new,data)
  }
}
contour(theta_true[7]+C_new,theta_true[8]+D_new,llbis,nlevels=20,xlab="C",ylab="D")
hist(llbis) ## what values are in there
min(llbis)
custom_levels=c(40000,20000,10000,5000,1000,500,100,50,0,-50,-100,-150,-200,-210,-220)
contour(theta_true[7]+C_new,theta_true[8]+D_new,llbis,nlevels=20,xlab="C",ylab="D")
contour(theta_true[7]+C_new,theta_true[8]+D_new,llbis,levels=custom_levels,xlab="C",ylab="D")

### Zooming in
interval = (niter/2-10):(niter/2+10)
custom_levels=c(1000,500,100,50,0,-50,-100,-150,-200,-210,-215,-220)
contour(theta_true[7]+C_new[interval],theta_true[8]+D_new[interval],llbis[interval,interval],levels=custom_levels,xlab="C",ylab="D")

### Let's make that zoom much more precise - small wiggles are likely discretization artefacts
DeltaC = 0.25 # C is +/- 2
DeltaD = 0.5 #D is +/-0.9
niter = 50
C_new=D_new=rep(0,niter)
llbis=matrix(0,nrow=niter,ncol=niter)
for (i in 1:niter){
  for (j in 1:niter){
    C_new[i] =2*DeltaC*i/niter-DeltaC
    D_new[j] = 2*DeltaD*j/niter-DeltaD
    theta_new=theta_true+c(0,0,0,0,0,0,C_new[i],D_new[j] ,0)
    llbis[i,j]=logLik(theta_new,data)
  }
}
hist(llbis) ## what values are in there
min(llbis)
custom_levels=c(500,100,0,-50,-200,-210,-215,-220,-225)
contour(theta_true[7]+C_new,theta_true[8]+D_new,llbis,nlevels=20,xlab="C",ylab="D")
contour(theta_true[7]+C_new,theta_true[8]+D_new,llbis,levels=custom_levels,xlab="C",ylab="D")
### Here we see quite well the optimum even though there is a correlation. 
### The reason is that the slope decreases gradually

### Try to vizualise this using a 3D surface plot
library(rgl)
interval = (niter/2-10):(niter/2+10)
custom_levels=c(1000,500,100,50,0,-50,-100,-150,-200,-210,-215,-220)
persp3d(theta_true[7]+C_new[interval],theta_true[8]+D_new[interval],-llbis[interval,interval], col="skyblue")
### Maximum not visible on the surface. 

################# Now with (r,gamma=1/K)
logLik(theta_true,data)
niter = 50
r_new=g_new=rep(0,niter)
llbis=matrix(0,nrow=niter,ncol=niter)
Deltar = 1.5 # rmax_V is +/- 1.5
Deltag = 0.9 #1/K is +/- 0.9
for (i in 1:niter){
  for (j in 1:niter){
    r_new[i] =2*Deltar*i/niter-Deltar
    g_new[j] = 2*Deltag*j/niter-Deltag
    theta_new=theta_true+c(r_new[i],g_new[j],0,0,0,0,0,0)
    llbis[i,j]=logLik(theta_new,data)
  }
}
contour(theta_true[1]+r_new,theta_true[2]+g_new,llbis,nlevels=50,xlab="r",ylab="1/K")

#changing the representation of the levels
hist(llbis) ## what values are in there
custom_levels=c(80000,60000,40000,20000,5000,1000,500,100,50,0,-10,-20-30)
contour(theta_true[1]+r_new,theta_true[2]+g_new,llbis,levels=custom_levels,xlab="r",ylab="1/K")

# Surface plotting on a restricted area
interval = (niter/2-10):(niter/2+10)
persp3d(theta_true[1]+r_new[interval],theta_true[2]+g_new[interval],-llbis[interval,interval], col="skyblue")


############################################# More analyses ##########################################
# Let's assume we know the max growth rates (most realistic assumption)
# Plot these now for many 1/K values - let's imagine we get this wrong -> is there still a well-defined MLE for (C,D)?
gamma = c(0.5,1,1.5,2)
negll=matrix(0,nrow=50,ncol=51)
custom_levels=c(40000,20000,10000,5000,1000,500,100,0,-50,-200,-210,-215,-220,-225)
par(mfrow=c(2,2))
for (kg in 1:length(gamma)){
  D_FR = seq(0.01,2.02,by=0.04)
    for (kd in 1:length(D_FR)){
      C_FR = seq(0.05,5,0.1) 
      for (kc in 1:length(C_FR)){
       theta = c(rmax_V,gamma[kg],sqrt(0.05),rmax_P,Q,sqrt(0.05),C_FR[kc],D_FR[kd],sqrt(0.05))
       negll[kc,kd]=logLik(theta,data)
      }
    }
  custom_levels=quantile(negll,probs=c(0.025,0.05,0.075,0.1,0.25,0.5,0.75),na.rm=T)
  contour(C_FR,D_FR,negll,levels=custom_levels,xlab="C",ylab="D", main=c("1/K = ",toString(gamma[kg])))
}
# computation time for 100 x 101 grid: 14:57 to 15:08
plot(C_FR,negll[,length(D_FR)/2])
plot(D_FR,negll[length(C_FR)/2,])

########### Other profiles #####################################
# We have (r_V, 1/K) and (C,D) - we also need (r_P,Q)
################# Now with (r,gamma=1/K)

niter = 50
rP_new=q_new=rep(0,niter)
llbis=matrix(0,nrow=niter,ncol=niter)
Deltar = 0.3 # rmax_P is +/- 0.3
Deltaq = 2 # Q is +/- 2
for (i in 1:niter){
  for (j in 1:niter){
    rP_new[i] =2*Deltar*i/niter-Deltar
    q_new[j] = 2*Deltaq*j/niter-Deltaq
    theta_new=theta_true+c(0,0,0,rP_new[i],q_new[j],0,0,0,0)
    llbis[i,j]=logLik(theta_new,data)
  }
}
contour(theta_true[4]+rP_new,theta_true[5]+q_new,llbis,nlevels=50,xlab="rP",ylab="Q")

#changing the representation of the levels
hist(llbis) ## what values are in there
min(llbis)
custom_levels=quantile(llbis,probs=c(0.025,0.05,0.075,0.1,0.25,0.5,0.75),na.rm=T)
contour(theta_true[4]+rP_new,theta_true[5]+q_new,llbis,levels=custom_levels,xlab="rP",ylab="Q")

# Surface plotting on a restricted area
interval = (niter/2-10):(niter/2+10)
persp3d(theta_true[4]+rP_new[interval],theta_true[5]+q_new[interval],-llbis[interval,interval], col="skyblue")

# Zoom
custom_levels=quantile(llbis,probs=c(0.001,0.005,0.025,0.075,0.1,0.25,0.5),na.rm=T)
contour(theta_true[4]+rP_new[interval],theta_true[5]+q_new[interval],llbis[interval,interval],levels=custom_levels,xlab="rP",ylab="Q")
###############################################################################   

### Working directly with the sum of squares. 

# basic check 
RSS(theta_true,data)

DeltaC = 2 # C is +/- 2
DeltaD = 0.9 #D is +/-0.9
niter = 50
rssbis=matrix(0,nrow=niter,ncol=niter)
C_new=D_new=rep(0,niter)

for (i in 1:niter){
  for (j in 1:niter){
    C_new[i] =2*DeltaC*i/niter-DeltaC
    D_new[j] = 2*DeltaD*j/niter-DeltaD
    theta_new=theta_true+c(0,0,0,0, C_new[i],D_new[j])
    rssbis[i,j]=RSS(theta_new,data)
  }
}
hist(rssbis)
custom_levels = quantile(rssbis,probs=c(0.025,0.05,0.075,0.1,0.25,0.5,0.75),na.rm=T)
contour(theta_true[5]+C_new,theta_true[6]+D_new,rssbis,levels=custom_levels,xlab="C",ylab="D")
## Very similar - which is to be expected!

### Do we see a likelihood ridge?
# do that with (r,gamma)
RSS(theta_true,data)
Deltar = 1.5 # rmax_V is +/- 1.5
Deltag = 0.9 #1/K is +/- 0.9

rssbis=matrix(0,nrow=niter,ncol=niter)
r_new=g_new=rep(0,niter)
for (i in 1:niter){
  for (j in 1:niter){
    r_new[i] =2*Deltar*i/niter-Deltar
    g_new[j] = 2*Deltag*j/niter-Deltag
    theta_new=theta_true+c(r_new[i],g_new[j],0,0,0,0)
    rssbis[i,j]=RSS(theta_new,data)
  }
}
custom_levels = quantile(rssbis,probs=c(0.025,0.05,0.075,0.1,0.25,0.5,0.75),na.rm=T)
contour(theta_true[1]+r_new,theta_true[2]+g_new,rssbis,levels=custom_levels,xlab="r",ylab="1/K")
### Interestingly, this does seem to provide quite a different profile from the likelihood, 
### although qualitatively the results are the same -- smthg to explore? 

# Surface plotting on a restricted area
interval = (niter/2-10):(niter/2+10)
persp3d(theta_true[1]+r_new[interval],theta_true[2]+g_new[interval],rssbis[interval,interval], col="skyblue")


### Zooming in -- on r especially
Deltar = 0.2 # rmax_V is +/- 0.5
Deltag = 0.2 #1/K is +/- 0.5
rssbis=matrix(0,nrow=niter,ncol=niter)
r_new=g_new=rep(0,niter)
for (i in 1:niter){
  for (j in 1:niter){
    r_new[i] =2*Deltar*i/niter-Deltar
    g_new[j] = 2*Deltag*j/niter-Deltag
    theta_new=theta_true+c(r_new[i],g_new[j],0,0,0,0)
    rssbis[i,j]=RSS(theta_new,data)
  }
}
custom_levels = quantile(rssbis,probs=c(0.025,0.05,0.1,0.5),na.rm=T)
contour(theta_true[1]+r_new,theta_true[2]+g_new,rssbis,levels=custom_levels,xlab="r",ylab="1/K")


################################# Without the FR data ############################
par(mfrow=c(2,2))

### Looking only at (C,D) since this is the part difficult to identify. 
theta_true  = c(rmax_V,1/K,sqrt(0.05),rmax_P,Q,sqrt(0.05),C,D)

DeltaC = 2.5 # C is +/- 2
DeltaD = 1.0 #D is +/-0.9
niter = 50
C_new=D_new=rep(0,niter)
llbis=matrix(0,nrow=niter,ncol=niter)
for (i in 1:niter){
  for (j in 1:niter){
    C_new[i] =2*DeltaC*i/niter-DeltaC
    D_new[j] = 2*DeltaD*j/niter-DeltaD
    theta_new=theta_true+c(0,0,0,0,0,0,C_new[i],D_new[j] ,0)
    llbis[i,j]=logLik_FRwoutNoise(theta_new,data)
  }
}
contour(theta_true[7]+C_new,theta_true[8]+D_new,llbis,nlevels=20,xlab="C",ylab="D")
custom_levels=quantile(llbis,probs=c(0.025,0.05,0.075,0.1,0.25,0.5,0.75),na.rm=T)
contour(theta_true[7]+C_new,theta_true[8]+D_new,llbis,levels=custom_levels,xlab="C",ylab="D")


persp3d(theta_true[7]+C_new,theta_true[8],llbis, col="skyblue")

### Let's make that zoom much more precise - around the true values
DeltaC = 0.25 # C is +/- 2
DeltaD = 0.5 #D is +/-0.9
niter = 50
C_new=D_new=rep(0,niter)
llbis=matrix(0,nrow=niter,ncol=niter)
for (i in 1:niter){
  for (j in 1:niter){
    C_new[i] =2*DeltaC*i/niter-DeltaC
    D_new[j] = 2*DeltaD*j/niter-DeltaD
    theta_new=theta_true+c(0,0,0,0,0,0,C_new[i],D_new[j] ,0)
    llbis[i,j]=logLik_FRwoutNoise(theta_new,data)
  }
}

contour(theta_true[7]+C_new,theta_true[8]+D_new,llbis,nlevels=20,xlab="C",ylab="D")
custom_levels=quantile(llbis,probs=c(0.025,0.05,0.075,0.1,0.25,0.5,0.75),na.rm=T)
contour(theta_true[7]+C_new,theta_true[8]+D_new,llbis,levels=custom_levels,xlab="C",ylab="D")

### Perhaps I should compare this ridge to the one with the FR data. 
### Would allow to see if it is (likely) flatter

### Perhaps we should take a larger view
### 4-figure panel for comparison

theta_center = c(rmax_V,1/K,sqrt(0.05),rmax_P,Q,sqrt(0.05),C+2,D+2)

DeltaC = 4.5 # C is 4.5 +/- 4.5
DeltaD = 3.0 # D is 3 +/- 3

niter = 50
C_new=D_new=rep(0,niter)
llbis=matrix(0,nrow=niter,ncol=niter)
for (i in 1:niter){
  for (j in 1:niter){
    C_new[i] =2*DeltaC*i/niter-DeltaC
    D_new[j] = 2*DeltaD*j/niter-DeltaD
    theta_new=theta_center+c(0,0,0,0,0,0,C_new[i],D_new[j] ,0)
    llbis[i,j]=logLik_FRwoutNoise(theta_new,data)
  }
}

custom_levels=quantile(llbis,probs=c(0.025,0.05,0.075,0.1,0.25,0.5,0.75),na.rm=T)
contour(theta_center[7]+C_new,theta_center[8]+D_new,llbis,levels=custom_levels,xlab="C",ylab="D",main="Without FR data")
abline(v=2.5,col="red")
abline(h=1,col="red")

### Let's make that zoom much more precise - around the putative MLE
theta_true= c(rmax_V,1/K,sqrt(0.05),rmax_P,Q,sqrt(0.05),C+0.55,D+1,sqrt(0.05))## not really true
DeltaC = 0.25 # C is +/- 2
DeltaD = 0.5 #D is +/-0.9
niter = 50
C_new=D_new=rep(0,niter)
llbis=matrix(0,nrow=niter,ncol=niter)
for (i in 1:niter){
  for (j in 1:niter){
    C_new[i] =2*DeltaC*i/niter-DeltaC
    D_new[j] = 2*DeltaD*j/niter-DeltaD
    theta_new=theta_true+c(0,0,0,0,0,0,C_new[i],D_new[j] ,0)
    llbis[i,j]=logLik_FRwoutNoise(theta_new,data)
  }
}

custom_levels=quantile(llbis,probs=c(0.025,0.05,0.075,0.1,0.25,0.5,0.75),na.rm=T)
contour(theta_true[7]+C_new,theta_true[8]+D_new,llbis,levels=custom_levels,xlab="C",ylab="D",main="Zoom without FR data")
abline(v=2.5,col="red")
abline(h=1,col="red")

### Let's compare just below with what the LL with FR data looks like. 

theta_center = c(rmax_V,1/K,sqrt(0.05),rmax_P,Q,sqrt(0.05),C+2,D+2,sqrt(0.05))

DeltaC = 4.5 # C is 4.5 +/- 4.5
DeltaD = 3.0 # D is 3 +/- 3

niter = 50
C_new=D_new=rep(0,niter)
llbis=matrix(0,nrow=niter,ncol=niter)
for (i in 1:niter){
  for (j in 1:niter){
    C_new[i] =2*DeltaC*i/niter-DeltaC
    D_new[j] = 2*DeltaD*j/niter-DeltaD
    theta_new=theta_center+c(0,0,0,0,0,0,C_new[i],D_new[j] ,0)
    llbis[i,j]=logLik(theta_new,data)
  }
}
custom_levels=quantile(llbis,probs=c(0.025,0.05,0.075,0.1,0.25,0.5,0.75),na.rm=T)
contour(theta_center[7]+C_new,theta_center[8]+D_new,llbis,levels=custom_levels,xlab="C",ylab="D",main="With FR data")
abline(v=2.5,col="red")
abline(h=1,col="red")

### Let's make that zoom much more precise - around the true values
theta_true= c(rmax_V,1/K,sqrt(0.05),rmax_P,Q,sqrt(0.05),C,D,sqrt(0.05))## not really true

DeltaC = 0.25 # C is +/- 2
DeltaD = 0.5 #D is +/-0.9
niter = 50
C_new=D_new=rep(0,niter)
llbis=matrix(0,nrow=niter,ncol=niter)
for (i in 1:niter){
  for (j in 1:niter){
    C_new[i] =2*DeltaC*i/niter-DeltaC
    D_new[j] = 2*DeltaD*j/niter-DeltaD
    theta_new=theta_true+c(0,0,0,0,0,0,C_new[i],D_new[j] ,0)
    llbis[i,j]=logLik(theta_new,data)
  }
}

custom_levels=quantile(llbis,probs=c(0.025,0.05,0.075,0.1,0.25,0.5,0.75),na.rm=T)
contour(theta_true[7]+C_new,theta_true[8]+D_new,llbis,levels=custom_levels,xlab="C",ylab="D",main="Zoom with FR data")
abline(v=2.5,col="red")
abline(h=1,col="red")

### Much steeper but difficult to see in that representation. 
### Maybe I can show these as two surfaces of different colors, in 3D? 
### Perhaps I need the likelihood profiles after all. 

#### Zooming out - you never know
par(mfrow=c(1,2))
theta_center = c(rmax_V,1/K,sqrt(0.05),rmax_P,Q,sqrt(0.05),C+22,D+22)

DeltaC = 24.5 # 
DeltaD = 23.0 # 

niter = 100
C_new=D_new=rep(0,niter)
llbis=matrix(0,nrow=niter,ncol=niter)
for (i in 1:niter){
  for (j in 1:niter){
    C_new[i] =2*DeltaC*i/niter-DeltaC
    D_new[j] = 2*DeltaD*j/niter-DeltaD
    theta_new=theta_center+c(0,0,0,0,0,0,C_new[i],D_new[j] ,0)
    llbis[i,j]=logLik_FRwoutNoise(theta_new,data)
  }
}

custom_levels=quantile(llbis,probs=c(0.025,0.05,0.075,0.1,0.25,0.5,0.75),na.rm=T)
contour(theta_center[7]+C_new,theta_center[8]+D_new,llbis,levels=custom_levels,xlab="C",ylab="D",main="Without FR data")
abline(v=2.5,col="red")
abline(h=1,col="red")

theta_center = c(rmax_V,1/K,sqrt(0.05),rmax_P,Q,sqrt(0.05),C+22,D+22,sqrt(0.05))

niter = 100
C_new=D_new=rep(0,niter)
llbis=matrix(0,nrow=niter,ncol=niter)
for (i in 1:niter){
  for (j in 1:niter){
    C_new[i] =2*DeltaC*i/niter-DeltaC
    D_new[j] = 2*DeltaD*j/niter-DeltaD
    theta_new=theta_center+c(0,0,0,0,0,0,C_new[i],D_new[j] ,0)
    llbis[i,j]=logLik(theta_new,data)
  }
}

custom_levels=quantile(llbis,probs=c(0.025,0.05,0.075,0.1,0.25,0.5,0.75),na.rm=T)
contour(theta_center[7]+C_new,theta_center[8]+D_new,llbis,levels=custom_levels,xlab="C",ylab="D",main="With FR data")
abline(v=2.5,col="red")
abline(h=1,col="red")

### Not very easy to view differences from contours. 

#################### Old comments  ##################
### Does suggest a ridge there -- a little like in the Polanski paper, no? 
### We do manage to identify the model though! #

### + analytic computations in Mathematica or similar ?
### (should be able to work out the RSS derivations with Deriv )


