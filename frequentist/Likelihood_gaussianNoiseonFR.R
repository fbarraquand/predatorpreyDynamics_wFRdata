### FB 13/11/2017 - predator-prey model with noisy functional response data
### Writes down the likelihood of the model under the assumption of Gaussian noise on the FR
### Designed to be the simplest form noise tested first
### FB 09/01/2019 - better likelihood profiles for pairs of potentially correlated params

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

##################### Helper function #####################################################
### By OG 
logprot <- function(v){
  eps <- 2.2204e-016
  u <- log(eps) * (1+vector(length=length(v)))
  index <- (v>eps)
  u[index] <- log(v[index])
  u
}

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
d1=dnorm(y[t,1], mu1, theta[3],log=T) ## directly asking for the log avoids some problems
d2=dnorm(y[t,2], mu2, theta[6],log=T)
d3=dnorm(y[t,3], mu3, theta[9],log=T)
ll=ll+d1+d2+d3
}
return(-ll)
}

################ Now optimize the LL ###########################################
theta_start = c(runif(1,0.5,2),runif(1,0.001,5),runif(1,0.1,1),runif(1,0.1,2),runif(1,0.1,5),runif(1,0.1,2),runif(1,0.01,0.1),runif(1,10,100),runif(1,0.1,5))
# theta_start=rep(1,9) # you never know...
# Put high values for D not C -- previously the reverse was done for bayesian est. (error?)
theta_start = c(runif(1,0.5,2),runif(1,0.001,5),runif(1,0.05,1),runif(1,0.1,1),runif(1,0.1,5),runif(1,0.05,1),runif(1,0.01,0.1),runif(1,1,5),runif(1,0.05,1))
theta_true  = c(rmax_V,1/K,sqrt(0.05),rmax_P,Q,sqrt(0.05),C,D,sqrt(0.05))


p_opt<-optim(theta_start, logLik, y=data,method="BFGS",hessian=T)
### warnings initially: log(1 + theta[2] * N_t) : NaNs produced --> use of logprot to avoid this. 
p_opt$convergence
p_opt$fn
p_opt$par
theta_true
theta_start
# > p_opt$convergence
# [1] 0
# > p_opt$fn
# NULL
# > p_opt$par
# [1]  2.0293915  1.0343923  0.2225042  0.5062073 10.2952243  0.2314227  2.4326657  0.8129131  0.2191251
# > theta_true
# [1]  2.0000000  1.0000000  0.2236068  0.5000000 10.0000000  0.2236068  2.5000000  1.0000000  0.2236068
# > theta_start
# [1] 1.31112804 1.25949471 0.32802543 0.66675616 2.84856121 0.90796690 0.09700391 2.57365722 0.81605342
## Estimation quality can depend on initial values, especially sigmas (difficult if chosen too large)

## Still warnings()
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
# [1] 0
# > p_opt$fn
# NULL
# > p_opt$par
# [1]  2.0289516  1.0338520  0.2225055  0.5060431 10.2908387  0.2314251  2.4333060  0.8145374  0.2191241
# > theta_true
# [1]  2.0000000  1.0000000  0.2236068  0.5000000 10.0000000  0.2236068  2.5000000  1.0000000  0.2236068
# > theta_start
# [1] 1.31112804 1.25949471 0.32802543 0.66675616 2.84856121 0.90796690 0.09700391 2.57365722 0.81605342

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
# [1] 1.3369039 0.4461682 0.2225209 0.1509868 2.0933996 0.3150897 3.2261805 1.6031505 0.8683464
# > theta_true
# [1]  2.0000000  1.0000000  0.2236068  0.5000000 10.0000000  0.2236068  2.5000000  1.0000000  0.2236068
# > theta_start
# [1] 1.31112804 1.25949471 0.32802543 0.66675616 2.84856121 0.90796690 0.09700391 2.57365722 0.81605342
# no warnings this time but perhaps we should stick to BFGS?

library(optimx)
p_opt<-optimx(theta_start, logLik, y=data,hessian=T)
#slow - still same warnings
p_opt
theta_true
theta_start
# > p_opt
# p1        p2        p3        p4       p5        p6       p7        p8        p9     value fevals
# Nelder-Mead 1.336904 0.4461682 0.2225209 0.1509868  2.09340 0.3150897 3.226180 1.6031505 0.8683464  941.1084    502
# BFGS        2.029391 1.0343923 0.2225042 0.5062073 10.29522 0.2314227 2.432666 0.8129131 0.2191251 -227.4338    291
# gevals niter convcode  kkt1  kkt2  xtime
# Nelder-Mead     NA    NA        1 FALSE FALSE 11.472
# BFGS            55    NA        0  TRUE  TRUE 39.348
# > theta_true
# [1]  2.0000000  1.0000000  0.2236068  0.5000000 10.0000000  0.2236068  2.5000000  1.0000000  0.2236068
# > theta_start
# [1] 1.31112804 1.25949471 0.32802543 0.66675616 2.84856121 0.90796690 0.09700391 2.57365722 0.81605342

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
# [1]  2.0294189  1.0344203  0.2225043  0.5062008 10.2950547  0.2314247  2.4326617  0.8129034  0.2191242
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
# [1]  1.9953701  0.9919645  0.2222523  0.5168936 10.5795043  0.2316653  2.4375966  0.8283733  0.2185120
# > theta_true
# [1]  2.0000000  1.0000000  0.2236068  0.5000000 10.0000000  0.2236068  2.5000000  1.0000000  0.2236068
### In this straightforward case, Nelder-Mead works similarly to BGFS 
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
plot(C_FR,negll[,2])

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


##################  NB Should I get rid of the sigmas? ############
### Working directly with the sum of squares ######################

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
    mu1 = y[t-1,1] + theta[1] - log(1+theta[2]*N_t) - y[t-1,3]*P_t/N_t # corrected
    mu2 = y[t-1,2] + theta[3] - log((1+theta[4]*P_t/N_t))
    mu3 = (theta[5]*N_t)/(theta[6] + N_t)
    rss=rss+(y[t,1] - mu1)^2+(y[t,2]-mu2)^2+(y[t,3]-mu3)^2
  }
  return(rss)
}

### New theta true
theta_true  = c(rmax_V,1/K,rmax_P,Q,C,D)
theta_init = theta_true + rnorm(6,0,sd=0.01)
p_opt<-optim(theta_init, RSS, y=data,hessian=T)
p_opt$par
theta_true
### Here no problems with  In log(1 + theta[2] * N_t) : NaNs produced

### More noise on starting conditions
theta_init = theta_true + rnorm(6,0,sd=0.1)
p_opt<-optim(theta_init, RSS, y=data,hessian=T)
p_opt$par
theta_true

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



#################### Old comments  ##################
### Does suggest a ridge there -- a little like in the Polanski paper, no? 
### We do manage to identify the model though! 

########## Previously with an error I condluded ######################
### The likelihood seems to favor unusually high r and unusally low C
### We need to check this for other datasets -- longer datasets too...
#######################################################################

### + analytic computations in Mathematica or similar 
### (should be able to work out the RSS derivations with Deriv )


