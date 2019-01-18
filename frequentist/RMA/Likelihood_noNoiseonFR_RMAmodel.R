############### Discrete-time RMA predator-prey model ###########################################
### FB 13/11/2017 - predator-prey model with noisy functional response data
### Writes down the likelihood of the model under the assumption of Gaussian noise on the FR
### 24/11/2017 Rosenzweig-MacArthur version -- should improve fitting but doesnt seem to 
### FB 12/03/2018 - back to classic predator-prey model fit with no noise on functional response data
### FB 16/01/2019 -- edited so that functional response timing match math model (no delays)

rm(list=ls())
graphics.off()

### Model used (use stored simulated data later on, just tryouts for now)
n.years<-1000  	# Number of years - 25 first, perhaps use 50 or 100 / worked very well with almost no process error on the real scale
N1<-1			# Initial pop size
P1<-0.1
K<-100			# threshold dd 1/K
beta<-1			# density-dependence exponent
rmax_V<-2		# Max AVERAGE growth rate (thus not a true max...)
rmax_P<-(-0.25)	        # Mortality of the predator #-0.25
sigma2.proc<-0.05	# worked well with 0.005
# Process sigma on the log-scale, use the Peretti et al. value. 0.005


### FR and predator parameters
C<-2.5
D<-1
epsilon<-0.1

### Simulation of data
#set.seed(42) 
#set.seed(41)
set.seed(40)

y<-N<-P<-FR<-numeric(n.years)
N[1]<-N1
P[1]<-P1
FR[1]<-C*N[1]/(D+N[1])

rV<-rnorm(n.years-1,rmax_V,sqrt(sigma2.proc))
rP<-rnorm(n.years-1,rmax_P,sqrt(0.01))
### FRnoise<-rnorm(n.years,0,sqrt(sigma2.proc)) #not present here

for (t in 1:(n.years-1)){
 
  N[t+1]<-N[t]*(exp(rV[t]-(N[t]/K)^beta -FR[t]*P[t]/N[t]) + 0.01) #+0.01 to maintain the prey pop positive
  P[t+1]<-P[t]*exp(rP[t] + epsilon*FR[t])
  FR[t+1]<-max(C*N[t+1]/(D+N[t+1]),0.0001) #+ FRnoise[t+1],0.0)
}
## Plotting time series of abundances and FR
par(mfrow=c(2,2))
plot(1:n.years,N,type="b")
plot(1:n.years,P,type="b")
#curve(dbeta(x,a,b),from=0, to=1)
plot(N,FR)
plot(log(N),log(P))
### Produce dataset to fit
#data = cbind(log(N),log(P),FR) 
data = cbind(log(N),log(P)) 

### May be interesting later on to compare the type of dynamical behaviour to identifiability

################# Define the likelihood ########################################
logLik=function(theta,y){
#### y is the data with n rows and 3 columns // log-abundance data + FR data

#### Parameters
# theta1 = r 
# theta2 = gamma = 1/K
# theta3 = sigma1
# theta4 = mu  #mortality
# theta5 = epsilon
# theta6 = sigma2
# theta7 = C
# theta8 = D
# theta9 = sigma3

n=nrow(y)
ll = 0.0 ### Or p_1(a_1|x_1) p(x_1)
for (t in 2:n){
N_t = exp(y[t-1,1])
P_t = exp(y[t-1,2]) 

########## Found error here -- N_t/P_t instead of P_t/N_t in multiplier of FR ########
### mu1 = y[t-1,1] + theta[1] - theta[2]*N_t - y[t-1,3]*N_t/P_t ######################
### was written instead of 
### mu1 = y[t-1,1] + theta[1] - theta[2]*N_t - y[t-1,3]*P_t/N_t ######################
######################################################################################
#mu1 = y[t-1,1] + theta[1] - theta[2]*N_t - ((theta[7]*N_t)/(theta[8] + N_t))*P_t/N_t
mu1 = y[t-1,1] + theta[1] - theta[2]*N_t - (theta[7]*P_t)/(theta[8] + N_t)
### The new DD should remove previous problems with log(1+theta[2]*N_t) 
mu2 = y[t-1,2] + theta[4] + theta[5]*(theta[7]*N_t)/(theta[8] + N_t)
# mu3 = not defined here 
#ll= ll + log(dnorm(y[t,1], mu1, theta[3])) +  log(dnorm(y[t,2], mu2, theta[6])) + log(dnorm(y[t,3], mean = mu3, sd = theta[9]))
# we have log(0) problem
d1=dnorm(y[t,1], mu1, theta[3],log=T) ## directly asking for the log avoids problems
d2=dnorm(y[t,2], mu2, theta[6],log=T)
#d3=dnorm(y[t,3], mu3, theta[9],log=T)
ll=ll+d1+d2#+d3
}
return(-ll)
}

################ Now optimize the LL ###########################################
theta_start = c(runif(1,0.5,2),runif(1,0.001,5),runif(1,0.1,1),runif(1,0.1,0.5),runif(1,0.1,0.5),runif(1,0.1,2),runif(1,0.01,0.1),runif(1,10,100),runif(1,0.1,5))
# Put high values for D not C -- previously the reverse was done for bayesian est. (error?)
# theta_start=rep(1,9) # you never know...

p_opt<-optim(theta_start, logLik, y=data,method="BFGS",hessian=T)
#warnings(): 44: In dnorm(y[t, 2], mu2, theta[6], log = T) : NaNs produced
p_opt<-optim(theta_start, logLik, y=data,hessian=T)
# does not work -- unclear why? Use L-BFGS? 
library(optimx)
p_opt<-optimx(theta_start, logLik, y=data,hessian=T)
## Always same warnings -> In log(1 + theta[2] * N_t) : NaNs produced

### Let's try something extreme 
theta_true  = c(rmax_V,1/K,sqrt(0.05),rmax_P,epsilon,sqrt(0.05),C,D,sqrt(0.05))
p_opt<-optim(theta_true, logLik, y=data,method="BFGS",hessian=T)
p_opt$par
theta_true

############ Create likelihood profiles ########################################

# basic check 
logLik(theta_true,data)
llbis=matrix(0,nrow=10,ncol=10)
for (i in 1:10){
  for (j in 1:10){
    theta_new=theta_true+c(0,0,0,0,0,0,0.05*i-0.5,0.5*j-2,0)
    llbis[i,j]=logLik(theta_new,data)
  }
}
contour(llbis)

# do that with (r,gamma)
logLik(theta_true,data)
llbis=matrix(0,nrow=10,ncol=10)
r_new=g_new=rep(0,10)
for (i in 1:10){
  for (j in 1:10){
    r_new[i] = 1*i-0.1
    g_new[j] = 0.2*j-0.9
    theta_new=theta_true+c(r_new[i],g_new[j],0,0,0,0,0,0)
    llbis[i,j]=logLik(theta_new,data)
  }
}
contour(theta_true[1]+r_new,theta_true[2]+g_new,llbis)

# Let's assume we know the max growth rates (most realistic assumption)

gamma = c(0.1,0.05,0.1,0.5,1,2)
negll=matrix(0,nrow=100,ncol=101)

par(mfrow=c(2,3))
for (kg in 1:length(gamma)){
  D_FR = seq(0.5,200.5,by=2)
    for (kd in 1:length(D_FR)){
      C_FR = seq(0.05,5,0.05) 
      for (kc in 1:length(C_FR)){
       theta = c(rmax_V,gamma[kg],sqrt(0.05),rmax_P,epsilon,sqrt(0.05),C_FR[kc],D_FR[kg],sqrt(0.05))
       negll[kc,kd]=-logLik(theta,data)
      }
    }
  contour(C_FR,D_FR,negll,xlab="C",ylab="D", main=c("Gamma = ",toString(gamma[kg])))
}
plot(C_FR,negll[,2])

### Zooming in on C (between 2 and 3)
par(mfrow=c(2,3))
for (kg in 1:length(gamma)){
  D_FR = seq(0.5,200.5,by=2)
  for (kd in 1:length(D_FR)){
    C_FR = seq(2.01,3,0.01) 
    for (kc in 1:length(C_FR)){
      theta = c(rmax_V,gamma[kg],sqrt(0.05),rmax_P,epsilon,sqrt(0.05),C_FR[kc],D_FR[kg],sqrt(0.05))
      negll[kc,kd]=-logLik(theta,data)
    }
  }
  contour(C_FR,D_FR,negll,xlab="C",ylab="D", main=c("Gamma = ",toString(gamma[kg])))
}

plot(C_FR,negll[,2])

plot(C_FR,negll[,100])

########### Other more realistic profiles #####################################
# todo
###############################################################################   

### NB Should I normalize everything and get rid of the sigmas? 

############### Working directly with the sum of squares ######################
RSS=function(theta,y){
  #### y is the data with n rows and 3 columns // log-abundance data + FR data
  
  #### Parameters
  # theta1 = r 
  # theta2 = gamma
  # theta3 = mu
  # theta4 = epsilon 
  # theta5 = C
  # theta6 = D
  # not the same theta
  
  n=nrow(y)
  rss = 0.0 ### Or p_1(a_1|x_1) p(x_1)
  for (t in 2:n){
    N_t = exp(y[t-1,1])
    P_t = exp(y[t-1,2]) 
    ### mu1 = y[t-1,1] + theta[1] - theta[2]*N_t - y[t-1,3]*N_t/P_t ## Error corrected. 
    mu1 = y[t-1,1] + theta[1] - theta[2]*N_t - (theta[5]*P_t)/(theta[6] + N_t)
    mu2 = y[t-1,2] + theta[3] + theta[4]*(theta[5]*N_t)/(theta[6] + N_t)
    #mu3 = (theta[5]*N_t)/(theta[6] + N_t)
    rss=rss+(y[t,1] - mu1)^2+(y[t,2]-mu2)^2#+(y[t,3]-mu3)^2
  }
  return(rss)
}

### New theta true
theta_true  = c(rmax_V,1/K,rmax_P,epsilon,D,C)
theta_init = theta_true + rnorm(6,0,sd=0.01)
p_opt<-optim(theta_init, RSS, y=data,hessian=T)
p_opt$par
theta_true
### Fairly good... 

### Let's try some more challenging starting values
theta_init = theta_true + rnorm(6,0,sd=0.1)
p_opt<-optim(theta_init, RSS, y=data,hessian=T)
p_opt$par
theta_true
### Fairly good for a number of repeats! 


# basic check 
RSS(theta_true,data)
rssbis=matrix(0,nrow=10,ncol=10)
for (i in 1:10){
  for (j in 1:10){
    theta_new=theta_true+c(0,0,0,0,0.05*i-0.5,0.5*j-2)
    rssbis[i,j]=RSS(theta_new,data)
  }
}
contour(rssbis)

# do that with (r,gamma)
RSS(theta_true,data)
rssbis=matrix(0,nrow=10,ncol=10)
r_new=g_new=rep(0,10)
for (i in 1:10){
  for (j in 1:10){
    r_new[i] = 1*i-0.1
    g_new[j] = 0.2*j-0.9
    theta_new=theta_true+c(r_new[i],g_new[j],0,0,0,0)
    rssbis[i,j]=RSS(theta_new,data)
  }
}
contour(theta_true[1]+r_new,theta_true[2]+g_new,rssbis)

### 

