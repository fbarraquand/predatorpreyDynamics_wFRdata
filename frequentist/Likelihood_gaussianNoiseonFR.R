### FB 13/11/2017 - predator-prey model with noisy functional response data
### Writes down the likelihood of the model under the assumption of Gaussian noise on the FR

rm(list=ls())
graphics.off()

### Model used (use stored simulated data later on, just tryouts for now)
n.years<-100  	# Number of years - 25 first, perhaps use 50 or 100 / worked very well with almost no process error on the real scale
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
mu1 = y[t-1,1] + theta[1] - log(1+theta[2]*N_t) - y[t-1,3]*N_t/P_t
mu2 = y[t-1,2] + theta[4] - log((1+theta[5]*P_t/N_t))
mu3 = (theta[7]*N_t)/(theta[8] + N_t)
ll= ll + log(dnorm(y[t,1], mu1, theta[3])) +  log(dnorm(y[t,2], mu2, theta[6])) +log(dnorm(y[t,3], mean = mu1, sd = theta[9]))
}
return(ll)
}

################ Now optimize the LL ###########################################
theta_start = c(runif(1,0.5,2),runif(1,0.001,0.01),runif(1,0.1,2),runif(1,0.1,2),runif(1,0.1,5),runif(1,0.1,2), runif(1,10,100),runif(1,0.01,0.1),runif(1,1,10))
#theta=theta_start
#y=data

p_opt<-optim(theta_start, logLik, y=data,method="BFGS",hessian=T)

