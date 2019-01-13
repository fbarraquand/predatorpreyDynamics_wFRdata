## FB 12/09/2019 - compute the Jacobian matrix around the fixed point to know 
## For Leslie-type predator-prey model with Holling type II response

rm(list=ls())
graphics.off()

## This is the Jacobian on the log-scale because we consider Gaussian small noise on that scale

#### Parameter sets definition
# Parameter set 1
rmax_V<-2		
K<-1			 
rmax_P<-0.5
Q<-10
### FR parameters
C<-2.5
D<-1
theta1 = c(rmax_V,K,rmax_P,Q,C,D)

# Parameter set 2
rmax_V<-2		
K<-1			
rmax_P<-0.5
Q<-10
### FR parameters
C<-15
D<-0.25
theta2 = c(rmax_V,K,rmax_P,Q,C,D)

### Function that defines the prey equilibrium
prey_density_func <- function(x,theta){

r = theta[1]
K = theta[2]
s = theta[3]
Q = theta[4]
C = theta[5]
D = theta[6]

A = C*(exp(s)-1)/Q

f = exp(r-A*x/(D+x)) - (1+x/K)

return(f)

}

#test function
prey_density_func(1,theta1)

x=seq(0.001,30,by=0.1)
plot(x,prey_density_func(x,theta1))
abline(h=0,col="red")

x=seq(4,6,by=0.01)
plot(x,prey_density_func(x,theta1))
abline(h=0,col="red")

#library(uniroot)
#N_star = uniroot(prey_density_func,c(0.01,80),theta=theta1)
### Obviously can't be installed...
#https://www.rdocumentation.org/packages/rootSolve/versions/1.7/topics/uniroot.all
library(rootSolve)
N_star = uniroot.all(prey_density_func,c(3,7),theta=theta1)
N_star

predator_density <- function(theta,N){
return(N*(exp(theta[3])-1)/theta[4])
}

P_star=predator_density(theta1,N_star)

### Now we define the Jacobian

Jacobian <- function(theta,N,P){
r = theta[1]
K = theta[2]
s = theta[3]
Q = theta[4]
C = theta[5]
D = theta[6]

J=matrix(0,2,2)

J[1,1] = 1 + N*( C*P/((D+N)^2) - 1/(K+N) )

J[1,2] = - C*P/(D+N)

J[2,1] = Q*(P/N)/(1+Q*(P/N))#(exp(s)-1)/exp(s)

J[2,2] = 1 - Q*(P/N)/(1+Q*(P/N))#(exp(s)-1)/exp(s)

return(J)
}

J=Jacobian(theta1,N_star,P_star)
eigen(J)
## oscillatory decay to equilibrium

### Modulus of the eigenvalues
Mod(eigen(J)$values)
# stable
### Arg
Arg(eigen(J)$values) 


#### New parameter set
### Checking parameter 2

## Finding prey equilibrium
x=seq(0.02,3,by=0.01)
plot(x,prey_density_func(x,theta2))
abline(h=0,col="red")
N_star = uniroot.all(prey_density_func,c(0.02,4),theta=theta2)
N_star

## Finding predator equilibrium
P_star=predator_density(theta2,N_star)
P_star

## Finding Jacobian
J=Jacobian(theta2,N_star,P_star)
eigen(J)
## Imaginary part non-zero
Mod(eigen(J)$values)
abs(eigen(J)$values)
## Just above 1! Unstable fixed point. Settling on the invariant loop. 


