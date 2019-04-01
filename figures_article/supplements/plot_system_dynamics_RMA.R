### FB 06/03/2019 - Plot illustrative time series 
### FB 27/03/2019 - RMA model

rm(list=ls())
graphics.off()

### Simulation

# ------------------------ FP dataset ---------------------------------- #

n.years<-200  	# Number of years - 25 first, perhaps use 50 or 100 / worked very well with almost no process error on the real scale
#N1<-1			# Initial pop size
#P1<-0.1 ### Too long transients
N1<-4
P1<-0.7
K<-1			# threshold dd 
beta<-1			# density-dependence exponent
rmax_V<-2			# Max AVERAGE growth rate (thus not a true max...)
rmax_P<-(-0.20) #
sigma2.proc<-0.05		# worked well with 0.005
# Process sigma on the log-scale, use the Peretti et al. value. 0.005

### FR and predator parameters
C<-2.5
D<-1.0 #0.6 = quasi-cycles
epsilon<-0.1

### Simulation of data
set.seed(42) 
#set.seed(41)
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
  P[t+1]<-P[t]*exp(rP[t]+epsilon*FR[t])
  FR[t+1]<-(C*N[t+1]/(D+N[t+1])) + FRnoise[t+1]
}

after_transients=101:200
data1=data.frame(Time=after_transients,logN=log(N[after_transients]),logP=log(P[after_transients]),FR=FR[after_transients])

# ------------------------ QC dataset ---------------------------------- #

### Parameters for simulation of Hassell model
n.years<-200  	# Number of years - 25 first, perhaps use 50 or 100 / worked very well with almost no process error on the real scale
#N1<-1			# Initial pop size
#P1<-0.1 ### Too long transients
N1<-2.5
P1<-0.9

K<-1			# threshold dd 
beta<-1			# density-dependence exponent
rmax_V<-2			# Max AVERAGE growth rate (thus not a true max...)
rmax_P<-(-0.20) #
sigma2.proc<-0.05		# worked well with 0.005
# Process sigma on the log-scale, use the Peretti et al. value. 0.005


### FR and predator parameters
C<-2.5
D<-0.5 #= quasi-cycles
epsilon<-0.1

### correction to avoid extinction and weird dynamics
eps <-0.1

### Simulation of data
set.seed(42) 
#set.seed(41)
#set.seed(40)

y<-N<-P<-FR<-numeric(n.years)
N[1]<-N1
P[1]<-P1
FR[1]<-C*N[1]/(D+N[1])

rV<-rnorm(n.years-1,rmax_V,sqrt(sigma2.proc))
rP<-rnorm(n.years-1,rmax_P,sqrt(sigma2.proc))
FRnoise<-rnorm(n.years-1,0,sqrt(sigma2.proc))

for (t in 1:(n.years-1)){
  N[t+1]<-N[t]*(exp(rV[t])/(1+(N[t]/K)^beta))*exp(-FR[t]*P[t]/(N[t]+eps))
  P[t+1]<-P[t]*exp(rP[t]+epsilon*FR[t])
  FR[t+1]<-(C*N[t+1]/(D+N[t+1])) + FRnoise[t+1]
}

after_transients=101:200
data2=data.frame(Time=after_transients,logN=log(N[after_transients]),logP=log(P[after_transients]),FR=FR[after_transients])

# ------------------------ LC dataset ---------------------------------- #

### Parameters for simulation of Hassell model
n.years<-200  	# Number of years - 25 first, perhaps use 50 or 100 / worked very well with almost no process error on the real scale
N1<-1			# Initial pop size
P1<-0.1 ### Too long transients
K<-1			# threshold dd 
beta<-1			# density-dependence exponent
rmax_V<-1.8			# Max AVERAGE growth rate (thus not a true max...)
rmax_P<-(-0.7) # 0.5 works with rmaxV=2
sigma2.proc<-0.05 #0.05		# worked well with 0.005
# Process sigma on the log-scale, use the Peretti et al. value. 0.005

### FR and predator parameters
C<-10
D<-0.6 #= quasi-cycles
epsilon<-0.1 # conversion efficiency
eps <-0.1 ### correction to avoid extinction and weird dynamics

### Simulation of data
set.seed(42)  #
#set.seed(41)
#set.seed(40)

y<-N<-P<-FR<-numeric(n.years)
N[1]<-N1
P[1]<-P1
FR[1]<-C*N[1]/(D+N[1])

rV<-rnorm(n.years-1,rmax_V,sqrt(sigma2.proc))
rP<-rnorm(n.years-1,rmax_P,sqrt(sigma2.proc))
FRnoise<-rnorm(n.years,0,sqrt(sigma2.proc))

for (t in 1:(n.years-1)){
  N[t+1]<-N[t]*(exp(rV[t])/(1+(N[t]/K)^beta))*exp(-FR[t]*P[t]/(N[t]+eps))
  P[t+1]<-P[t]*exp(rP[t]+epsilon*FR[t])
  FR[t+1]<-(C*N[t+1]/(D+N[t+1])) + FRnoise[t+1]
}

after_transients=101:200
data3=data.frame(Time=after_transients,logN=log(N[after_transients]),logP=log(P[after_transients]),FR=FR[after_transients])

# ----------------------------------------------------------------------- #

# Checking datasets

names(data1)=c("Time","n","p","KR")
head(data1)

names(data2)=c("Time","n","p","KR")
head(data2)

names(data3)=c("Time","n","p","KR")
head(data3)

source('../figlabel.R')

pdf(file = "system-dynamics_RMA.pdf",height=20,width=18)
par(mfrow=c(3,3),pch=19,cex=1.5)
# Better than homemade plotyy 
# https://www.r-bloggers.com/r-single-plot-with-two-different-y-axes/
cex_labels=1.2
with(data1,plot(Time,exp(n),col="blue",type="b",ylab="N"))
fig_label("(A)",pos="topleft",cex=cex_labels)
par(new = TRUE)
with(data1, plot(Time, exp(p), axes=F, col="red",xlab=NA, ylab=NA,type="b"))
axis(side = 4)
#mtext(side = 4, line = 3, 'P')
plot(exp(data1$n),exp(data1$p),type="b",xlab="N",ylab="P")
fig_label("(B)",cex=cex_labels)
plot(exp(data1$n),data1$KR,xlab="N",ylab="Kill rate",ylim=c(0,4))
holling2 <- function(x,C=2.5,D=1){C*x/(D+x)}
curve(holling2,col="blue",xlab=NA,ylab=NA,ylim=c(0,4),add=TRUE,lwd=2)
fig_label("(C)",cex=cex_labels)

with(data2,plot(Time,exp(n),col="blue",type="b",ylab="N"))
fig_label("(D)",cex=cex_labels)
par(new = TRUE)
with(data2, plot(Time, exp(p), axes=F, col="red",xlab=NA, ylab=NA,type="b"))
axis(side = 4)
#mtext(side = 4, line = 3, 'P')
plot(exp(data2$n),exp(data2$p),type="b",xlab="N",ylab="P")
fig_label("(E)",cex=cex_labels)
plot(exp(data2$n),data2$KR,xlab="N",ylab="Kill rate",ylim=c(0,4))
holling2_2 <- function(x,C=2.5,D=0.5){C*x/(D+x)}
curve(holling2_2,col="blue",xlab=NA,ylab=NA,ylim=c(0,4),add=TRUE,lwd=2)
fig_label("(F)",cex=cex_labels)

with(data3,plot(Time,exp(n),col="blue",type="b",ylab="N"))
fig_label("(G)",cex=cex_labels)
par(new = TRUE)
with(data3, plot(Time, exp(p), axes=F, col="red",xlab=NA, ylab=NA,type="b"))
axis(side = 4)
#mtext(side = 4, line = 3, 'P')
plot(exp(data3$n),exp(data3$p),type="b",xlab="N",ylab="P")
fig_label("(H)",cex=cex_labels)
plot(exp(data3$n),data3$KR,xlab="N",ylab="Kill rate",ylim=c(0,12))
holling2_3 <- function(x,C=10,D=0.6){C*x/(D+x)}
curve(holling2_3,col="blue",xlab=NA,ylab=NA,ylim=c(0,12),add=TRUE,lwd=2)
fig_label("(I)",cex=cex_labels)

dev.off()
