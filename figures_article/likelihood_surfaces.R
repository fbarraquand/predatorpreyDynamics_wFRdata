### Plot likelihood surfaces (contours) for the various pairs of parameters. 

rm(list=ls())
graphics.off()

# Loading one dataset of length T=100

# taking simuls from 
file_path1 = "../simulations/small_noise_on_FR/parameter_sets/perturbed_fixed_point/T=100/"
file_path2 = "../simulations/small_noise_on_FR/parameter_sets/noisy_cycles/T=100/"

### Datasets of the form
### data = cbind(log(N),log(P),FR) 

krep=1 # sample chosen
data1 = read.csv(paste(file_path1,"predatorPrey_withGaussianFR",krep,".csv",sep=""))
names(data1)=c("Time","n","p","KR")
data1[,1]<-NULL
head(data1)

# switch to noisy LC dataset
data2 = read.csv(paste(file_path2,"predatorPrey_withGaussianFR",krep,".csv",sep=""))
names(data2)=c("Time","n","p","KR")
data2[,1]<-NULL
head(data2)

## Sourcing the functions to compute the likelihood
source('../frequentist/LikelihoodFunctions.R')
source('figlabel.R')
cex_labels = 1.2

### Parameters of the model - fixed point parameter set
K<-1			# threshold dd 
rmax_V<-2			# Max AVERAGE growth rate (thus not a true max...)
rmax_P<-0.5
sigma2.proc<-0.05	
### FR and predator parameters
C<-2.5
D<-1
Q<-10

theta_true1  = c(rmax_V,1/K,sqrt(0.05),rmax_P,Q,sqrt(0.05),C,D,sqrt(0.05))
theta_true2  = c(rmax_V,1/K,sqrt(0.05),rmax_P,Q,sqrt(0.05),C,D)

pdf(file="likelihood_surfaces_FP.pdf",width=12,height=8)

par(mfrow=c(2,3))

### r, gamma
niter = 50
r_new=g_new=rep(0,niter)
llbis=matrix(0,nrow=niter,ncol=niter)
Deltar = 0.75 # rmax_V is +/- 1.5
Deltag = 0.9 #1/K is +/- 0.9
for (i in 1:niter){
  for (j in 1:niter){
    r_new[i] =2*Deltar*i/niter-Deltar
    g_new[j] = 2*Deltag*j/niter-Deltag
    theta_new=theta_true1+c(r_new[i],g_new[j],0,0,0,0,0,0,0)
    llbis[i,j]=logLik(theta_new,data1)
  }
}
custom_levels=quantile(llbis,probs=c(0.025,0.05,0.075,0.1,0.25,0.5,0.75),na.rm=T)
contour(theta_true1[1]+r_new,theta_true1[2]+g_new,llbis,levels=custom_levels,xlab="r",ylab=expression(gamma))
fig_label("(A)",pos="topleft",cex=cex_labels)

### s,Q
niter = 50
rP_new=q_new=rep(0,niter)
llbis=matrix(0,nrow=niter,ncol=niter)
Deltar = 0.2 # rmax_P is +/- 0.3
Deltaq = 5 # Q is +/- 2
for (i in 1:niter){
  for (j in 1:niter){
    rP_new[i] =2*Deltar*i/niter-Deltar
    q_new[j] = 2*Deltaq*j/niter-Deltaq
    theta_new=theta_true1+c(0,0,0,rP_new[i],q_new[j],0,0,0,0)
    llbis[i,j]=logLik(theta_new,data1)
  }
}
custom_levels=quantile(llbis,probs=c(0.025,0.05,0.075,0.1,0.25,0.5,0.75),na.rm=T)
contour(theta_true1[4]+rP_new,theta_true1[5]+q_new,llbis,levels=custom_levels,xlab="s",ylab="Q")
fig_label("(B)",pos="topleft",cex=cex_labels)

### C,D
DeltaC = 0.25 # C is +/- 1 -> produces not well defined enough maximum
DeltaD = 0.5 #D is +/-0.9
niter = 50
C_new=D_new=rep(0,niter)
llbis=matrix(0,nrow=niter,ncol=niter)
for (i in 1:niter){
  for (j in 1:niter){
    C_new[i] =2*DeltaC*i/niter-DeltaC
    D_new[j] = 2*DeltaD*j/niter-DeltaD
    theta_new=theta_true1+c(0,0,0,0,0,0,C_new[i],D_new[j] ,0)
    llbis[i,j]=logLik(theta_new,data1)
  }
}
#hist(llbis) ## what values are in there
min(llbis)
custom_levels=quantile(llbis,probs=c(0.025,0.05,0.075,0.1,0.25,0.5,0.75),na.rm=T)
#contour(theta_true[7]+C_new,theta_true[8]+D_new,llbis,nlevels=20,xlab="C",ylab="D")
contour(theta_true1[7]+C_new,theta_true1[8]+D_new,llbis,levels=custom_levels,xlab="C",ylab="D")
fig_label("(C)",pos="topleft",cex=cex_labels)

### --------------- Without the KR data ------------------------ 

### r, gamma
niter = 50
r_new=g_new=rep(0,niter)
llbis=matrix(0,nrow=niter,ncol=niter)
Deltar = 0.75 # rmax_V is +/- 1.5
Deltag = 0.9 #1/K is +/- 0.9
for (i in 1:niter){
  for (j in 1:niter){
    r_new[i] =2*Deltar*i/niter-Deltar
    g_new[j] = 2*Deltag*j/niter-Deltag
    theta_new=theta_true2+c(r_new[i],g_new[j],0,0,0,0,0,0)
    llbis[i,j]=logLik_FRwoutNoise(theta_new,data1)
  }
}
custom_levels=quantile(llbis,probs=c(0.025,0.05,0.075,0.1,0.25,0.5,0.75),na.rm=T)
contour(theta_true2[1]+r_new,theta_true2[2]+g_new,llbis,nlevels=50,xlab="r",ylab=expression(gamma))
fig_label("(D)",pos="topleft",cex=cex_labels)

### s,Q
niter = 50
rP_new=q_new=rep(0,niter)
llbis=matrix(0,nrow=niter,ncol=niter)
Deltar = 0.2 # rmax_P is +/- 0.3
Deltaq = 5 # Q is +/- 2
for (i in 1:niter){
  for (j in 1:niter){
    rP_new[i] =2*Deltar*i/niter-Deltar
    q_new[j] = 2*Deltaq*j/niter-Deltaq
    theta_new=theta_true2+c(0,0,0,rP_new[i],q_new[j],0,0,0,0)
    llbis[i,j]=logLik_FRwoutNoise(theta_new,data1)
  }
}
custom_levels=quantile(llbis,probs=c(0.025,0.05,0.075,0.1,0.25,0.5,0.75),na.rm=T)
contour(theta_true2[4]+rP_new,theta_true2[5]+q_new,llbis,levels=custom_levels,xlab="s",ylab="Q")
fig_label("(E)",pos="topleft",cex=cex_labels)

### C,D
DeltaC = 2 #0.25 # C is +/- 1 -> produces not well defined enough maximum
DeltaD = 0.95 #0.5 #D is +/-0.9
niter = 50
C_new=D_new=rep(0,niter)
llbis=matrix(0,nrow=niter,ncol=niter)
for (i in 1:niter){
  for (j in 1:niter){
    C_new[i] =2*DeltaC*i/niter-DeltaC
    D_new[j] = 2*DeltaD*j/niter-DeltaD
    theta_new=theta_true2+c(0,0,0,0,0,0,C_new[i],D_new[j] ,0)
    llbis[i,j]=logLik_FRwoutNoise(theta_new,data1)
  }
}
#hist(llbis) ## what values are in there
min(llbis)
custom_levels=quantile(llbis,probs=c(0.025,0.05,0.075,0.1,0.25,0.5,0.75),na.rm=T)
#contour(theta_true[7]+C_new,theta_true[8]+D_new,llbis,nlevels=20,xlab="C",ylab="D")
contour(theta_true2[7]+C_new,theta_true2[8]+D_new,llbis,levels=custom_levels,xlab="C",ylab="D")
fig_label("(F)",pos="topleft",cex=cex_labels)
dev.off()

#######################################################################################
#### ------------------- Now dataset 2 - noisy limit cycles ---------------------- ####
#######################################################################################

### Parameters of the model - fixed point parameter set
K<-1			# threshold dd 
rmax_V<-2			# Max AVERAGE growth rate (thus not a true max...)
rmax_P<-0.5
sigma2.proc<-0.05	
### FR and predator parameters
C<-15
D<-0.25
Q<-10

theta_true1  = c(rmax_V,1/K,sqrt(0.05),rmax_P,Q,sqrt(0.05),C,D,sqrt(0.05))
theta_true2  = c(rmax_V,1/K,sqrt(0.05),rmax_P,Q,sqrt(0.05),C,D)

pdf(file="likelihood_surfaces_LC.pdf",width=12,height=8)

par(mfrow=c(2,3))

### r, gamma
niter = 50
r_new=g_new=rep(0,niter)
llbis=matrix(0,nrow=niter,ncol=niter)
Deltar = 0.75 # rmax_V is +/- 1.5
Deltag = 0.9 #1/K is +/- 0.9
for (i in 1:niter){
  for (j in 1:niter){
    r_new[i] =2*Deltar*i/niter-Deltar
    g_new[j] = 2*Deltag*j/niter-Deltag
    theta_new=theta_true1+c(r_new[i],g_new[j],0,0,0,0,0,0,0)
    llbis[i,j]=logLik(theta_new,data2)
  }
}
custom_levels=quantile(llbis,probs=c(0.025,0.05,0.075,0.1,0.25,0.5,0.75),na.rm=T)
contour(theta_true1[1]+r_new,theta_true1[2]+g_new,llbis,levels=custom_levels,xlab="r",ylab=expression(gamma))
fig_label("(A)",pos="topleft",cex=cex_labels)

### s,Q
niter = 50
rP_new=q_new=rep(0,niter)
llbis=matrix(0,nrow=niter,ncol=niter)
Deltar = 0.2 # rmax_P is +/- 0.3
Deltaq = 5 # Q is +/- 2
for (i in 1:niter){
  for (j in 1:niter){
    rP_new[i] =2*Deltar*i/niter-Deltar
    q_new[j] = 2*Deltaq*j/niter-Deltaq
    theta_new=theta_true1+c(0,0,0,rP_new[i],q_new[j],0,0,0,0)
    llbis[i,j]=logLik(theta_new,data2)
  }
}
custom_levels=quantile(llbis,probs=c(0.025,0.05,0.075,0.1,0.25,0.5,0.75),na.rm=T)
contour(theta_true1[4]+rP_new,theta_true1[5]+q_new,llbis,levels=custom_levels,xlab="s",ylab="Q")
fig_label("(B)",pos="topleft",cex=cex_labels)

### C,D
DeltaC = 3 # C is +/- 1 
DeltaD = 0.25 #D is +/-0.9
niter = 50
C_new=D_new=rep(0,niter)
llbis=matrix(0,nrow=niter,ncol=niter)
for (i in 1:niter){
  for (j in 1:niter){
    C_new[i] =2*DeltaC*i/niter-DeltaC
    D_new[j] = 2*DeltaD*j/niter-DeltaD
    theta_new=theta_true1+c(0,0,0,0,0,0,C_new[i],D_new[j] ,0)
    llbis[i,j]=logLik(theta_new,data2)
  }
}
#hist(llbis) ## what values are in there
min(llbis)
custom_levels=quantile(llbis,probs=c(0.025,0.05,0.075,0.1,0.25,0.5,0.75),na.rm=T)
#contour(theta_true[7]+C_new,theta_true[8]+D_new,llbis,nlevels=20,xlab="C",ylab="D")
contour(theta_true1[7]+C_new,theta_true1[8]+D_new,llbis,levels=custom_levels,xlab="C",ylab="D")
fig_label("(C)",pos="topleft",cex=cex_labels)

### --------------- Without the KR data ------------------------ 

### r, gamma
niter = 50
r_new=g_new=rep(0,niter)
llbis=matrix(0,nrow=niter,ncol=niter)
Deltar = 0.75 # rmax_V is +/- 1.5
Deltag = 0.9 #1/K is +/- 0.9
for (i in 1:niter){
  for (j in 1:niter){
    r_new[i] =2*Deltar*i/niter-Deltar
    g_new[j] = 2*Deltag*j/niter-Deltag
    theta_new=theta_true2+c(r_new[i],g_new[j],0,0,0,0,0,0)
    llbis[i,j]=logLik_FRwoutNoise(theta_new,data2)
  }
}
custom_levels=quantile(llbis,probs=c(0.025,0.05,0.075,0.1,0.25,0.5,0.75),na.rm=T)
contour(theta_true2[1]+r_new,theta_true2[2]+g_new,llbis,nlevels=50,xlab="r",ylab=expression(gamma))
fig_label("(D)",pos="topleft",cex=cex_labels)

### s,Q
niter = 50
rP_new=q_new=rep(0,niter)
llbis=matrix(0,nrow=niter,ncol=niter)
Deltar = 0.2 # rmax_P is +/- 0.3
Deltaq = 5 # Q is +/- 2
for (i in 1:niter){
  for (j in 1:niter){
    rP_new[i] =2*Deltar*i/niter-Deltar
    q_new[j] = 2*Deltaq*j/niter-Deltaq
    theta_new=theta_true2+c(0,0,0,rP_new[i],q_new[j],0,0,0,0)
    llbis[i,j]=logLik_FRwoutNoise(theta_new,data2)
  }
}
custom_levels=quantile(llbis,probs=c(0.025,0.05,0.075,0.1,0.25,0.5,0.75),na.rm=T)
contour(theta_true2[4]+rP_new,theta_true2[5]+q_new,llbis,levels=custom_levels,xlab="s",ylab="Q")
fig_label("(E)",pos="topleft",cex=cex_labels)

### C,D
DeltaC = 3 #0.25 # C is +/- 1 -> produces not well defined enough maximum
DeltaD = 0.20 #0.5 #D is +/-0.9
niter = 50
C_new=D_new=rep(0,niter)
llbis=matrix(0,nrow=niter,ncol=niter)
for (i in 1:niter){
  for (j in 1:niter){
    C_new[i] =2*DeltaC*i/niter-DeltaC
    D_new[j] = 2*DeltaD*j/niter-DeltaD
    theta_new=theta_true2+c(0,0,0,0,0,0,C_new[i],D_new[j] ,0)
    llbis[i,j]=logLik_FRwoutNoise(theta_new,data2)
  }
}
#hist(llbis) ## what values are in there
min(llbis)
custom_levels=quantile(llbis,probs=c(0.025,0.05,0.075,0.1,0.25,0.5,0.75),na.rm=T)
#contour(theta_true[7]+C_new,theta_true[8]+D_new,llbis,nlevels=20,xlab="C",ylab="D")
contour(theta_true2[7]+C_new,theta_true2[8]+D_new,llbis,levels=custom_levels,xlab="C",ylab="D")
fig_label("(F)",pos="topleft",cex=cex_labels)
dev.off()