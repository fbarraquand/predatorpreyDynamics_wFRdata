### FB 06/03/2019 - Plot illustrative time series 

rm(list=ls())
graphics.off()

# taking simuls from 
file_path1 = "../simulations/small_noise_on_FR/parameter_sets/perturbed_fixed_point/T=100/"
file_path2 = "../simulations/small_noise_on_FR/parameter_sets/noisy_cycles/T=100/"

### Datasets of the form
### data = cbind(log(N),log(P),FR) 

krep=1 # sample chosen
data1 = read.csv(paste(file_path1,"predatorPrey_withGaussianFR",krep,".csv",sep=""))
names(data1)=c("Time","n","p","KR")
head(data1)

# switch to noisy LC dataset
data2 = read.csv(paste(file_path2,"predatorPrey_withGaussianFR",krep,".csv",sep=""))
names(data2)=c("Time","n","p","KR")
head(data2)

source('figlabel.R')

pdf(file = "system-dynamics.pdf",height=12,width=18)
par(mfrow=c(2,3),pch=19,cex=1.5)
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
plot(exp(data1$n),data1$KR,xlab="N",ylab="Kill rate",ylim=c(0,3))
holling2 <- function(x,C=2.5,D=1){C*x/(D+x)}
curve(holling2,col="blue",xlab=NA,ylab=NA,ylim=c(0,3),add=TRUE,lwd=2)
fig_label("(C)",cex=cex_labels)

with(data2,plot(Time,exp(n),col="blue",type="b",ylab="N"))
fig_label("(D)",cex=cex_labels)
par(new = TRUE)
with(data2, plot(Time, exp(p), axes=F, col="red",xlab=NA, ylab=NA,type="b"))
axis(side = 4)
#mtext(side = 4, line = 3, 'P')
plot(exp(data2$n),exp(data2$p),type="b",xlab="N",ylab="P")
fig_label("(E)",cex=cex_labels)
plot(exp(data2$n),data2$KR,xlab="N",ylab="Kill rate",ylim=c(0,20))
holling2_2 <- function(x,C=15,D=0.25){C*x/(D+x)}
curve(holling2_2,col="blue",xlab=NA,ylab=NA,ylim=c(0,20),add=TRUE,lwd=2)
fig_label("(F)",cex=cex_labels)
dev.off()

### NB checking the deterministic dynamics -- it is not very damped... 
file_path1 = "../simulations/small_noise_on_FR/parameter_sets/perturbed_fixed_point/"
data1 = read.csv(paste(file_path1,"predatorPrey_withGaussianFRdata_longSim_noNoise.csv",sep=""))
names(data1)=c("Time","n","p","KR")
head(data1)
data1=subset(data1,data1$Time<101)
nrow(data1)

plot(data1$Time,exp(data1$n))

par(mfrow=c(2,3),pch=19,cex=1.2)
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
plot(exp(data1$n),data1$KR,xlab="N",ylab="Kill rate",ylim=c(0,3))
holling2 <- function(x,C=2.5,D=1){C*x/(D+x)}
curve(holling2,col="blue",xlab=NA,ylab=NA,ylim=c(0,3),add=TRUE,lwd=2)
fig_label("(C)",cex=cex_labels)
