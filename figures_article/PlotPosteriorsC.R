### FB 14/03/2019 -- Analysis of data rarefaction experiment


rm(list=ls())
graphics.off()

file_path1 = "/media/frederic/DATA/Simuls_wOlivier/predatorprey_KRdata/perturbed_fixed_point/sigma=0.05/"
file_path2 = "/media/frederic/DATA/Simuls_wOlivier/predatorprey_KRdata/noisy_limit_cycles/sigma=0.05/"


timemax = c(100,50,25)
p_KR = c(1,0.25,0)
nrep = 100

par(mfrow=c(3,3))

for (lt in 1:length(timemax))
  {
  print(lt)
  
  for (lkr in 1:length(p_KR))
    {
    print(lkr)
    estim_mean<-NULL
    for (krep in 1:nrep){
      
      l=try(load(paste(file_path1,"T=",as.character(timemax[lt]),"/p_KR=",as.character(p_KR[lkr]),"/predpreyJAGS",krep,".RData",sep="")))
      if (class(l)=="try-error"){next}
      
      estim_mean <- rbind(estim_mean,out$BUGSoutput$mean$C)
      if (krep==1){
        y=density(out$BUGSoutput$sims.list$C,bw=0.05)
        plot(y,ylim=c(0,10),lwd=0.5,xlim=c(0,5),ylab="Pr(C|data)", xlab=paste("T = ",timemax[lt],", p_KR = ",p_KR[lkr],sep=""),main="",col="grey")
        abline(v=2.5,col="red",lwd=3)
      } else {
        lines(density(out$BUGSoutput$sims.list$C,bw=0.05),lwd=0.2,col="grey")
      }
    }
  }
}

## Now the limit cycle case. Will go to an Appendix. 


par(mfrow=c(3,3))

for (lt in 1:length(timemax))
{
  print(lt)
  
  for (lkr in 1:length(p_KR))
  {
    print(lkr)
    estim_mean<-NULL
    for (krep in 1:nrep){
      
      l=try(load(paste(file_path2,"T=",as.character(timemax[lt]),"/p_KR=",as.character(p_KR[lkr]),"/predpreyJAGS",krep,".RData",sep="")))
      if (class(l)=="try-error"){next}
      
      estim_mean <- rbind(estim_mean,out$BUGSoutput$mean$C)
      if (krep==1){
        y=density(out2$BUGSoutput$sims.list$C,bw=0.2)
        plot(y,ylim=c(0,4),lwd=0.5,xlim=c(0,20),ylab="Pr(C|data)", xlab=paste("T = ",timemax[lt],", p_KR = ",p_KR[lkr],sep=""),main="",col="grey")
        abline(v=15,col="red",lwd=3)
      } else {
        lines(density(out2$BUGSoutput$sims.list$C,bw=0.2),lwd=0.2,col="grey")
      }
    }
  }
}

## Now for D?


par(mfrow=c(3,3))

for (lt in 1:length(timemax))
{
  print(lt)
  
  for (lkr in 1:length(p_KR))
  {
    print(lkr)
    estim_mean<-NULL
    for (krep in 1:nrep){
      
      l=try(load(paste(file_path1,"T=",as.character(timemax[lt]),"/p_KR=",as.character(p_KR[lkr]),"/predpreyJAGS",krep,".RData",sep="")))
      if (class(l)=="try-error"){next}
      
      estim_mean <- rbind(estim_mean,out$BUGSoutput$mean$D)
      if (krep==1){
        y=density(out$BUGSoutput$sims.list$D,bw=0.05)
        plot(y,ylim=c(0,10),lwd=0.5,xlim=c(0,3),ylab="Pr(D|data)", xlab=paste("T = ",timemax[lt],", p_KR = ",p_KR[lkr],sep=""),main="",col="grey")
        abline(v=1,col="red",lwd=3)
      } else {
        lines(density(out$BUGSoutput$sims.list$D,bw=0.05),lwd=0.2,col="grey")
      }
    }
  }
}

## Now the limit cycle case. Will go to an Appendix. 


par(mfrow=c(3,3))

for (lt in 1:length(timemax))
{
  print(lt)
  
  for (lkr in 1:length(p_KR))
  {
    print(lkr)
    estim_mean<-NULL
    for (krep in 1:nrep){
      
      l=try(load(paste(file_path2,"T=",as.character(timemax[lt]),"/p_KR=",as.character(p_KR[lkr]),"/predpreyJAGS",krep,".RData",sep="")))
      if (class(l)=="try-error"){next}
      
      estim_mean <- rbind(estim_mean,out$BUGSoutput$mean$D)
      if (krep==1){
        y=density(out2$BUGSoutput$sims.list$D,bw=0.01)
        plot(y,ylim=c(0,10),lwd=0.5,xlim=c(0,1),ylab="Pr(D|data)", xlab=paste("T = ",timemax[lt],", p_KR = ",p_KR[lkr],sep=""),main="",col="grey")
        abline(v=0.25,col="red",lwd=3)
      } else {
        lines(density(out2$BUGSoutput$sims.list$D,bw=0.01),lwd=0.2,col="grey")
      }
    }
  }
}

## Do we need to consider only the mean? 

