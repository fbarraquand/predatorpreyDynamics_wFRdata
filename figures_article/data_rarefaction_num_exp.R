### FB 07/03/2019 -- Loop to estimate predator-prey models with different time series lengths and amounts of KR data

rm(list=ls())
graphics.off()

library("R2jags")

file_path1 = "../simulations/small_noise_on_FR/parameter_sets/perturbed_fixed_point/T=100/"
file_path2 = "../simulations/small_noise_on_FR/parameter_sets/noisy_cycles/T=100/"
file_save = "/media/frederic/DATA/Simuls_wOlivier/predatorprey_KRdata/"

timemax = c(100,50,25)
p_KR = c(1,0.25,0)
nrep = 100

source('BayesianModels.R')
## code with the functions that are well encoded

    
    for (krep in 1:nrep){
      
      print(krep)
      # Pick one of the 100 repeats and load data
      # taking simuls from 
      
      ### Datasets of the form
      ### data = cbind(log(N),log(P),FR) 

      data1 = read.csv(paste(file_path1,"predatorPrey_withGaussianFR",krep,".csv",sep=""))
      names(data1)=c("Time","n","p","KR")
      head(data1)
      # switch to noisy LC dataset
      data2 = read.csv(paste(file_path2,"predatorPrey_withGaussianFR",krep,".csv",sep=""))
      names(data2)=c("Time","n","p","KR")
      head(data2)
      
      
      for (lt in 1:length(timemax)){
        print(lt)
        
        data1 = subset(data1,data1$Time<(timemax[lt]+1))
        data2 = subset(data2,data2$Time<(timemax[lt]+1))
        
        for (lkr in 1:length(p_KR)){
          print(lkr)
          
          if (p_KR[lkr]==0){ # Estimate the model without kill rate data, for both param sets
            ### [missing jags.data]
            n.years = timemax[lt]
            jags.data = list(T=n.years,logN=data1$n,logP=data1$p)
            while(TRUE){
              out <- try(jags(jags.data, inits, parameters, "predprey_without_sepFR.txt", n.chains = nc, n.iter = ni, n.burnin = nb, working.directory = getwd()), silent=TRUE)
              if(!is(out, 'try-error')) break
            }
            
          } else { # Estimate the model with kill rate data, for both param sets
            n.years = timemax[lt]
            nKR = round(p_KR[lkr]*timemax[lt]) # Number of points for the functional response
            FR_times = sample(data1$Time,nKR)
            masked = data1$Time[!data1$Time %in% FR_times]
            data1$KR[masked]<-NA  # [check this works well for the full vector]
            jags.data = list(T=n.years,logN=data1$n,logP=data1$p,FR=data1$KR)
            
            while(TRUE){
              out <- try(jags(jags.data, inits, parameters, "predprey.txt", n.chains = nc, n.iter = ni, n.burnin = nb, working.directory = getwd())) #, silent=TRUE
              if(!is(out, 'try-error')) break
            }
            print(out)
            save(out,file=paste(file_save,"perturbed_fixed_point/sigma=0.05/T=",as.character(timemax[lt]),"/p_KR=",as.character(p_KR[lkr]),"/predpreyJAGS",krep,".RData",sep=""))
            # Store in the appropriate repo 
            # Within sigma = 0.05, T=100,T=50,T=25; within these p_KR=1;p_KR=0.25,p_KR=0. 
            } #end of condition on p_KR
            
            ### Same instructions for data2 - the noisy limit cycle parameter set [OK, code could be more concise...]
            
            if (p_KR[lkr]==0){ # Estimate the model without kill rate data, 
              ### [missing jags.data]
              n.years = timemax[lt]
              jags.data = list(T=n.years,logN=data2$n,logP=data2$p)
              while(TRUE){
                out2 <- try(jags(jags.data, inits, parameters, "predprey_without_sepFR.txt", n.chains = nc, n.iter = ni, n.burnin = nb, working.directory = getwd()), silent=TRUE)
                if(!is(out2, 'try-error')) break
              }
              
            } else { # Estimate the model with kill rate data, for both param sets
             
              data2$KR[masked]<-NA  # [check this works well for the full vector]
              jags.data = list(T=n.years,logN=data2$n,logP=data2$p,FR=data2$KR)
              
              while(TRUE){
                out2 <- try(jags(jags.data, inits, parameters, "predprey.txt", n.chains = nc, n.iter = ni, n.burnin = nb, working.directory = getwd())) #, silent=TRUE
                if(!is(out2, 'try-error')) break
              }
              print(out2)
              save(out2,file=paste(file_save,"noisy_limit_cycles/sigma=0.05/T=",as.character(timemax[lt]),"/p_KR=",as.character(p_KR[lkr]),"/predpreyJAGS",krep,".RData",sep=""))
              # Store in the appropriate repo 
              # Within sigma = 0.05, T=100,T=50,T=25; within these p_KR=1;p_KR=0.25,p_KR=0. 
            
              } #end of condition on p_KR
            
    } #end of loop on p_KR
    
  } #end of loop on timemax
  
} #end of loop on krep


  
  