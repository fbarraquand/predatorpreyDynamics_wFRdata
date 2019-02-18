
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



################# Define the likelihood ########################################
logLik_reparam=function(theta,y){
  #### y is the data with n rows and 3 columns // log-abundance data + FR data
  
  #### Parameters
  # theta1 = r 
  # theta2 = K
  # theta3 = sigma1
  # theta4 = s
  # theta5 = q 
  # theta6 = sigma2
  # theta7 = a
  # theta8 = h
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
    mu1 = y[t-1,1] + theta[1] - logprot(1+(exp(theta[1])-1)*N_t/theta[2]) - y[t-1,3]*P_t/N_t #this is correct timing
    mu2 = y[t-1,2] + theta[4] - logprot((1+theta[5]*(exp(theta[4])-1)*P_t/N_t))
    mu3 = (theta[7]*N_tplus1)/(1 + theta[7]*theta[8]*N_tplus1) #easier to update all variables simultaneously, minimizes errors
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

logLik_FRwoutNoise_reparam=function(theta,y){
  #### y is the data with n rows and 3 columns // log-abundance data + FR data
  
  #### Parameters
  # theta1 = r 
  # theta2 = K
  # theta3 = sigma1
  # theta4 = s
  # theta5 = q 
  # theta6 = sigma2
  # theta7 = a
  # theta8 = h
  # theta9 = sigma3
  
  n=nrow(y)
  ll = 0.0 ### Or p_1(a_1|x_1) p(x_1)
  for (t in 2:n){
    N_t = exp(y[t-1,1])
    P_t = exp(y[t-1,2]) 
    
    ######## Error that was previously there!! N_t/P_t instead P_t/N_t ###
    ### mu1 = y[t-1,1] + theta[1] - log(1+theta[2]*N_t) - y[t-1,3]*N_t/P_t
    ######################################################################
    mu1 = y[t-1,1] + theta[1] - logprot(1+(exp(theta[1])-1)*N_t/theta[2]) - (theta[7]*P_t)/(1 + theta[7]*theta[8]*N_t)
    mu2 = y[t-1,2] + theta[4] - logprot((1+theta[5]*(exp(theta[4])-1)*P_t/N_t))
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
