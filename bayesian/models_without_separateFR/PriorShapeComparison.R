### Investigating the shape of the various priors considered

### Gamma priors most often used
curve(dgamma(x, shape = 0.01, rate = 0.01, log = FALSE),from=0,to=10) #shape rate param as in JAGS
curve(dgamma(x, shape = 0.1, rate = 0.1, log = FALSE),from=0,to=10)
curve(dgamma(x, shape = 0.001, rate = 0.001, log = FALSE),from=0,to=10)
# No peak anywhere -- perhaps difficult to then find a max in a posteriori distrib
### In JAGS dgamma(shape, rate too) p. 45 of JAGS Version 4.3.0 user manual 

mu = 2.5
sigma = 5
curve(dgamma(x, shape = mu^2/sigma, rate = mu/sigma, log = FALSE),from=0,to=10)
mu = 2.5
sigma = 2.5 
curve(dgamma(x, shape = mu^2/sigma, rate = mu/sigma, log = FALSE),from=0,to=10)
# Now a peak appears
curve(dgamma(x, shape = 5, rate = 1, log = FALSE),from=0,to=10)

curve(dlnorm(x, meanlog = 0, sdlog = 1, log = FALSE),from=0,to=10)

### Multipanel to compare 

### Priors gamma
par(mfrow=c(2,2))
curve(dgamma(x, shape = 0.01, rate = 0.01),from=0,to=50)
abline(v=2.5,col="red")
curve(dgamma(x, shape = 0.2, rate = 0.2),from=0,to=50)
abline(v=2.5,col="red")
curve(dgamma(x, shape = 2, rate = 0.1),from=0,to=50)
abline(v=2.5,col="red")
curve(dgamma(x, shape = 2, rate = 1),from=0,to=50)
abline(v=2.5,col="red")

### Priors log-normal
curve(dlnorm(x, 1, 100),from=0,to=50)
abline(v=2.5,col="red")
curve(dlnorm(x, 1, 10),from=0,to=50)
abline(v=2.5,col="red")
curve(dlnorm(x, 1, 2),from=0,to=50)
abline(v=2.5,col="red")
curve(dlnorm(x, 1, 0.5),from=0,to=50)
abline(v=2.5,col="red")
