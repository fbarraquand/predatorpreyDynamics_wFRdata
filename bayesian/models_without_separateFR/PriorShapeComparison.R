### Investigating the shape of the various priors considered

### Gamma priors most often used
curve(dgamma(x, shape = 0.01, rate = 0.01, log = FALSE),from=0,to=10) #shape rate param as in JAGS
curve(dgamma(x, shape = 0.1, rate = 0.1, log = FALSE),from=0,to=10)
curve(dgamma(x, shape = 0.001, rate = 0.001, log = FALSE),from=0,to=10)
# No peak anywhere -- perhaps difficult to then find a max in a posteriori distrib

mu = 2.5
sigma = 5
curve(dgamma(x, shape = mu^2/sigma, rate = mu/sigma, log = FALSE),from=0,to=10)
mu = 2.5
sigma = 2.5 
curve(dgamma(x, shape = mu^2/sigma, rate = mu/sigma, log = FALSE),from=0,to=10)
# Now a peak appears
curve(dgamma(x, shape = 5, rate = 1, log = FALSE),from=0,to=10)

curve(dlnorm(x, meanlog = 0, sdlog = 1, log = FALSE),from=0,to=10)