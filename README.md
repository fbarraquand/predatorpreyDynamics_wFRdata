# Predator-prey dynamics with kill rate / functional response data 
Fits dynamic predator-prey models with (and without) auxiliary kill rate (KR) data, which together with prey density N constitutes the functional response KR = G(N). 

* ``bayesian`` fits stochastic predator-prey models using JAGS. These are argely discrete-time Leslie-May models, though Rosenzweig-MacArthur models are considered as well. We fit models to densities-only datasets and coupled datasets = (densities, kill rates). 
* ``frequentist`` attempts maximum likelihood estimation of the same models. RMA denotes the Rosenzweig-MacArthur model. 
* ``simulations`` contains specific simulations of mechanistic models (including some a little more complex than the fitted models).
* ``reports`` contain an early progress report. The ``CHANGELOG.md`` file tracks important changes to the code during development. 
* ``figures_article`` tracks the most important, up-to-date files that produced the figures of the corresponding article.
