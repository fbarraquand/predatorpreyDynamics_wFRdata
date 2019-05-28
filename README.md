# Predator-prey dynamics with kill rate/functional response data 
Fits dynamic predator-prey models with (and without) auxillary kill rate (KR) data, which together with prey density N constitutes the functional response KR = G(N). 

* ``bayesian`` fits stochastic predator-prey models using JAGS (largely discrete-time Leslie-May models though Rosenzweig-MacArthur are considered as well). We fit models both with KR data and without KR data. 
* ``frequentist`` attempts maximum likelihood on the same models. RMA denotes the Rosenzweig-MacArthur model. 
* ``simulations`` contains specific simulations of mechanistic models (including some a little more complex than the fitted models).
* ``reports`` contain the progress report. The ``CHANGELOG.md`` file tracks important changes to the code during development. 
* ``figures_article`` tracks up-to-date, files that served to produce the figures of the corresponding article.
