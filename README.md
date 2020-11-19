[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4281411.svg)](https://doi.org/10.5281/zenodo.4281411)

# Fitting stochastic predator-prey models using both population density and kill rate data
The code provided here fits dynamic predator-prey models to data, contrasting cases where we have access to auxiliary kill rate data vs cases where we do not. The predator kill rate (KR), together with prey density N, constitutes the functional response KR = g(N). Here we make the assumption that the kill rate is subjected to process noise (i.e., random temporal variability) rather than observational error. Companion code to [our paper with the same title](https://arxiv.org/abs/1904.02145)

* ``bayesian`` fits the stochastic predator-prey models using JAGS. These are largely discrete-time Leslie-May models, though Rosenzweig-MacArthur models are considered as well. We fit models to densities-only datasets and coupled datasets = (densities, kill rates). 
* ``frequentist`` attempts maximum likelihood estimation of the same models. RMA denotes the Rosenzweig-MacArthur model. 
* ``simulations`` contains the specific simulations of mechanistic models used for all model fitting.
* ``reports`` contain an early progress report. The ``CHANGELOG.md`` file tracks important changes to the code during development. 
* ``figures_article`` tracks the most important, up-to-date files that produced the figures of the corresponding article.
