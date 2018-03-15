# predatorpreyDynamics_wFRdata
Fits dynamic predator-prey models with (and without) auxillary functional response (FR) data

* ``bayesian`` fits predator-prey models using JAGS (largely discrete-time Leslie-May models). With FR data and without FR data. 
* ``frequentist`` attempts maximum likelihood on the same models. RMA denote the Rosenzweig-MacArthur model. 
* ``simulations`` contains specific simulations of mechanistic models (including some a little more complex than the fitted models)

### Progress report FB 25/11/2017

* It seems that the predator-prey model with noisy FR is identifiable - using JAGS - provided data on (N,P,FR(N,P)) when the simulation is made with a noisy FR. In contrast, a predator-prey model with deterministic FR fitted without using the FR data does not provide reliable estimates of the functional response. **This suggests that FR data is paramount to identifiability of predator-prey models**

* This could be simply the noisyness of the functional response that impedes the inference. However, further investigation onto the predator-prey models with deterministic FR suggests that it may not be the case. Indeed, fitting a predator-prey model (in JAGS) with deterministic FR to data simulated with a deterministic FR provides flawed estimates of maximum intake. Hence there might be identifiability issues even with the simplest predator-prey models. I previously simulated host-parasitoid models and found no such issues (again in JAGS), so the problem might be model-structure-dependent (i.e., depending on the precise functional forms used for the attack rate). **There is a need for a detailed investigation of the identifiability properties of predator-prey models**. 

* Identifiability might be looked at using the tools developed by Catchpole and Morgan (e.g., exhaustive summaries) :: TODO

* I tried to investigate this identifiability in a simple frequentist context, by producing likelihood profiles and optimization to find min/max. So far it has not been very fruitful, see ``frequentist``, many errors of optim and difficult to do likelihood profiles that make sense. 

### Progress report FB 12/03/2018

* **But this was due to an error in the code, N/P instead of P/N** (damn it!). The profiles are now much easier to produce and the optimization code to find the max likelihood frequently produces warnings but works. 

* Next: fits predator-prey with deterministic FR to predator-prey with noisy FR data? Or check what the discrepancy in JAGS might be due to (-> check parameters used). Also improve the contour plots (some have shitty axes). 

### Progress report FB 14/03/2018

* The discrepancy between ``JAGS`` and ``optim`` indeed stems from differences in parameter values (was too much noise, too low C). 
* We need to now look at the Hessian and make sure the model is identifiable (check Olivier's old paper as well). 
