# predatorpreyDynamics_wFRdata
Fits dynamic predator-prey models with (and without) auxillary functional response (FR) data

* ``bayesian`` fits predator-prey models using JAGS (largely discrete-time Leslie-May models). With FR data and without FR data. 
* ``frequentist`` attempts maximum likelihood on the same models. RMA denote the Rosenzweig-MacArthur model. 
* ``simulations`` contains specific simulations of mechanistic models (including some a little more complex than the fitted models)

### Progress report FB 25/11/2017

* It seems that the predator-prey model with noisy FR is identifiable - using JAGS - provided data on (N,P,FR(N,P)) when the simulation is made with a noisy FR. In contrast, a predator-prey model with deterministic FR fitted without using the FR data does not provide reliable estimates of the functional response. **This suggests that FR data is paramount to identifiability of predator-prey models**

* This could be simply the noisyness of the functional response that impedes the inference. However, further investigation onto the predator-prey models with deterministic FR suggests that it may not be the case. Indeed, fitting a predator-prey model with deterministic FR to data simulated with a deterministic FR provides flawed estimates of maximum intake. Hence there might be identifiability issues even with the simplest predator-prey models. I previously simulated host-parasitoid models and found no such issues (again in JAGS), so the problem might be model-structure-dependent (i.e., depending on the precise functional forms used for the attack rate). **There is a need for a detailed investigation of the identifiability properties of predator-prey models**. 

* Identifiability might be looked at using the tools developed by Catchpole and Morgan (e.g., exhaustive summaries) :: TODO

* I tried to investigate this simply by producing likelihood profiles. So far it has been fruitful, see ``frequentist``. The profile are difficult to produce and the optimization code to find the max likelihood frequently fails. 
