
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

* Long time series (T=1000, painfully long in JAGS) reveal that the identifiability (identification rather?) improves for the predator-prey model without the FR. Hence perhaps all the models are identifiable in a structural sense (infinite TS) but adding the FR info improves greatly convergence whenever there is limited data. 

* We need to now look at the Hessian and make sure the model is identifiable (check Olivier's old paper as well). Remember I saw a ridge in the likelihood regarding (r,K) as previously highlighted by Polanski (on the other hand, won't this strongly depend on parameter space -- I had no such problems in Barraquand et al. 2014). Morality: Use frequentist code for identifiability, the point in parameter space previously used + a point in the predator-prey noisy limit cycle region of parameter space, as well as T=1000 for near perfect conditions first. 

### Progress report FB 15/03/2018

* Using the point in parameter space previously used, and T=1000, no identifiability problems even with the model without noisy functional response/with deterministic FR (see ``frequentist/Hessian.R``). As evaluated through the rank of the Hessian. 

* I'm unclear as to whether the Hessian matrices are always positive definite (for T=1000), which is important as well for concavity! (some eigenvalues can be negative, but always when these are close to zero). NB Perhaps I can look at this with a Cholesky decomposition as well. 

* There can be change in several orders of magnitude in the eigenvalues (e.g. from $10^3$ to $10^{-2}$). Does this mean I should consider everyting below 2 orders of magnitude from the max to zero? In which case there would not be full rank matrices... Or switch to symbolic computation, but if I can avoid that I'd rather do so... 

* NB Important: the fact that the "basic" parameter set is identifiable may be due to the fact that it exhibits fluctuation when noise is added, which in turn might be due to a pair of complex eigenvalues -- dampened oscillations to equilibrium. To check. Although we have r=2 so the basic prey model might exhibit dampening oscillations as well. A Jacobian computation might be really needed!

### Progress report FB 14/01/2019

* Identifiability assessed for the model and parameter correlations highlighed, for both bayesian and frequentist versions of the code. Code in ``Hessian.R`` notably improved. Added (over the last few weeks) more diagnostics of posterior distributions in the bayesian ``longTS_T=1000`` folder. Profile log-likelihood code improved as well. 

* Jacobian-based stability analysis of the deterministic model allows to classify noisy cycles vs perturbed focus parameter sets. In ``analytical``

* A delayed FR, depending on densities at the previous time step, rather than depending on current density was found in the simulation code in the ``bayesian`` folder (hence model simulated was slightly different than fitted). The ``simulations`` folder was corrected accordingly. In the ``frequentist`` folder, we had a delayed FR that was both simulated and fitted. Therefore all previous results (e.g., on identifiability and parameter correlations) may hold, even though the final model that makes sense biologically and that we want to analyse has a FR<sub>t</sub> that depends on N<sub>t</sub> and P<sub>t</sub>. Suprisingly, the delayed FR model was not much more instable than the nondelayed FR model, which probably explain why this error was not found until now. Now corrected. 

* The reparameterization was corrected (after looking up the equations in the math derivations). It looks like the (r,K) and (s,Q) correlations are now mostly removed, but the (C,D) or (a,h) correlations are more "resistant" to reparameterization...
