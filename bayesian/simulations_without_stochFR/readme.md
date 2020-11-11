## Investigation of the poor performance of the model without explicit FR data *even when fitted to data simulated under this model* 

* ``woutFR`` refers here to a model without explicit FR data & an explicitly stochastic "ground truth" model for the functional response. 

Could we blame the following factors for the low performance? 
* Environmental noise level (perhaps much better at low noise levels) -> NO
* Gamma priors on the functional response parameters -- would say, lognormal instead change results? ->  Not really
* Time series - is 100 in fact too low? -> NO
* Boils down to *whether the deterministic data describes a limit cycle or not*. Very good performance, even for T=100 of the noisy FR data. 
* Precision of the priors is key for identifiability of the perturbed fixed point. 

Here slightly more precise priors re-allow identifiability (``PredPreyStoch_woutFR _sigmaLow_CHigh_dgammaLessVague.R``, dgamma(2,1)) in the perturbed fixed point case. Unlike what occurs when we have a separate, stochastic FR (i.e., including envtal noise on the FR), where slightly more precise priors do not re-allow identifiability (``../TS=100/PredPreyStoch_analysisBis_smallNoise_morePrecisePriors.R``, where the prior information still dominates over the likelihood). 







