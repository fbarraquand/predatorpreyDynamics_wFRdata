## Things to investigate to explain (very) poor performance of the model without explicit FR data *even when fitted to data simulated under this model* 

* ``woutFR`` refers here to model without explicit FR data/stochastic FR model. 

Could we blame those factors for the low performance?  
* Environmental noise level (perhaps much better at low noise levels) -> NO
* Gamma priors on the functional response parameters -- would say, lognormal change results? ->  Not really
* Time series - is 100 in fact too low? -> NO
* Boils down to *whether the deterministic data describes a limit cycle or not*. Very good performance, even for T=100 of the noisy FR data. Would be interesting to compare to the case where T=100 but we simulate a stochastic FR, for the noisy limit cycle parameters.  

(NB perhaps in some cases very close to the limit cycle but "under" the model works)




