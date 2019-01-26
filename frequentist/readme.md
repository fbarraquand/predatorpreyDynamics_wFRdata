I am more and more convinced that we should not compute identifiability properties (through the hessian) at the MLE:
* this makes sense only if the MLE (for that particular dataset) or endpoint of the algo is very close to the true value
* if not, we compte whether we can attain a point which we may not want to attain anyway -> makes no sense? 

The expected FIM (e.g., over 100 time series of length 100?) may be more useful for some cases. 
