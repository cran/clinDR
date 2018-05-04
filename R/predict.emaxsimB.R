"predict.emaxsimB" <-
function(object,dose, dref=0, ...){
	warning(paste("Predicted values for doses included in the study",
								"can be obtained from the fitpredv and associated",
								"standard errors sepredv, sedifv, stored in the emaxsimB",
								"object.  For other doses,",
								"emaxsimB must be re-run and the predicted values",
								"computed using custom code."))	
  return(invisible())
}

