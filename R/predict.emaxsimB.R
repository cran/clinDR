"predict.emaxsimB" <-
function(object,dose, dref=0, ...){
	warning(paste("\nPredicted values for doses included in the study\n",
								"can be obtained from the fitpredv and associated\n",
								"standard errors sepredv, sedifv, stored in the emaxsimB\n",
								"object.  For other doses, emaxsimB must be re-run and\n",
								"the predicted values computed using custom code.\n",sep=''))	
  return(invisible())
}

