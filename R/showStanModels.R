showStanModels<-function(emod=c('basemodel.stan','mrmodel.stan')){
	

	emod<-match.arg(emod)
	emod<-system.file(package="clinDR", "models", emod)
	
	if(file.access(emod,mode=0)<0)stop(paste('The rstan model',
										'could not be accessed.  You must',
										'run compileStanModels once before using',
										'the Bayesian functions'))	
	
	file.show(emod)	
	
	return(invisible(emod))
}
		