selEstan<-function(emod=c('basemodel.rds','mrmodel.rds')){
	

	emod<-match.arg(emod)
	
	if(isTRUE(grep("64",Sys.getenv("R_ARCH"))>0)){
		ebase<-'comp64'
	}else ebase<-'comp32'
	
	emod<-file.path(ebase,emod)
	
	emod<-system.file(package="clinDR", "models", emod)
	
	if(file.access(emod,mode=0)<0)stop(paste('The compiled rstan model',
										'could not be accessed.  You must',
										'run compileStanModels once before using',
										'the Bayesian functions'))	
	
	estan<-readRDS(emod)	
	
	if(!inherits(estan,'stanmodel'))stop('unable to create estan model')		
	
	return(estan)
}
		