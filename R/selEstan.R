selEstan<-function(){
	

	if(isTRUE(grep("64",Sys.getenv("R_ARCH"))>0)){
		emod<-'comp64'
	}else emod<-'comp32'
	
	emod<-file.path(emod,'basemodel.rds')
	
	emod<-system.file(package="clinDR", "models", emod)
	
	if(file.access(emod,mode=0)<0)stop(paste('The compiled rstan model',
										'could not be accessed.  You must',
										'run compileStanModels once before using',
										'the Bayesian functions'))	
	
	estan<-readRDS(emod)	
	
	if(class(estan)!='stanmodel')stop('unable to create estan model')		
	
	return(estan)
}
		