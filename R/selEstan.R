selEstan<-function(modType=3,binary=FALSE,pboAdj=FALSE){
	
	if(!binary & pboAdj)stop('placebo-adjusted binary data not implemented')

	if(isTRUE(grep("64",Sys.getenv("R_ARCH"))>0)){
		emod<-'comp64'
	}else emod<-'comp32'
	
	
	if(binary){
		emod<-file.path(emod,'estanbin')
		emod<-paste(emod,modType,sep='')
	####################
	#### continuous data
	}else{
		emod<-file.path(emod,'estancont')
		if(!pboAdj){
			emod<-paste(emod,modType,sep='')
		}else{  ### pbo excluded
			emod<-paste(emod,modType,'NoPbo',sep='')		
		}	
	}
	emod<-paste(emod,'.rds',sep='')
	emod<-system.file(package="clinDR", "models", emod)
	
	if(file.access(emod,mode=0)<0)stop(paste('The compiled rstan model',
										'could not be accessed.  You must',
										'run compileStanModels once before using',
										'the Bayesian functions'))	
	
	estan<-readRDS(emod)	
	
	if(class(estan)!='stanmodel')stop('unable to create estan model')		
	
	return(estan)
}
		