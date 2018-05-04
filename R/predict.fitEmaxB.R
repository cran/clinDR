"predict.fitEmaxB"<-
function(object,dosevec,clev=0.9,int=1,dref=0, ...){
	
	if(missing(dosevec))stop('dosevec must be specified')

	modType<-object$modType
	binary<-object$binary
	pboAdj<-object$pboAdj
	nprot<-max(as.numeric(object$prot))

	if(! int%in%(c(1:nprot)))stop('The intercept specification is invalid')

	clevup<- 0.5+clev/2
	clevlow<-0.5-clev/2
	
	if(!binary) sigsim<-sigma(object) else sigsim<-NULL
	if(pboAdj){
		parms<-coef(object)
		parms<-cbind(parms,rep(0,nrow(parms)))		
	}else parms<-coef(object)[,c(1:(modType-1),modType+(int-1))]
	
	predout<- emaxfun(dosevec,parms)
	predref<- as.vector(emaxfun(dref,parms))
	if(binary){
		predout<-plogis(predout)
		predref<-plogis(predref)
	}
	pred<-  apply(predout,2,mean) 
	sepred<- sqrt( apply(predout,2,var) )
	lb<- apply(predout,2,quantile,probs=clevlow)
	ub<- apply(predout,2,quantile,probs=clevup)
	names(pred)<-dosevec
	names(sepred)<-dosevec
	names(lb)<-dosevec
	names(ub)<-dosevec
	
	fitdif<-apply(predout-predref,2,mean)
	sedif<- sqrt( apply(predout-predref,2,var) )
	lbdif<- apply(predout-predref,2,quantile,probs=clevlow)
	ubdif<- apply(predout-predref,2,quantile,probs=clevup)
	names(fitdif)<-dosevec
	names(sedif)<-dosevec
	names(lbdif)<-dosevec
	names(ubdif)<-dosevec
	
	colnames(predout)<-dosevec
	
	return(list(pred=pred,lb=lb,ub=ub,se=sepred,
		 fitdif=fitdif,lbdif=lbdif,
		 ubdif=ubdif,sedif=sedif,simResp=predout,sigsim=sigsim))
}


