"predict.fitEmaxB"<-
function(object,dosevec,clev=0.9,int=1,dref=0, xvec=NULL, ...){
	
	if(missing(dosevec))stop('dosevec must be specified')

	modType<-object$modType
	binary<-object$binary
	pboAdj<-object$pboAdj
	prot<-object$prot
	nprot<-max(as.numeric(prot))
	nbase<-object$nbase
	if(nbase)xbase<-object$xbase[prot==int,]
	
	if(! int%in%(c(1:nprot)))stop('The intercept specification is invalid')

	clevup<- 0.5+clev/2
	clevlow<-0.5-clev/2
	
	parms<-coef(object)
	nparm<-ncol(parms)
	if(nbase>0)bparms<-parms[,(1+nparm-nbase):nparm]
	parms<-parms[,1:(nparm-nbase)]
	if(pboAdj){
		parms<-cbind(parms,rep(0,nrow(parms)))		
	}else parms<-coef(object)[,c(1:(modType-1),modType+(int-1))]
	if(!binary) sigsim<-sigma(object) else sigsim<-NULL
	
	predout<- emaxfun(dosevec,parms)
	predref<- as.vector(emaxfun(dref,parms))
	if(nbase>0){
		if(!is.null(xvec)){
			xvec<-matrix(xvec,ncol=1)
			bcont<-bparms%*%xvec
		}else bcont<-bparms%*%t(xbase)
	}
	if(binary){
		if(nbase>0){
			predref<-apply(plogis(predref+bcont),1,mean)	
			for(i in 1:length(dosevec)){
				predout[,i]<-apply(plogis(predout[,i]+bcont),1,mean)	
			}
		}else{
			predout<-plogis(predout)
			predref<-plogis(predref)
		}
	}else{
		if(nbase>0)badd<-apply(bcont,1,mean) else badd<-0
		predout<-predout+badd
		predref<-predref+badd
	}
	pred<-  apply(predout,2,mean) 
	predMed<-  apply(predout,2,median) 
	sepred<- sqrt( apply(predout,2,var) )
	lb<- apply(predout,2,quantile,probs=clevlow)
	ub<- apply(predout,2,quantile,probs=clevup)
	names(pred)<-dosevec
	names(predMed)<-dosevec
	names(sepred)<-dosevec
	names(lb)<-dosevec
	names(ub)<-dosevec
	
	fitdif<-apply(predout-predref,2,mean)
	fitdifMed<-apply(predout-predref,2,median)
	sedif<- sqrt( apply(predout-predref,2,var) )
	lbdif<- apply(predout-predref,2,quantile,probs=clevlow)
	ubdif<- apply(predout-predref,2,quantile,probs=clevup)
	names(fitdif)<-dosevec
	names(fitdifMed)<-dosevec
	names(sedif)<-dosevec
	names(lbdif)<-dosevec
	names(ubdif)<-dosevec
	
	colnames(predout)<-dosevec
	
	return(list(pred=pred,predMed=predMed,lb=lb,ub=ub,se=sepred,
		 fitdif=fitdif,fitdifMed=fitdifMed,lbdif=lbdif,
		 ubdif=ubdif,sedif=sedif,simResp=predout,sigsim=sigsim))
}


