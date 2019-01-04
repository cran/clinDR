"checkMonoEmax" <-
function(y,dose,parm,sigma2,nvec=rep(1,length(dose)),xbase=NULL,
			modelFun=emaxfun,trend='positive',binary=FALSE,logit=binary)
{
	if(any(is.na(c(y,dose))))stop('Missing data are not allowed')
	if(!is.matrix(parm))stop('parm must be a matrix')
	
	if(!missing(logit)){
		if(logit!=binary)stop('binary and logit contradict.  logit is deprecated, use binary')
		warning('logit is deprecated, use binary')	
	}
	
	nparm<-ncol(parm)
	nbase<-0
	if(!is.null(xbase)){
		if(any(nvec!=1))stop('nvec!=1 not allowed with covariates')
		if(is.vector(xbase))xbase<-as.matrix(xbase,ncol=1)
		if(!is.matrix(xbase))stop('xbase must be a matrix with rows matching y')
		if(nrow(xbase)!=length(y))stop('rows of xbase not equal to length of y')
		nbase<-ncol(xbase)
		if(nbase>=nparm)stop('Number of parameters specified must exceed columns of xbase')
		bparm<-parm[,(1+nparm-nbase):nparm]
		parm<-parm[,1:(nparm-nbase)]
		nparm<-nparm-nbase
	}	
	
	tol<-.Machine$double.eps
	
	
	nsim<-nrow(parm)
	ntot<-length(y)
	
	### assure ordering of doses
	dord<-order(dose)
	dose<-dose[dord]
	y<-y[dord]
	nvec<-nvec[dord]
	if(nbase)xbase<-xbase[dord,]
	
	dvec<-sort(unique(dose))
	ndose<-length(dvec)
	if(min(dvec)<tol)nld<-2 else nld<-1  ## omit lowest dose if pbo
	if(nld>ndose-1)return(Bpval=NA)     ## there must be at least 2 non-pbo doses
	
	ym<-numeric(ndose)
	nd<-numeric(ndose)
	for(i in 1:ndose){
		ym[i]<-weighted.mean(y[dose==dvec[i]],nvec[dose==dvec[i]])
		nd[i]<-sum(nvec[dose==dvec[i]])
	}
	
	if(trend=='positive'){
	    mdif<-max(ym[nld:(ndose-1)]-ym[ndose])
	}else{
	    mdif<-min(ym[nld:(ndose-1)]-ym[ndose])
	}
	
	if(!nbase){
		ypred<-modelFun(dvec,parm)
	}else{
		ypred<-modelFun(dose,parm)
		ypred<-ypred+bparm%*%t(xbase)
	}
	if(binary)ypred<-plogis(ypred)
	
	
	Bpval<-numeric(nsim)
	ymsim<-numeric(ndose)
	for(i in 1:nsim){
		if(!nbase){
	    if(!binary){
	        ymsim[1:ndose]<-rnorm(ndose,ypred[i,],sqrt(sigma2[i]/nd))
	    }else ymsim[1:ndose]<-rbinom(ndose,nd,ypred[i,])/nd
		}else{
		    if(!binary){
	        ysim<-rnorm(ntot,ypred[i,],sqrt(sigma2[i]))
		    }else ysim<-rbinom(ntot,1,ypred[i,])		
		    ymsim<-tapply(ysim,dose,mean)
		}
	    
    if(trend=='positive'){
        msimdif<-max(ymsim[nld:(ndose-1)]-ymsim[ndose])
        Bpval[i]<-1*(msimdif>=mdif)
    }else{
        msimdif<-min(ymsim[nld:(ndose-1)]-ymsim[ndose])
        Bpval[i]<-1*(msimdif<=mdif)
    }
	}
	
	return(Bpval=mean(Bpval))
}



