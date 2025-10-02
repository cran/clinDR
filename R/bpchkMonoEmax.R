"bpchkMonoEmax" <-
function(x, trend='positive', protSel=1)
{
	if(!inherits(x,'fitEmaxB'))stop('first input must be fitEmaxB output object')
	
	prot<-as.numeric(factor(x$prot))         #convert to 1,2,3,....
  nprot<-max(prot)
	if(protSel<0 || protSel>nprot)stop('invalid protSel specified')
  
	y<-x$y[prot==protSel]
	dose<-x$dose[prot==protSel]
	parm<-coef(x)
	modType<-x$modType
	nvec<-x$count
	nbase<-0
	xbase<-x$xbase
	if(length(xbase)){
		nbase<-ncol(xbase)
		xbase<-xbase[prot==protSel,]
	}
	binary<-x$binary
	dimFit<-x$dimFit
	vcest<-x$vcest
	if(modType==3) parm<-cbind(parm[,1],rep(1,nrow(parm)),parm[,2:ncol(parm)])
	if(nbase){
		parm<-cbind(parm[,1:3],parm[,3+protSel],
								parm[,(4+nprot):(3+nprot+nbase)])
	}else parm<-cbind(parm[,1:3],parm[,3+protSel])
	nparm<-ncol(parm)
	
	if(!binary && !dimFit)sigma2<-(sigma(x)^2)else sigma2<-NULL
	
	if(nbase){
		bparm<-parm[,(1+nparm-nbase):nparm]
		parm<-parm[,1:(nparm-nbase)]
		nparm<-nparm-nbase
	}	
	
	tol<-sqrt(.Machine$double.eps)
	
	
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
		ypred<-emaxfun(dvec,parm)
	}else{
		ypred<-emaxfun(dose,parm)
		ypred<-ypred+bparm%*%t(xbase)
	}
	if(binary && !dimFit)ypred<-plogis(ypred)
	
	
	Bpval<-numeric(nsim)
	ymsim<-numeric(ndose)
	for(i in 1:nsim){
		if(!nbase){
			if(dimFit){
					ymsim[1:ndose]<-ypred[i,]+rmvnorm(1,rep(0,dimFit),vcest)
	    }else if(!binary){
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



