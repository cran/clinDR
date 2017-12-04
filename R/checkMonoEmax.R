"checkMonoEmax" <-
function(y,dose,parm,sigma2,nvec=rep(1,length(dose)),
			modelFun=emaxfun,trend='positive',logit=FALSE)
	{
    if(any(is.na(c(y,dose))))stop('Missing data is not allowed')
  
    tol<-.Machine$double.eps
  
  
    nsim<-nrow(parm)

    ### assure ordering of doses
    dord<-order(dose)
    dose<-dose[dord]
    y<-y[dord]
    nvec<-nvec[dord]

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

    ypred<-modelFun(dvec,parm)
    if(logit)ypred<-plogis(ypred)


    Bpval<-numeric(nsim)
    ymsim<-numeric(ndose)
    for(i in 1:nsim){
        if(!logit){
            ymsim[1:ndose]<-rnorm(ndose,ypred[i,],sqrt(sigma2[i]/nd))
        }else ymsim[1:ndose]<-rbinom(ndose,nd,ypred[i,])/nd
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



