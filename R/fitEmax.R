"fitEmax"<-
function(y,dose,iparm,modType=3,prot=rep(1,length(y)),count=rep(1,length(y)),
		 binary=FALSE,diagnostics=TRUE,msSat=NULL,pboAdj=FALSE,optObj=TRUE){

	tol<-.Machine$double.eps
	protorig<-factor(prot)             #ensure that prot is a factor
	prot<-as.numeric(protorig)         #convert to 1,2,3,....
    nprot<-max(prot)
	dlev<-sort(unique(dose))

	### check input consistency
	if(isTRUE(pboAdj && binary))stop('PBO adjustment not available with binary data')

	### check for missing data
	if(any(is.na(y)))stop('Missing y data should be imputed or deleted')

	### check input consistency
	lengthvec<-c(length(y),length(dose),length(prot),length(count))
	if(any(abs(lengthvec-lengthvec[1])>tol))stop('Length of y and dose must be equal')

    if(missing(iparm)){
        iparm<-startEmax(y,dose,count=count,modType=modType,binary=binary)
        ### replicate intercept for multiple protocols
    }else{
        if(modType==3){
			if(length(iparm)!=3)stop('Incorrect number of starting parameters')
			names(iparm)<-c('led50','emax','e0')
		}else if(modType==4){
			if(length(iparm)!=4)stop('Incorrect number of starting parameters')
			names(iparm)<-names(iparm)<-c('led50','lambda','emax','e0')
		}
    }
	### replicate intercept for multiple protocols
    svec<-c(iparm,rep(iparm[modType],nprot-1))
    for(i in modType:(modType+nprot-1))names(svec)[i]<-paste('e',i+1-modType,sep='')

	###############################################################################
	#### binary data
	if(binary){
		### y must be coded 0/1
		if(any(y!=0 & y!=1))stop("y must be 0/1")


		if(diagnostics){
			print.level<-2
			fit<-nlm(f=nllogis, p=svec, hessian=TRUE, print.level=print.level,
				y=y,dose=dose,prot=prot,count=count)
		}else{
			print.level<-0
			fit<-suppressWarnings(nlm(f=nllogis, p=svec, hessian=TRUE, print.level=print.level,
				y=y,dose=dose,prot=prot,count=count))
		}

		if(fit$code>2){
			if(diagnostics)warning('nlm failed to converge to a stable solution')
			return(NULL)
		}else {
			#### var-cov computation
			vcfail<-FALSE
			ev<-try(eigen(fit$hessian,symmetric=TRUE),silent=!diagnostics)
			if(class(ev)=='try-error'){
			   vcfail<-TRUE
			}else{
				evals<-ev$values
				if(any(evals<=0)){vcfail<-TRUE
				}else{
					evec<-ev$vectors
					vc<-evec%*%diag(1/evals)%*%t(evec)
				}
				if(vcfail){
					if(diagnostics)warning('var-cov matrix was not positive definite')
					return(NULL)
				}
			}
		}
	   ### negative loglihood likelihoods
		nllMod<- fit$minimum
		nllSat<-0
		for(i in 1:nprot){
			for(d in dlev){
				npd<-sum(dose==d & prot==i)
				if(npd>0){
					npd0<-sum(count[y==0 & prot==i & dose==d])
					npd1<-sum(count[y==1 & prot==i & dose==d])
					phat<-(npd1+.5)/(npd0+npd1+1)   ## add cont correction
					nllSat<-nllSat-npd0*log(1-phat)-npd1*log(phat)
				}
			}
		}

		totdp<-nrow(unique(cbind(prot,dose)))
		if(totdp-(nprot+modType-1)>0){
			dfMod<-totdp-(nprot+modType-1)
			gofTest<- pchisq(2*(nllMod-nllSat),dfMod,lower.tail=FALSE)
		}else{
			gofTest<-NA
			dfMod<-NA
			if(diagnostics)warning('No df for testing gof.')
		}

		rownames(vc)<-names(svec)
		colnames(vc)<-names(svec)
		names(fit$estimate)<-names(svec)
		fitl<-list(estimate=fit$estimate,vc=vc)
		if(optObj){optObj<-fit
		}else optObj<-NULL

        fitl<-list(fit=fitl,y=y,dose=dose,modType=modType,prot=protorig,count=count,
				  nprot=nprot,binary=binary,pboAdj=FALSE,
				  residSD=NULL,gofTest=gofTest,nll=c(nllMod=nllMod,nllSat=nllSat),
				  df=c(dfMod=dfMod,dfSat=0),optObj=optObj)
		class(fitl)<-"fitEmax"
		return(fitl)


	###############################################################################
	#### continuous data
	}else{
    ### form strata indicators
    for(i in 1:nprot)assign(paste("I",i,sep=""),1*(prot==i))

    intCode<-'+e1*I1'
    if(nprot>1){
        for(i in 2:nprot)intCode<-paste(intCode,' + e',i,'*I',i,sep='')
    }
		if(pboAdj){
			intCode<-''
			svec<-svec[1:(modType-1)]
		}
 
    if(modType==3){chform<-paste('y ~','(emax*dose)/(dose + exp(led50))',intCode,sep='') 
    }else{chform<-paste('y ~','(emax*dose^lambda)/(dose^lambda + (exp(led50))^lambda)',intCode,sep='') 
    }

    nlsW.fit<-try(nls(as.formula(chform),
                           start=svec, weights=count,
                           control = nls.control(
                           maxiter = 100),trace=diagnostics,na.action=na.omit),
									silent=!diagnostics) 

    if(class(nlsW.fit)!='nls'){
            
        intCode<-',e1=I1'
        if(nprot>1){
            for(i in 2:nprot)intCode<-paste(intCode,',e',i,'=I',i,sep='')
        }

        if(modType==3){
					if(!pboAdj)chform<-paste('y ~ cbind(emax=dose/(dose + exp(led50))',intCode,')',sep='')
   				if(pboAdj)chform<-paste('y ~ cbind(emax=dose/(dose + exp(led50))',')',sep='')
            svsub<-svec[1]
        }else{
					if(!pboAdj)chform<-paste('y ~ cbind(emax=(dose^lambda)/(dose^lambda + (exp(led50))^lambda)',
											 intCode,')',sep='')
					if(pboAdj)chform<-paste('y ~ cbind(emax=(dose^lambda)/(dose^lambda + (exp(led50))^lambda)', ')',
											sep='')
		                       svsub<-svec[1:2]
        } 
        nlsW.fit <- try(nls(as.formula(chform), 
						            start = c(svsub), weights=count, control = nls.control(
    							      maxiter = 100), algorithm="plinear", trace = diagnostics),
										silent= !diagnostics)
    }
	badvar<-FALSE
	if(class(nlsW.fit)=='nls'){
		vc<-try(vcov(nlsW.fit),silent=!diagnostics)
		if(class(vc)=='try-error'){ badvar<-TRUE
		}else if(any(is.na(vc)))badvar<-TRUE
	}

    if(class(nlsW.fit)!='nls'){
        if(diagnostics)warning('nls failed to converge to a stable solution')
        return(NULL)
	}else if(badvar){
        if(diagnostics)warning('var-cov matrix was not positive definite')
        return(NULL)
    }else{
		totdp<-nrow(unique(cbind(prot,dose)))
		if(any(count>1) & is.null(msSat)){  ### msSat not supplied with aggregate data
			gofTest<-NA
			residSD<-summary(nlsW.fit)$sigma
			if(pboAdj){
				dfMod<-length(y)-(nprot+modType-2)
			}else dfMod<-length(y)-(nprot+modType-1)
			msMod<-residSD^2
			nllMod<-msMod*dfMod
			dfSat<-NA
			nllSat<-NA
		}else if(any(count>1) & !is.null(msSat)){  ## user supplied msSat 
			### the ss are decomposition into within-saturated, saturated 
			### within protocol, protocol within common intercept
			dfSat<-sum(count)-totdp
			ssSat<-dfSat*msSat

			### compute ss from the emax model because the nls weighted ss is not 
			### appropriate if there are replicated dose group means within a protocol
			ssMod<-0
			for(i in 1:nprot){
				for(d in dlev){
					if(sum(dose==d & prot==i)>0){
						ind<-(dose==d & prot==i)
						ecoef<-coef(nlsW.fit)
						if(pboAdj){
							ecoef<-c(ecoef,0)
						}else{
							ecoef[modType]<-ecoef[modType+i-1]
							ecoef<-ecoef[1:modType]
						}
						ehat<-emaxfun(d,ecoef) 
						ym<-weighted.mean(y[ind],count[ind])
						nm<-sum(count[ind])
						ssMod<-ssMod+nm*(ym-ehat)^2
					}
				}
			}
			### decomposition into within-saturated model, and saturated model within emax model
			if(pboAdj){
				dfMod<-totdp-(nprot+modType-2)
			}else dfMod<-totdp-(nprot+modType-1)
			if(dfMod>0){
				msMod<-ssMod/dfMod
				gofTest<-pf(msMod/msSat,
					 dfMod,dfSat,lower.tail=FALSE)
				ssMod<-ssMod+ssSat  ### improve error estimate using ss from sat also
				dfMod<-sum(count)-(nprot+modType-(1+pboAdj*1))
				msMod<-ssMod/dfMod
				residSD<-sqrt(msMod)
				nllMod<-ssMod
				nllSat<-ssSat
			}else{
				gofTest<-NA
				dfMod<-NA
				dfSat<-NA
				msMod<-NA
				residSD<-NA
				nllMod<-NA
				nllSat<-NA
				if(diagnostics)warning('No df for GOF testing')
			}
		}else{   ### patient-level data are available
			### compute ss from saturated model
			ssSat<-0
			for(i in 1:nprot){
				for(d in dlev){
					if(sum(dose==d & prot==i)>0){
						ind<-(dose==d & prot==i)
						ymsat<-mean(y[ind])  
						ssSat<-ssSat+sum((y[ind]-ymsat)^2)
					}
				}
			}
			dfSat<-sum(count)-totdp
			if(pboAdj){
				dfMod<-sum(count)-(nprot+modType-2)
			}else dfMod<-sum(count)-(nprot+modType-1)
			msMod<-(summary(nlsW.fit)$sigma)^2
			if((dfSat>0) & (dfMod-dfSat>0)){
				msSat<-ssSat/dfSat
				gofTest<-pf((dfMod*msMod-dfSat*msSat)/((dfMod-dfSat)*msSat),
					 dfMod-dfSat,dfSat,lower.tail=FALSE)
				nllMod<-dfMod*msMod
				nllSat<-dfSat*msSat
			}else{
				gofTest<-NA
				dfMod<-NA
				msMod<-NA
				msSat<-NA
				nllMod<-NA
				nllSat<-NA
				if(diagnostics)warning('No df for GOF testing')
			}
			residSD<-sqrt(msMod)
		}
	}

      if(any(count>1) & !is.null(msSat)){  ### adjust vc to use most df as possible
			vc<-vc*msMod/((summary(nlsW.fit)$sigma)^2)
		}

		fit<-list(estimate=coef(nlsW.fit),vc=vc)
		if(optObj){optObj<-nlsW.fit
		}else optObj<-NULL
        fit<-list(fit=fit,y=y,dose=dose,modType=modType,prot=protorig,count=count,
				  binary=binary,pboAdj=pboAdj,
				  residSD=residSD,gofTest=gofTest,nll=c(nllMod=nllMod,nllSat=nllSat),
				  df=c(dfMod=dfMod,dfSat=dfSat),optObj=optObj)
        class(fit)<-'fitEmax'
        return(fit)
    }
}


