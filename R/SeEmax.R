"SeEmax" <-
function(fit,doselev,modType,dref=0,nbase=0,x=NULL,
				 binary=FALSE, clev=0.9){

	if(! modType %in% c(3,4))stop("modType must be 3 or 4")

	if(class(fit)=='nls'){
		if(!is.null(x))stop(paste('Baseline covariates and multiple',
			'protocols are not supported with objects of class NLS.  Use list input.'))
		parm<-coef(fit)
		vfcov <- vcov(fit)[1:modType,1:modType]
	}else{
		parm<-fit[[1]]
		vfcov<-fit[[2]]
	}
	nparm<-length(parm)
	if(nparm!=(modType+nbase))stop('Incorrect number of parameters')
	if(any(dim(vfcov)!=nparm))stop('Variance-covariance matrix does not match parmeters')
	if(nbase){
		if (is.null(x))stop('x must be specified when nbase>0')
		if(is.vector(x)){
			if(length(x)!=nbase)stop('number of covariates does not match nbase')
			nx<-1
			x<-matrix(x,nrow=1)
		}else if(is.matrix(x)){
			if(ncol(x)!=nbase)stop('number of columns in x does not match nbase')
			nx<-nrow(x)
		}
		bparm<-parm[(1+modType):nparm]
	}

	ed50<-exp(parm[1])
	ndose<-length(doselev)

	predout <- emaxfun(doselev, parm[1:modType])
	predR <- emaxfun(dref, parm[1:modType])
	
	if(nbase){
		if(!binary)xmean<-apply(x,2,mean)
		covc<-as.vector(x%*%bparm)	

		pw<-numeric(ndose)
		pwx<-matrix(numeric(ndose*nbase),nrow=ndose)
		for(i in 1:ndose){
			phold<-predout[i]+covc
			if(binary){
				phold<-plogis(phold)
				ph2<-phold*(1-phold)
				pw[i]<-mean(ph2)
				pwx[i,]<-apply(ph2*x,2,mean)
			}else{
				pw[i]<-1
				pwx[i,]<-xmean
			}
			predout[i]<-mean(phold)
		}
		phold<-predR+covc
		if(binary){
			phold<-plogis(phold)
			ph2<-phold*(1-phold)
			pw0<-mean(ph2)
			pwx0<-apply(ph2*x,2,mean)
		}else{
			pw0<-1
			pwx0<-xmean
		}
		predR<-mean(phold)
	}else{
		if(binary){
			predout<-plogis(predout)
			pw<-predout*(1-predout)
			predR<-plogis(predR)
			pw0<-predR*(1-predR)
		}else{
			pw<-rep(1,ndose)
			pw0<-1			
		}
	}

	fitdif <- predout - predR
	
	if(modType==4){
		dlam<-doselev^parm[2]
		edlam<-ed50^parm[2]
		drat<-dlam/(dlam+edlam)
		g<-parm[3]*drat

		dere0<-   rep(1,ndose)
		deremax<- drat
		dered50<- -parm[2]*g*edlam/(dlam+edlam) 
		derlam<-  ifelse(doselev==0,0,
			 g*(log(doselev)- (dlam*log(doselev)+
			 edlam*log(ed50))/(dlam+edlam)))

		L<-pw*cbind(dered50,derlam,deremax,dere0)

		### repeat for reference dose
		dlamR<-dref^parm[2]
		edlamR<-ed50^parm[2]
		dratR<-dlamR/(dlamR+edlamR)
		gR<-parm[3]*dratR

		dere0R<-   1
		deremaxR<- dratR
		dered50R<- -parm[2]*gR*edlamR/(dlamR+edlamR) 
		derlamR<-  ifelse(dref==0,0,
			 gR*(log(dref)- (dlamR*log(dref)+
			 edlamR*log(ed50))/(dlamR+edlamR)))

		LR<-pw0*c(dered50R,derlamR,deremaxR,dere0R)
	}else{
		dere0<-   rep(1,ndose)
		deremax<- doselev/(ed50+doselev)
		dered50<- -parm[2]*deremax*ed50/(ed50+doselev) 

		L<-pw*cbind(dered50,deremax,dere0)

		### repeat for reference dose
		dere0R<-   1
		deremaxR<- dref/(ed50+dref)
		dered50R<- -parm[2]*deremaxR*ed50/(ed50+dref) 

		LR<-pw0*c(dered50R,deremaxR,dere0R)
	}
	
	if(nbase) L<-cbind(L,pwx)
	sepred<- sqrt(diag(L%*%vfcov%*%t(L)))
	
	if(nbase) LR<-c(LR,pwx0)
	seref<-sqrt(t(LR)%*%vfcov%*%LR)
	
	Ldif<-t(t(L)-LR)
	vdif<-diag((Ldif)%*%vfcov%*%t(Ldif))
	vdif[vdif<0]<-0                     ### rounding with near colinearity
	sedif<- sqrt(vdif)
	
	lb<- predout-qnorm(clev)*sepred
	ub<- predout+qnorm(clev)*sepred
	
	lbdif<- fitdif-qnorm(clev)*sedif
	ubdif<- fitdif+qnorm(clev)*sedif
	
	names(predout)<-doselev
	names(sepred)<-doselev
	names(lb)<-doselev
	names(ub)<-doselev
	names(fitdif)<-doselev
	names(sedif)<-doselev
	names(lbdif)<-doselev
	names(ubdif)<-doselev
		
		return(list(pred=predout,lb=lb,ub=ub,se=sepred,
								fitdif=fitdif,lbdif=lbdif,
								ubdif=ubdif,sedif=sedif))
}


