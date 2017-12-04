"predict.fitEmax"<-
function(object,dosevec,clev=0.9,int=1,dref=0, ...){

	if(missing(dosevec))stop('dosevec must be specified')

	modType<-object$modType
	binary<-object$binary
	pboAdj<-object$pboAdj
	parms<-object$fit$estimate
	vc<-object$fit$vc
	nprot<-max(as.numeric(object$prot))

	if(! int%in%(c(1:nprot)))stop('The intercept specification is invalid')

	clev<- 0.5+clev/2
	if(pboAdj)parms<-c(parms,rep(0,nprot))  ## put into standard form
	nparm<-length(parms)

	nsub<-nparm-nprot
	if(! nsub%in%c(2,3) )stop('parm has invalid length')
	parms<-parms[c(1:nsub,nsub+int)]

	### put vc into standard form before applying subsetting code for one int
	if(pboAdj){
		if(modType==3){
			tmpcov<-diag(2+nprot)
			diag(tmpcov)<-0
			tmpcov[1:2,1:2]<-vc
			vc<-tmpcov
		}else{
			tmpcov<-diag(3+nprot)
			diag(tmpcov)<-0
			tmpcov[1:3,1:3]<-vc
			vc<-tmpcov
		}
	}

   ### create vc matrix including the single selected intercept
	vsub<-as.vector(vc[1:nsub,nsub+int])
	vc<-cbind(rbind(vc[1:nsub,1:nsub],vsub),c(vsub,vc[nsub+int,nsub+int]))

	mout<-SeEmax(list(parms,vc),dosevec,modType=nsub+1,dref=dref)

	predout<-mout$fitpred
	sepred<-mout$sepred
	lb<- predout-qnorm(clev)*sepred
	ub<- predout+qnorm(clev)*sepred
	if(!binary){
		fitdif<-mout$fitdif
		sedif<-mout$sedif
		lbdif<- fitdif-qnorm(clev)*sedif
		ubdif<- fitdif+qnorm(clev)*sedif
	}else{
		predout<-plogis(predout)
		lb<-plogis(lb)
		ub<-plogis(ub)
		sepred<-predout*(1-predout)*sepred  ### se for proportions

		predref<-plogis(mout$predref)
		seref<-predref*(1-predref)*mout$seref
		fitdif<-predout-predref
		xhold<-sepred^2+seref^2-
				   2*predout*(1-predout)*predref*(1-predref)*mout$covref
		xhold[xhold<0]<-0  ### <0 possible with dref in doselev(rounding)
		sedif<-sqrt(xhold)
		lbdif<-fitdif-qnorm(clev)*sedif
		ubdif<-fitdif+qnorm(clev)*sedif
	}
	
		names(predout)<-dosevec
		names(sepred)<-dosevec
		names(lb)<-dosevec
		names(ub)<-dosevec
		names(fitdif)<-dosevec
		names(sedif)<-dosevec
		names(lbdif)<-dosevec
		names(ubdif)<-dosevec
	
	return(list(pred=predout,lb=lb,ub=ub,se=sepred,
		 fitdif=fitdif,lbdif=lbdif,
		 ubdif=ubdif,sedif=sedif))
}


