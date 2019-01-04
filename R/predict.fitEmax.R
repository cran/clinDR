"predict.fitEmax"<-
function(object,dosevec,clev=0.9,int=1,dref=0, xvec=NULL, ...){

	if(missing(dosevec))stop('dosevec must be specified')

	modType<-object$modType
	binary<-object$binary
	pboAdj<-object$pboAdj
	parmstot<-object$fit$estimate
	vc<-object$fit$vc
	nprot<-max(as.numeric(object$prot))
	prot<-as.numeric(object$prot)
	nbase<-object$nbase
	if(nbase && is.null(xvec))xvec<-object$xbase[prot==int,]

	if(! int%in%(c(1:nprot)))stop('The intercept specification is invalid')

	clev<- 0.5+clev/2
	if(pboAdj)parmstot<-c(parmstot,rep(0,nprot))  ## put into standard form
	
	nparm<-length(parmstot)
	nsub<-nparm-nprot-nbase
	if(! nsub%in%c(2,3) )stop('parm has invalid length')
	parms<-parmstot[c(1:nsub,nsub+int)]
	if(nbase>0)parms<-c(parms,parmstot[c((1+nparm-nbase):nparm)])
	
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
	if(!nbase){
		vc<-vc[c(1:nsub,nsub+int),c(1:nsub,nsub+int)]
	}else  vc<-vc[c(1:nsub,nsub+int,(1+nparm-nbase):nparm),
					c(1:nsub,nsub+int,(1+nparm-nbase):nparm)]
	
	mout<-SeEmax(list(parms,vc),dosevec,modType=nsub+1,dref=dref,
							 nbase=nbase,x=xvec, binary=binary, clev=clev)

	return(mout)	
}


