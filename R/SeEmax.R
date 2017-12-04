"SeEmax" <-
function(fit,doselev,modType,dref=0){

	if(! modType %in% c(3,4))stop("modType must be 3 or 4")

		if(class(fit)=='nls'){
			parm<-coef(fit)
			vfcov <- vcov(fit)[1:modType,1:modType]
		}else{
			parm<-fit[[1]][1:modType]
			vfcov<-fit[[2]][1:modType,1:modType]
		}

		ed50<-exp(parm[1])

		if(modType==4){
			fitpred <- emaxfun(doselev, parm)
			fitR <- emaxfun(dref, parm)
			fitdif <- fitpred - fitR

			dlam<-doselev^parm[2]
			edlam<-ed50^parm[2]
			drat<-dlam/(dlam+edlam)
			g<-parm[3]*drat

			dere0<-   rep(1,length(doselev))
			deremax<- drat
			dered50<- -parm[2]*g*edlam/(dlam+edlam) 
			derlam<-  ifelse(doselev==0,0,
				 g*(log(doselev)- (dlam*log(doselev)+
				 edlam*log(ed50))/(dlam+edlam)))

			L<-cbind(dered50,derlam,deremax,dere0)
			sepred<- sqrt(diag(L%*%vfcov%*%t(L)))


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

			LR<-c(dered50R,derlamR,deremaxR,dere0R)
		}else{
			fitpred <- emaxfun(doselev, parm)
			fitR <- emaxfun(dref, parm)
			fitdif <- fitpred - fitR

			dere0<-   rep(1,length(doselev))
			deremax<- doselev/(ed50+doselev)
			dered50<- -parm[2]*deremax*ed50/(ed50+doselev) 

			L<-cbind(dered50,deremax,dere0)
			sepred<- sqrt(diag(L%*%vfcov%*%t(L)))

			### repeat for reference dose
			dere0R<-   1
			deremaxR<- dref/(ed50+dref)
			dered50R<- -parm[2]*deremaxR*ed50/(ed50+dref) 

			LR<-c(dered50R,deremaxR,dere0R)
		}

		Ldif<-t(t(L)-LR)
		vdif<-diag((Ldif)%*%vfcov%*%t(Ldif))
		vdif[vdif<0]<-0                     ### rounding with near colinearity
		sedif<- sqrt(vdif)

		### add variance for reference and covariance
		### between reference and dose estimates for
		### potential use with transformed variables
		### (e.g., logistic regression)

		seref<-sqrt(t(LR)%*%vfcov%*%LR)
		covref<-L%*%vfcov%*%LR
		seref<-as.numeric(seref)
		covref<-as.vector(covref)

	return(list(doselev=doselev,dref=dref,fitpred=fitpred,sepred=sepred,
		fitdif=fitdif,sedif=sedif,
		predref=fitR,seref=seref,covref=covref))
}


