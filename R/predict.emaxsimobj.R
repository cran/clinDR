"predict.emaxsimobj" <-
function(object,dose, dref=0,  ...){

	binary<-object$binary

   predE<-predictEmax(dose,object$estA,object$est3,object$est4,
                      object$vc,object$fitType,dref,
                      maxE=max(object$dose),binary=binary)

   return(predE)
}



"predictEmax" <-
function(dose,estA,est3,est4,vc,fitType,dref=0,maxE,binary){

	Ndose<-length(dose)

   if(fitType=='L'){  # linear
       parm<-estA
       fitpred<-parm[1]+parm[2]*dose  
       sepred<-sqrt(vc[1,1]+vc[2,2]*dose^2 + 2*vc[1,2]*dose)
       fitref<-parm[1]+parm[2]*dref  
       seref<-sqrt(vc[1,1]+vc[2,2]*dref^2 + 2*vc[1,2]*dref)
	   if(!binary){
		   fitdif<-parm[2]*(dose-dref)
		   sedif<-abs(dose-dref)*sqrt(vc[2,2])
	   }else{
			dreg<-cbind(rep(1,Ndose),dose)		
			d0<-c(1,dref)
	   }
   }else if(fitType=='LL'){
   	parm<-estA
   	fitpred<-parm[1]+parm[2]*log(dose+1.0)
   	sepred<-sqrt(vc[1,1]+vc[2,2]*(log(dose+1.0))^2 + 
   							 	2*vc[1,2]*log(dose+1.0))
   	fitref<-parm[1]+parm[2]*log(dref+1.0)
   	seref<-sqrt(vc[1,1]+vc[2,2]*(log(dref+1.0))^2 + 
   								2*vc[1,2]*log(dref+1.0))
   	if(!binary){
   		fitdif<-parm[2]*(log(dose+1.0)-log(dref+1.0))
   		sedif<-abs(log(dose+1.0)-log(dref+1.0))*sqrt(vc[2,2])
   	}else{
   		dreg<-cbind(rep(1,Ndose),log(dose+1.0))		
   		d0<-c(1,log(dref+1.0))
   	}
   }else if(fitType=='E'){
   	parm<-estA
   	fitpred<-parm[1]+parm[2]*exp(dose/maxE)
   	sepred<-sqrt(vc[1,1]+vc[2,2]*exp(2*dose/maxE) + 
   							 	2*vc[1,2]*exp(dose/maxE))
   	fitref<-parm[1]+parm[2]*exp(dref/maxE)
   	seref<-sqrt(vc[1,1]+vc[2,2]*exp(2*dref/maxE) + 
                    2*vc[1,2]*exp(dref/maxE))
	   if(!binary){
		   fitdif<-parm[2]*(exp(dose/maxE)-exp(dref/maxE))
		   sedif<-abs(exp(dose/maxE)-exp(dref/maxE))*sqrt(vc[2,2])
	   }else{
			dreg<-cbind(rep(1,Ndose),exp(dose/maxE))		
			d0<-c(1,exp(dref/maxE))
		}
   }else{
       nparm<-as.numeric(fitType)
       if(nparm==3){est<-est3}else est<-est4
       fit<-list(estimate=est,vc=vc)
	   ### create fitEmax object to re-use predict code
	   ### inputs not used by predict are set to null
	   fit<-list(fit=fit,y=NULL,dose=NULL,modType=nparm,prot=1,
				 count=NULL,binary=binary,pboAdj=FALSE,residSD=NULL,
				 gofTest=NULL,nllmod=NULL,nbase=0)
	   class(fit)<-'fitEmax'
       seout<-predict(fit,dose,dref=dref,binary=binary)
       fitpred<-seout$pred
       fitdif<-seout$fitdif
       sepred<-seout$se
       sedif<-seout$sedif
   }
	if(binary && fitType%in%c('L','LL','E')){
		fitpred<-plogis(fitpred)
		fitref<-plogis(fitref)
		fitdif<-fitpred-fitref
		sdmult<-fitpred*(1-fitpred)
		sepred<- sdmult*sepred
		refmult<-fitref*(1-fitref)
		seref<-refmult*seref
		lcov<-dreg%*%vc%*%d0
		sedif<-sepred^2+seref^2-2*sdmult*refmult[1]*lcov
		sedif[sedif<0]<-0
		sedif<-sqrt(sedif)
	}

   
   names(fitpred)<-dose
   names(fitdif)<-dose
   names(sepred)<-dose
   names(sedif)<-dose
   return(list(fitpred=fitpred,fitdif=fitdif,sepred=sepred,sedif=sedif))
}

"predict.emaxsimBobj" <-
function(object,dose, dref=0, clev=0.9,  ...){
	
	predout<-predict(object$bfit, dosevec=dose, dref=dref, clev=clev, ... )
	
	return(list(pred=predout$pred,lb=predout$lb,ub=predout$ub,se=predout$se,
		 fitdif=predout$fitdif,lbdif=predout$lbdif,
		 ubdif=predout$ubdif,sedif=predout$sedif))

}
