  targetCI<-function(object,target,dgrid,cilev=0.80,high= TRUE){

	if(!inherits(object,'emaxsim'))stop('Input must be of class emaxsim')
	if(length(dgrid)==1){
	 doselev<-object$genObj$genP$doselev
	 dgrid<-seq(0,max(doselev),length=dgrid)
	}
	ngrid<-length(dgrid)
	nsim<-length(object$fitType)

	cnorm<-qnorm(cilev)
	predout<-predict(object,dgrid)
	fitdif<-predout$fitdif
	sedif<-predout$sedif

	### reverse comparison if lower is better
	if(!high){
	 fitdif<- -fitdif
	 target<- -target
	}

	lb<-(fitdif-cnorm*sedif>target)

	targetDose<-apply(lb,1,findFirst)
	targetDose[is.finite(targetDose)]<-dgrid[targetDose[is.finite(targetDose)]]

	return(targetDose)

}

findFirst<-function(lbvec){
	return(suppressWarnings(min(which(lbvec))))
}

