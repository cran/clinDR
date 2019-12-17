coef.fitEmax<-function(object, ...){
	return(object$fit$estimate)
}

####
coef.emaxsim<-function(object, ...){
	fitType<-object$fitType
	est<-object$est4
	colnames(est)<-NULL
	t3<- fitType=='3'
	est[t3,1:3]<-object$est3[t3,]
	ta<- (!t3) & (!(fitType=='4'))
	est[ta,1:2]<-object$estA[ta,]
	
	return(list(fitType=fitType, est=est))
}

####
coef.fitEmaxB<-function(object, local=FALSE, ...){
	binary<-object$binary
	matm<-as.matrix(object$estanfit)	
	localParm<-object$localParm
	if(is.null(localParm))localParm<-FALSE
	if(local & !localParm)stop('local=TRUE but no local parameters fit')
	
	if(local){
		nc<-1
		matm<-matm[,'difTarget',drop=FALSE]
	}else{
		if(localParm){ 
			if(binary)nc<-ncol(matm)-3 else nc<-ncol(matm)-4
		}else  if(binary)nc<-ncol(matm)-1 else nc<-ncol(matm)-2
	}
	
	return(matm[,1:nc,drop=FALSE])
}

####
coef.emaxsimB<-function(object, local=FALSE, ...){
	localParm<-object$localParm
	if(is.null(localParm))localParm<-FALSE
	if(local & !localParm)stop('local=TRUE but no local parameters fit')
	
	est<-object$est
	nc<-ncol(est)
	if((localParm) && local){
		est<-est[,nc,drop=FALSE]
	}else if(localParm)est<-est[,1:(nc-1)]
	
	return(est)
}

####
sigma.fitEmax<-function(object, ...){
	if(object$binary){
		warning("sigma returns NA for binary endpoint")
		return(NA)
	}else return(object$residSD)
}

####
sigma.emaxsim<-function(object, ...){
	if(object$binary){
		warning("sigma returns NA for binary endpoint")
		return(NA)
	}else return(object$residSD)
}

####

sigma.fitEmaxB<-function(object, ...){
	if(object$binary){
		warning("sigma returns NA for binary endpoint")
		return(NA)
	}else{
		matm<-as.matrix(object$estanfit)	
		return(as.vector(matm[,'sigma[1]']))
	}
}
#####
sigma.emaxsimB<-function(object, ...){
	if(object$binary){
		warning("sigma returns NA for binary endpoint")
		return(NA)
	}else return(object$residSD)
}

####
vcov.fitEmax<-function(object, ...){
	return(object$fit$vc)
}

####
vcov.emaxsim<-function(object, ...){
	return(object$vc)
}

