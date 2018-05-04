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
coef.fitEmaxB<-function(object, ...){
	prot<-length(unique(object$prot))
	if(object$pboAdj)prot<-0
	modType<-object$modType
	return(as.matrix(object$estanfit)[,1:(modType+(prot-1))])
}

####
coef.emaxsimB<-function(object, ...){
	return(object$est)
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
	prot<-length(unique(object$prot))
	if(object$pboAdj)prot<-0
	modType<-object$modType
	if(object$binary){
		warning("sigma returns NA for binary endpoint")
		return(NA)
	}else{
		return(as.vector(as.matrix(object$estanfit)[,modType+prot]))
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

