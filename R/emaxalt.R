emaxalt<-function(y,dose,modType=3,binary=FALSE,iparm=NA,
					ed50cutoff=2.5*max(doselev),
					ed50lowcutoff=doselev[2]/1000,
					switchMod= TRUE,truncLambda=6){

	if(! modType%in%c(3,4))stop("modType must be 3 or 4")


	### sort by dose to ensure standardizaton of output
	dord<-order(dose)
	dose<-dose[dord]
	y<-   y[dord]
	### summary statistics by dose group
	dm  <- tapply(y, dose, mean)
	if(!binary){
		dsd <- sqrt(tapply(y, dose, var))
	}else dsd<-NULL

	### utility variables
	intercept<-rep(1,length(dose))
	doselev<-sort(unique(dose))
	Ndose<-length(doselev)

	### initialize variables that may not have assigned values
	est4<-rep(NA,4)
	est3<-rep(NA,3)
	estA<-rep(NA,2)

	### assign initial parameter values
	if(!any(is.na(iparm))){
		if(is.null(names(iparm))){
		   if(modType==3){names(iparm)<-c("led50","emax","e0")
		   }else{names(iparm)<-c("led50","lambda","emax","e0")}
		}else if(!all(names(iparm)%in%c("led50","lambda","emax","e0"))){
		   stop("Invalid names for new parm")
		}
		Sparm<-iparm
	}else{
		### 
		### find starting values for the Emax model
		### per SSfp1 and Jim Rogers
		###
		Sparm<-startEmax(y,dose,modType=modType,binary=binary)

		if(modType==4 && Sparm[2]<=0){
				warning(paste("A negative lambda starting",
				   " value was generated and replaced by 1"))
				Sparm[2]<-1
		}
	}

	Sparmcopy<-Sparm  ### save for return in case changed later

	bigC<-F
	negC<-F
	AcceptConv<-F

	if(modType==4){
		fit<-fitEmax(y,dose,iparm=Sparm,modType=4,binary=binary,diagnostics=FALSE,
					 optObj=FALSE)
		
		### check for acceptable convergence and compute se's if ok
		if(!is.null(fit) && (fit$fit$estimate[2]<=truncLambda)){
			est4<-fit$fit$estimate
			ed50<-exp(est4[1])
		
			if(ed50<ed50lowcutoff){negC<-TRUE
			}else if(ed50>ed50cutoff)bigC<-TRUE

			if(!(negC&switchMod) & !(bigC&switchMod)){
				vc<-as.vector(vcov(fit))
				seout<-predict(fit,doselev)
				fitpred <-seout$pred
				sepred <- seout$se
				sedif <-  seout$sedif
				residSD<-fit$residSD
				AcceptConv<-TRUE
				fitType<-"4"
			}else est4<-rep(NA,4)
		}
	}


	### if 4 parm fit not found, select new 3 parm starting values
	if(modType==4 & AcceptConv== FALSE){ 
		Sparm<-startEmax(y,dose,modType=3,binary=binary)
		fit<-fitEmax(y,dose,iparm=Sparm,modType=3,binary=binary,diagnostics=FALSE,
					 optObj=FALSE)
	}else if(modType==3)fit<-fitEmax(y,dose,iparm=Sparm,modType=3,binary=binary,
											diagnostics=FALSE,optObj=FALSE)

	if(!AcceptConv){   
		### check for acceptable convergence and compute se's if ok
		if(!is.null(fit)){
			### check that info matrix is invertible
			est3<-fit$fit$estimate
			ed50<-exp(est3[1])

			negC3<-FALSE
			bigC3<-FALSE
			if(ed50<ed50lowcutoff){negC3<-T
			}else if(ed50>ed50cutoff)bigC3<-T
			if(modType==3){
				negC<-negC3
				bigC<-bigC3
			}

			if(!(negC3&&switchMod) && !(bigC3&&switchMod)){
				vc<-as.vector(vcov(fit))
				seout<-predict(fit,doselev)
				fitpred <-seout$pred
				sepred <- seout$se
				sedif <-  seout$sedif
				residSD<-fit$residSD
				AcceptConv<-TRUE
				fitType<-"3"
			}else est3<-rep(NA,3)
		}
	}

  ### alternative model linear in the parameters if unacceptable Emax fit
	if(!AcceptConv){
		### for binary sig is the deviance (-2ll)

		### linear model
		dfun<-dose
		if(binary){
			fitL <- glm(y ~ dfun,family=binomial())
			sigL<-fitL$deviance 
		}else{
			fitL <- lm(y ~ dfun)
			sigL<-summary(fitL)$sigma
		}
		### log-linear model
		if(binary){
			dfun<- log(dose+1.0)
			fitLL<- glm(y ~ dfun,family=binomial())
			sigLL<-fitLL$deviance
		}else{
			dfun<- log(dose+1.0)
			fitLL<- lm(y ~ dfun)
			sigLL<-summary(fitLL)$sigma
		}
		### (scaled) exp model
		if(binary){
			dfun<-exp(dose/max(doselev))
			fitE<-glm(y~dfun,family=binomial())
			sigE<-fitE$deviance
		}else{
			dfun<-exp(dose/max(doselev))
			fitE<-lm(y~dfun)
			sigE<-summary(fitE)$sigma
		}

		### select best fit
		if(sigL<=sigLL & sigL<=sigE){
			fitType<-"L"
			dfun<-dose
			dfunlev<-doselev
			fit<-fitL
			if(!binary)residSD<-sigL else residSD<-NULL
		}else if(sigLL<=sigE){
			fitType<-"LL"
			dfun<-log(dose+1.0)
			dfunlev<-log(doselev+1.0)
			fit<-fitLL
			if(!binary)residSD<-sigLL else residSD<-NULL
		}else{
			fitType<-"E"
			dfun<-exp(dose/max(doselev))
			dfunlev<-exp(doselev/max(doselev))
			fit<-fitE
			if(!binary)residSD<-sigE else residSD<-NULL
		}

		vc<-vcov(fit)
		predobj<- predict(fit, data.frame(dfun = dfunlev), se.fit = TRUE)
		if(binary){
			fitpred<- plogis(predobj$fit)        ## convert to probabilities
			sdmult<-fitpred*(1-fitpred)
			sepred<- sdmult*predobj$se.fit
			dreg<-cbind(rep(1,Ndose),dfunlev)		
			d0<-c(1,dfunlev[1])
			lcov<-dreg%*%vc%*%d0
			sedif<-sepred^2+(sepred[1])^2-2*sdmult*sdmult[1]*lcov
			sedif[sedif<0]<-0
			sedif<-sqrt(sedif)
		}else{
			fitpred<- predobj$fit
			sepred<- predobj$se.fit
			sdd<-sqrt(vc[2, 2])
			sedif <- dfunlev *sdd 
			if(fitType=="E"){sedif<-sedif- sdd}  
		}
		estA<-coef(fit)
		vc<-as.vector(vc)
	}

	names(fitpred)<-doselev
	names(sepred)<-doselev
	names(sedif)<-doselev
	return(structure(list( dm=dm,dsd=dsd,Sparm=Sparmcopy,
		binary=binary,fitType=fitType,residSD=residSD,
		vc=vc,fitpred=fitpred,sepred=sepred,sedif=sedif,bigC=bigC,
		negC=negC,est4=est4,est3=est3,estA=estA),class='emaxalt') )

}

