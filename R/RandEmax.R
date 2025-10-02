"RandEmax"<-
function(n,doselev,parmEmax,parmE0,
			p50,parmED50=c(3,0.79,0.6),parmLambda=c(3.03,18.15,0,6),
			resSD,dfSD=Inf,binary=FALSE)
{

		.Deprecated("RandEmax",package="clinDR",
		msg=paste0('This function has been replaced by randomEmax ',
							'that inputs a prior object created by function ',
							'emaxPrior.control.  RandEmax will be removed in a  ',
							'function version of package clinDR'))

	if(length(parmE0)!=2)stop('parme0 input must be length 2')
	if(length(parmEmax)!=2)stop('parmemax input must be length 2')
	
	dord<-order(doselev)
	doselev<-doselev[dord]
	dose<-rep(doselev,n[dord])

	parm<-NULL
	if(binary){
		resSD<-NULL
	}else{
		if(missing(resSD))stop('resSD must be specified for normal data')
		if(!is.numeric(dfSD)){'dfSD is invalid'}
	}

	###genP will be passed to genFun
	###genP must contain n,resSD,doselev,dose + any other inputs 
	###used by genfun

	genP<-list(n=n,resSD=resSD,dfSD=dfSD,doselev=doselev,dose=dose,
				binary=binary,
				parmE0=parmE0,p50=p50,parmED50=parmED50,parmEmax=parmEmax,
				parmLambda=parmLambda)

	genFun <- function(obj){
		###
		###         genFun must return population mean level for each dose
		###                group, parameters defining the DR curve
		###                (application specific), and the response
		###                data, in the order listed
		
		binary<-obj$binary
		### e0
		e0<-rnorm(1,obj$parmE0[1],obj$parmE0[2])

		### emax
		emax<-rnorm(1,obj$parmEmax[1],obj$parmEmax[2])

		### ed50
		parmED50<-obj$parmED50
		led50<-log(p50)+parmED50[2]+parmED50[3]*rt(1,parmED50[1])

		### lambda
		parmL<-obj$parmLambda
		if(length(parmL)==1 & isTRUE(all.equal(parmL,1))){
			p3<-TRUE
		}else{
			p3<-FALSE
		 lambda<-parmL[3] + 
				(parmL[4]-parmL[3])*rbeta(1,parmL[1],parmL[2])
		}
		
		if(!binary){
			dfSD<-obj$dfSD
			resSD<-obj$resSD
			if(is.finite(dfSD))resSD<-resSD*sqrt(dfSD/rchisq(1,dfSD))
		}

		if(p3){parm<-c(led50,emax,e0)
		}else parm<-c(led50,lambda,emax,e0)
		meanlev<-emaxfun(obj$doselev,parm)

		if(binary){
			meanlev<-plogis(meanlev)
			meanR<-rep(meanlev,obj$n)
			y<-rbinom(sum(obj$n),1,meanR)
		}else{
			meanR<-rep(meanlev,obj$n)
			y<-rnorm(sum(obj$n), mean = meanR, sd = resSD)
		}

		return(list(meanlev=meanlev,parm=parm,resSD=resSD,y=y))
	}

	return( list(genP=genP,genFun=genFun)  )
}
