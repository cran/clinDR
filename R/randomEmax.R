"randomEmax"<-
function(x,n,doselev,modType=c('4','3'))
{

	modType<-match.arg(modType,c('4','3'))

	dord<-order(doselev)
	doselev<-doselev[dord]
	dose<-rep(doselev,n[dord])

	parm<-NULL
	resSD<-NULL

	###genP will be passed to genFun
	###genP must contain n,resSD,doselev,dose + any other inputs 
	###used by genfun

	genP<-list(n=n,doselev=doselev,dose=dose,x=x,modType=modType)

	genFun <- function(obj){
		###
		###         genFun must return population mean level for each dose
		###                group, parameters defining the DR curve
		###                (application specific), and the response
		###                data, in the order listed
		
		
		modType<-obj$modType
		x<-obj$x
		if(length(x$effDF)>1){
			epDF<-x$effDF[1]
			diftDF<-x$effDF[2]
		}else{
			epDF<-x$effDF
			diftDF<-x$diftDF
		}

		### e0
		e0<-x$epmu+x$epsca*rt(1,epDF)

		### difTarget 
		difTarget<-x$difTargetmu+x$difTargetsca*rt(1,epDF)
		
		if(modType=='3'){
			led50<-x$loged50mu+log(x$p50)+x$loged50sca*rt(1,x$parmDF)	
			lambda<-1
		}else{
			sdiag<-diag(c(x$loged50sca,x$loglamsca))
			scamat<- matrix(c(1,x$parmCor,x$parmCor,1),ncol=2)
			scamat<-sdiag%*%scamat%*%sdiag
			holder<-rmvt(n=1,sigma=scamat,df=x$parmDF)	
			led50<-x$loged50mu+log(x$p50)+holder[1]
			lambda<-exp(x$loglammu+holder[2])
		}

		if(!x$binary){
			resSD<-runif(1,x$sigmalow,x$sigmaup)
		}
		
		emax<-solveEmax(difTarget,x$dTarget,led50,lambda,e0)

		if(modType=='3'){parm<-c(led50,emax,e0)
		}else parm<-c(led50,lambda,emax,e0)
		meanlev<-emaxfun(obj$doselev,parm)

		if(x$binary){
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
