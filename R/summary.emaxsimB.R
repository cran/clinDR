"summary.emaxsimB" <-
function(object,testalpha=0.05,clev=c('0.9','0.95','0.8'),seSim= FALSE,...)
{
	clev<-match.arg(clev)
	jsel<-match(clev,c('0.95','0.9','0.8'))
	clev<-as.numeric(clev)
	
	binary<-object$binary
  modType<-object$modType
  doselev<-object$genObj$genP$doselev
  Ndose<-length(doselev)
  n<-object$genObj$genP$n
  negEmax<-object$negEmax
  dirEff<- 1-2*(negEmax)
  ub<-object$ub
  lb<-object$lb
  fitdifP<-object$predpop[,2:Ndose]-object$predpop[,1]
  fitdifv<-object$fitdifv[,2:Ndose]
  sedifv<-object$sedifv[,2:Ndose]
  mv<-object$mv
  sdv<-object$sdv
  pVal<-object$pVal
  idmax<-object$idmax
  nsim<-length(pVal)
  
  ## select best pairwise function (return index to dif)
  if(!negEmax){f<-function(x){xmax<-max(x); which(x==xmax)[1]}
  }else f<-function(x){xmin<-min(x); which(x==xmin)[1]}
  
  cat(paste("\n",object$description,"\n\n"))
	cat(paste("Number of simulations: ",nsim,'\n',sep=''))

	### pairwise mean+sv	ndose<-length(doselev)
	mdifv <- mv[,2:Ndose] - mv[, 1]
	if(!binary){
		vp<-scale((n[1] - 1) * sdv[, 1]^2 + 
			  scale(sdv[,2:Ndose]^2,center=FALSE,scale=1/(n[2:Ndose]-1)),center=FALSE,
				scale=(n[1]+n[2:Ndose]-2))
		semdifv<-sqrt(scale(vp,center=FALSE,scale=1/(1/n[1] + 1/n[2:Ndose])))
	}else{
		### add 1/2 y/n if 0/1 rate for se
		tol<-sqrt(.Machine$double.eps)
		tmpm<-mv
		for(i in 1:Ndose){
			tmpm[mv[,i]<tol,i]<- 0.5/(n[i]+1)
			tmpm[mv[,i]>1-tol,i]<-(n[i]+0.5)/(n[i]+1)
		}
		semdifv<-sqrt(tmpm[,1]*(1-tmpm[,1])/n[1] + 
								scale(tmpm[,2:Ndose]*(1-tmpm[,2:Ndose]),
											center=FALSE,scale=n[2:Ndose]))
	}

	### power
	cat(paste("\nPower for 1-sided tests at level ", testalpha," :",sep=""))
	pow<- mean(pVal<=testalpha)
	mean.pow <- mean(dirEff*mdifv[,idmax-1]/semdifv[,idmax-1] >qnorm(1-testalpha))
	if(seSim== TRUE){
	cat(paste("\n  Global null test based on MCP-mod:                          ",round(pow,3),
        "(",round(sqrt(mean(pow)*(1-mean(pow))/nsim),4),")",
        "\n",sep=""))
	cat(paste("  Pairwise comparison (simple, no model) of high dose vs pbo: ",round(mean.pow,3),
        "(",round(sqrt(mean(mean.pow)*(1-mean(mean.pow))/nsim),4),")",
        "\n",sep=""))
	}else{
	cat(paste("\n  Global null test based on MCP-mod:                           ",round(pow,3),
        "\n",sep=""))
	cat(paste("  Pairwise comparison (simple, no model) of high dose vs pbo:  ",round(mean.pow,3),
        "\n",sep=""))
	}
	
	### confidence interval coverage and bias
	actlev <- apply(t(lb[,,jsel])<=fitdifP & t(ub[,,jsel])>=fitdifP,2,mean)
	
	mean.actlev<- apply(abs(mdifv - fitdifP)/semdifv <=qnorm(clev+(1-clev)/2),2,mean)
	seldif<-apply(mdifv,1,f)
	seldif<-cbind(seq_along(seldif),seldif)
	bestm<-mdifv[seldif]
	bestpop<-fitdifP[seldif]
	sel.actlev <- mean(abs(bestm - bestpop)/semdifv[seldif]<=
										 	qnorm(clev+(1-clev)/2))
	sel.up.err<- mean((bestm - bestpop)/semdifv[seldif]>=
											qnorm(clev+(1-clev)/2))
	sel.low.err<- mean((bestm - bestpop)/semdifv[seldif]<=
										 	-qnorm(clev+(1-clev)/2))	
	se.actlev<-sqrt( actlev*(1-actlev)/nsim )
	se.mean.actlev<-sqrt( mean.actlev*(1-mean.actlev)/nsim )
	
	cat(paste("\n\nCoverage probabilities for nominal ",
        clev," intervals [Dose-PBO]:",sep=""))
	names(actlev)<-doselev[2:Ndose]
	names(mean.actlev)<-doselev[2:Ndose]
	names(se.actlev)<-doselev[2:Ndose]
	names(se.mean.actlev)<-doselev[2:Ndose]
	cat(paste("\nBayesian Dose response modeling posterior intervals:\n"))
	print(round(actlev,3))
	if(seSim== TRUE){
		cat(paste("Simulation standard errors:\n"))
		print(round(se.actlev,4))
	}
	cat(paste("\nPairwise comparisons:\n"))
	print(round(mean.actlev,3))
	if(seSim== TRUE){
		cat(paste("Simulation standard errors:\n"))
		print(round(se.mean.actlev,4))
	}
	cat(paste("Most favorable pairwise comparison:\n"))	
	cat(paste(round(sel.actlev,3),' Intervals too low(',
						round(sel.low.err,3),
						') Intervals too high(',round(sel.up.err,3)),')\n',sep='')
	
	bias<-apply(fitdifv-fitdifP,2,mean)
	se.bias<-sqrt( apply(fitdifv,2,var)/nsim )
	names(bias)<-doselev[2:Ndose]
	names(se.bias)<-doselev[2:Ndose]
	cat(paste("\n\nBias from Bayesian dose response modeling [DOSE-PBO, EST-POP, EST=posterior median]:\n"))
	print(round(bias,2))
	if(seSim== TRUE){
		cat(paste("Simulation standard errors:\n"))
		print(round(se.bias,3))
	}
	sel.mbias<-mean(bestm-bestpop)
	cat(paste("Bias in the most favorable pairwise comparison:\n"))	
	cat(paste(round(sel.mbias,2),'\n'))
	
	sd.sedifv <- apply(sedifv, 2, mean)
	names(sd.sedifv)<-doselev[2:Ndose]
	sd.semdifv <- apply(semdifv, 2, mean)
	names(sd.semdifv )<-doselev[2:Ndose]

	### mean squared errors
	sd.sedifv <- apply(sedifv, 2, mean)
	names(sd.sedifv)<-doselev[2:Ndose]
	sd.semdifv <- apply(semdifv, 2, mean)
	names(sd.semdifv )<-doselev[2:Ndose]
	
  mse.sedifv <- sqrt( apply((fitdifv-fitdifP)^2,2,mean) )
	names(mse.sedifv)<-doselev[2:Ndose]
	mse.pair<-sqrt( apply((mdifv-fitdifP)^2,2,mean)  )
	names(mse.pair)<-doselev[2:Ndose]
	sel.mse.pair<-sqrt( mean((bestm-bestpop)^2) )

	cat(paste("\n\nSquare Root Mean Squared Error [Dose-PBO]:\n"))
	cat(paste("Bayesian dose response modeling (EST=posterior median) :","\n"))
	print(round(mse.sedifv,3))
	cat(paste("\n","Pairwise comparisons:","\n"))
	print(round(mse.pair,3))
	cat(paste("\nMost favorable pairwise comparison:","\n"))
	cat(paste(round(sel.mse.pair,3),"\n"))

	if(seSim== TRUE){
		cat(paste("\nNote:  (Standard errors) are simulation errors based ",
         "on ",nsim," simulations","\n",sep=""))
		cat(paste("Note:  Distinguish (Standard errors) from ",
     	          "Bayesian posterior SD" ,"\n",sep=""))
	}
	
	mdiverge<-mean(object$divergence)
	cat(paste("\nThe mean of the proportion of diverent MCMC iterations: ",
						round(mdiverge,4),"\n"))	
	
	return(invisible(list(powMCPMOD=pow,PowMean=mean.pow,
												covMod=actlev,
												covMean=mean.actlev,covSelMean=sel.actlev,
												errLowSelMean=sel.low.err,errUpSelMean=sel.up.err,
												biasMod=bias,biasSelMean=sel.mbias,
												summarySEmod=sd.sedifv,summarySEmean=sd.semdifv,
												mseMod=mse.sedifv,mseMean=mse.pair,
												mseSelMean=sel.mse.pair,divergence=mdiverge)))
}

