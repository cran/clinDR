"summary.emaxsimB" <-
function(object,testalpha=0.05,clev=c('0.95','0.9','0.8'),seSim= FALSE,...)
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
  fitdifv<-object$fitpredv[,2:Ndose]-object$fitpredv[,1]
  sedifv<-object$sedifv[,2:Ndose]
  mv<-object$mv
  sdv<-object$sdv
  pVal<-object$pVal
  idmax<-object$idmax
  nsim<-length(pVal)
  
  cat(paste("\n",object$description,"\n\n"))
	cat(paste("Number of simulations:                     ",nsim,'\n',sep=''))

	### pairwise mean+sv	ndose<-length(doselev)
	mdifv <- mv[,2:Ndose] - mv[, 1]
	if(!binary){
		vp<-scale((n[1] - 1) * sdv[, 1]^2 + 
			  scale(sdv[,2:Ndose]^2,center=FALSE,scale=1/(n[2:Ndose]-1)),center=FALSE,
				scale=(n[1]+n[2:Ndose]-2))
		semdifv<-sqrt(scale(vp,center=FALSE,scale=1/(1/n[1] + 1/n[2:Ndose])))
	}else{
		semdifv<-sqrt(mv[,1]*(1-mv[,1])/n[1] + 
								scale(mv[,2:Ndose]*(1-mv[,2:Ndose]),
											center=FALSE,scale=n[2:Ndose]))
	}
	

	### power
	cat(paste("\nPower for 1-sided tests at level ", testalpha," :",sep=""))
	pow<- mean(pVal<=testalpha)
	mean.pow <- mean(dirEff*mdifv[,idmax-1]/semdifv[,idmax-1] >qnorm(1-testalpha))
	if(seSim== TRUE){
	cat(paste("\nPower for global null test based on MCP-mod:      ",round(pow,3),
        "(",round(sqrt(mean(pow)*(1-mean(pow))/nsim),4),")",
        "\n",sep=""))
	cat(paste("Power for pairwise comparion of high dose vs pbo: ",round(mean.pow,3),
        "(",round(sqrt(mean(mean.pow)*(1-mean(mean.pow))/nsim),4),")",
        "\n",sep=""))
	}else{
	cat(paste("\nPower for global null test based on MCP-mod:       ",round(pow,3),
        "\n",sep=""))
	cat(paste("Power for pairwise comparion of high dose vs pbo:  ",round(mean.pow,3),
        "\n",sep=""))
	}
	
	### confidence interval coverage and bias
	actlev <- apply(t(lb[,,jsel])<=fitdifP & t(ub[,,jsel])>=fitdifP,2,mean)
	
	mean.actlev<- apply(abs(mdifv - fitdifP)/semdifv <=qnorm(clev+(1-clev)/2),2,mean)
	
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
	
	bias<-apply(fitdifv-fitdifP,2,mean)
	se.bias<-sqrt( apply(fitdifv,2,var)/nsim )
	names(bias)<-doselev[2:Ndose]
	names(se.bias)<-doselev[2:Ndose]
	cat(paste("\n\nBias from Bayesian dose response modeling [DOSE-PBO, est=posterior mean]:\n"))
	print(round(bias,2))
	if(seSim== TRUE){
		cat(paste("Simulation standard errors:\n"))
		print(round(se.bias,3))
	}

	### summarize standard errors
	cat(paste("\n\nReported SEs by dose group [Dose-PBO]:\n"))
	cat(paste("(Bayesian dose response modeling (posterior SD):","\n"))
	sd.sedifv <- apply(sedifv, 2, mean)
	names(sd.sedifv)<-doselev[2:Ndose]
	print(sd.sedifv)

	cat(paste("\nPairwise comparisons:","\n"))
	sd.semdifv <- apply(semdifv, 2, mean)
	names(sd.semdifv )<-doselev[2:Ndose]
	print(sd.semdifv)

	### mean squared errors
  mse.sedifv <- sqrt( apply((fitdifv-fitdifP)^2,2,mean) )
	names(mse.sedifv)<-doselev[2:Ndose]
	mse.pair<-sqrt( apply((mdifv-fitdifP)^2,2,mean)  )
	names(mse.pair)<-doselev[2:Ndose]

	cat(paste("\n\nSquare Root Mean Squared Error [Dose-PBO]:\n"))
	cat(paste("Bayesian dose response modeling (est=posterior mean) :","\n"))
	print(round(mse.sedifv,3))
	cat(paste("\n","Pairwise comparisons:","\n"))
	print(round(mse.pair,3))

	if(seSim== TRUE){
		cat(paste("\nNote:  (Standard errors) are simulation errors based ",
         "on ",nsim," simulations","\n",sep=""))
		cat(paste("Note:  Distinguish (Standard errors) from ",
     	          "Bayesian posterior SD" ,"\n",sep=""))
	}
	return(invisible())
}

