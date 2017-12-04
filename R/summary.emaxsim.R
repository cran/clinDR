"summary.emaxsim" <-
function(object,testalpha=0.05,clev=0.9,seSim= FALSE,...)
{

	binary<-object$binary
  nsim<-length(object$fitType)
  modType<-object$modType
  doselev<-object$genObj$genP$doselev
  n<-object$genObj$genP$n
  negEmax<-object$negEmax
  dirEff<- 1-2*(negEmax)
  fitdifv<-object$fitpredv-object$fitpredv[,1]
  fitdifP<-object$predpop-object$predpop[,1]
  noFit<- (modType==4 & apply(is.na(object$est4),1,any)) | (modType==3 & apply(is.na(object$est3),1,any))
  propFit<-table(object$fitType)/nsim

  cat(paste("\n",object$description,"\n\n"))
	cat(paste("Number of simulations:                     ",nsim,'\n',sep=''))
	if(seSim== TRUE){
        cat(paste("Proportion failing to converge:            ",round(mean(noFit),3),
              "(",round(sqrt(mean(noFit)*(1-mean(noFit))/nsim),4),")",
              "\n",sep=""))
        if(object$switchMod){
            cat(paste("Proportion converged but ED50>upper limit: ",round(mean(object$bigC),3),
                  "(",round(sqrt(mean(object$bigC)*(1-mean(object$bigC))/nsim),4),")",
                  "\n",sep=""))
            cat(paste("Proportion converged but ED50<lower limit: ",round(mean(object$negC),3),
                  "(",round(sqrt(mean(object$negC)*(1-mean(object$negC))/nsim),4),")",
                  "\n",sep=""))
        }
        cat(paste("Proportions with modType:"))
        print(round(propFit,2))
        cat(paste("Simulation standard errors:"))
        print(round(sqrt(propFit*(1-propFit)/nsim),2))
	}else{
        cat(paste("Proportion failing to converge:            ",round(mean(noFit),3),
              "\n",sep=""))
        if(object$switchMod){
            cat(paste("Proportion converged but ED50>upper limit: ",round(mean(object$bigC),3),
                  "\n",sep=""))
            cat(paste("Proportion converged but ED50<lower limit: ",round(mean(object$negC),3),
                  "\n",sep=""))
        }
        cat(paste("Proportions with modType:"))
        print(round(propFit,2))
	}

	### pairwise mean+sv	ndose<-length(doselev)
	mdifv <- object$mv - object$mv[, 1]
	if(!binary){
		vp<-scale((n[1] - 1) * object$sdv[, 1]^2 + 
			  scale(object$sdv^2,center=FALSE,scale=1/(n-1)),center=FALSE,
				scale=(n[1]+n-2))
		semdifv<-sqrt(scale(vp,center=FALSE,scale=1/(1/n[1] + 1/n)))
	}else{
		semdifv<-sqrt(object$mv[,1]*(1-object$mv[,1])/n[1] + 
								scale(object$mv*(1-object$mv),center=FALSE,scale=n))
	}

	### power
	cat(paste("\nPower for 1-sided tests at level ", testalpha," :",sep=""))
	pow<- mean(object$pVal<=testalpha)
	mean.pow <- mean(dirEff*mdifv[,object$idmax]/semdifv[,object$idmax] >qnorm(1-testalpha))
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
	actlev <- apply(abs(fitdifv - fitdifP)/object$sedifv <=qnorm(clev+(1-clev)/2),2,mean)
	mean.actlev<- apply(abs(mdifv - fitdifP)/semdifv <=qnorm(clev+(1-clev)/2),2,mean)
	se.actlev<-sqrt( actlev*(1-actlev)/nsim )
	se.mean.actlev<-sqrt( mean.actlev*(1-mean.actlev)/nsim )
	cat(paste("\n\nCoverage probabilities for nominal ",
        clev," intervals [Dose-PBO]:",sep=""))
	names(actlev)<-doselev
	names(mean.actlev)<-doselev
	names(se.actlev)<-doselev
	names(se.mean.actlev)<-doselev
	cat(paste("\nDose response modeling:\n"))
	print(round(actlev[-1],3))
	if(seSim== TRUE){
		cat(paste("Simulation standard errors:\n"))
		print(round(se.actlev[-1],4))
	}
	cat(paste("Pairwise comparisons:\n"))
	print(round(mean.actlev[-1],3))
	if(seSim== TRUE){
		cat(paste("Simulation standard errors:\n"))
		print(round(se.mean.actlev[-1],4))
	}
	bias<-apply(fitdifv-fitdifP,2,mean)
	se.bias<-sqrt( apply(fitdifv,2,var)/nsim )
	names(bias)<-doselev
	names(se.bias)<-doselev
	cat(paste("\n\nMean bias from dose response modeling [EST-POP]:\n"))
	print(round(bias[-1],2))
	if(seSim== TRUE){
		cat(paste("Simulation standard errors:\n"))
		print(round(se.bias[-1],3))
	}
	bias<-apply(fitdifv-fitdifP,2,median)
	se.bias<-sequant(p=.5, sig2=(.75*apply(fitdifv,2,iqr))^2, n=nsim)
	names(bias)<-doselev
	names(se.bias)<-doselev
	cat(paste("\nMedian bias from dose response modeling [EST-POP]:\n"))
	print(round(bias[-1],2))
	if(seSim== TRUE){
		cat(paste("Simulation standard errors:\n"))
		print(round(se.bias[-1],3))
	}

	### summarize standard errors

	cat(paste("\n\nReported SEs by dose group [Dose-PBO]:\n"))
	cat(paste("Dose response modeling:","\n"))
	sd.sedifv <- apply(object$sedifv, 2, summary)
	colnames(sd.sedifv)<-doselev
	print(sd.sedifv[,-1])

	cat(paste("\nPairwise comparisons:","\n"))
	sd.semdifv <- apply(semdifv, 2, summary)
	colnames(sd.semdifv )<-doselev
	print(sd.semdifv[,-1])

	### mean squared errors
    mse.sedifv <- sqrt( apply((fitdifv-fitdifP)^2,2,mean) )
	names(mse.sedifv)<-doselev
	mse.pair<-sqrt( apply((mdifv-fitdifP)^2,2,mean)  )
	names(mse.pair)<-doselev

	cat(paste("\nSquare Root Mean Squared Error:\n"))
	cat(paste("Dose response modeling:","\n"))
	print(round(mse.sedifv[-1],3))
	cat(paste("\n","Pairwise comparisons:","\n"))
	print(round(mse.pair[-1],3))

	if(seSim== TRUE){
		cat(paste("\nNote:  (Standard errors) are simulation errors based ",
         "on ",nsim," simulations","\n",sep=""))
		cat(paste("Note:  Distinguish (Standard errors) from ",
     	          "reported non-linear estmation SE" ,"\n",sep=""))
	}
	return(invisible())
}
