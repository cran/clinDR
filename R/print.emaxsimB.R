"print.emaxsimB" <- function(x,nprint=min(nsim,20),id=x$idmax,digits=3,...){

	modType<-x$modType	
	doselev<-x$genObj$genP$doselev
	fitdifv<-x$fitpredv-x$fitpredv[,1]
	fitdifP<-x$predpop-x$predpop[,1]
	sedifv<-x$sedifv
	pval<-round(x$pVal,digits)
	nsim<-length(pval)
	
	est<-x$est
	est[,1]<-exp(est[,1])
	if(modType==4)estname<-c("ED50","lambda","Emax","E0")
	else estname<-c("ED50","Emax","E0")


	if(length(nprint)==1)nprint<-c(1,nprint)
	zval<-round((fitdifv[, id]-fitdifP[,id])/sedifv[, id],digits)
	
	tmp<-data.frame(1:nsim,
									round(est,3),zval,pval,row.names=1)[nprint[1]:nprint[2],]
	colnames(tmp)<-c(estname,"StdBias","P-val")
	
	cat(paste(x$description,"\n"))
	cat(paste("IDs printed are ",nprint[1],":",nprint[2]," out of ",nsim,"\n\n",sep=""))

	cat("Parameter estimates are posterior medians\n")
	cat(paste("StdBias: (posterior mean-population)/(posterior SD) for dose= ",doselev[id]," vs PBO","\n",sep=""))
	cat(paste("P-val  : MCP-Mod P-value for test of no drug effect","\n\n",sep=""))
	print(tmp,quote= FALSE) 
	
	
	return(invisible())
}

