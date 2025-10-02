"print.emaxsim" <-
function(x,nprint=min(length(x$fitType),20),id=x$idmax,digits=3,...){

    nsim<-length(x$fitType)
    doselev<-x$genObj$genP$doselev
    fitdifv<-x$fitpredv-x$fitpredv[,1]
    fitdifP<-x$predpop-x$predpop[,1]
    negC<-x$negC
    bigC<-x$bigC
    pval<-round(x$pVal,digits)

    est4<-x$est4
    est3<-x$est3
    if(x$modType==4){
        est4[,1]<-exp(est4[,1])
    }
    est3[,1]<-exp(est3[,1])
    

    noFit<- (x$modType==4 & apply(is.na(est4),1,any)) | (x$modType==3 & apply(is.na(est3),1,any))

	if(length(nprint)==1)nprint<-c(1,nprint)
	zval<-round((fitdifv[, id]-fitdifP[,id])/x$sedifv[, id],digits)

    if(x$modType==4){
	   tmp<-data.frame(1:length(x$fitType),x$fitType,noFit,negC,
               bigC,round(est4,3),zval,pval,row.names=1)[nprint[1]:nprint[2],]
       colnames(tmp)<-c("fitType","noFit","negC","bigC","ED50","lambda","Emax","E0","StdBias","P-val")
    }else{
	   tmp<-data.frame(1:length(x$fitType),x$fitType,noFit,negC,
               bigC,round(est3,3),zval,pval,row.names=1)[nprint[1]:nprint[2],]
       colnames(tmp)<-c("fitType","noFit","negC","bigC","ED50","Emax","E0","StdBias","P-val")
    }

	cat(paste(x$description,"\n"))
	cat(paste("IDs printed are ",nprint[1],":",nprint[2]," out of ",nsim,"\n\n",sep=""))
	cat(paste("noFit  : No convergence","\n",sep=""))
	cat(paste("negC   : Converged with ED50<lower limit","\n",sep=""))
	cat(paste("bigC   : Converged with ED50>upper limit","\n",sep=""))
	cat(paste("StdBias: (estimate-population)/SE for dose= ",doselev[id]," vs PBO","\n",sep=""))
	cat(paste("P-val  : MCP-Mod P-value for test of no drug effect","\n\n",sep=""))

	print(tmp,quote= FALSE)

	return(invisible(tmp))
}

