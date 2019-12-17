solveEmax<-function(target,dose,led50,lambda,e0,pboadj=TRUE){
	if(pboadj)padj<-0 else padj<-e0
	dl<-dose^lambda
	return((dl+exp(led50*lambda))*(target-padj)/dl)
}

solveDose<-function(target,led50,lambda,emax,e0,pboadj=TRUE){
	if(pboadj)tar<-target else tar<-target-e0
	if(sign(tar)!=sign(emax))stop('The Emax parameter and targetted pbo-adjusted effect must have the same sign')
	if(abs(tar)>abs(emax))stop('The pbo-adjusted target exceeds the Emax parameter')	
	return(exp(led50)*(tar/(emax-tar))^(1/lambda))
}
