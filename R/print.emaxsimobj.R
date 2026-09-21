"print.emaxsimobj" <-
function(x,nprint=min(length(x$y),20),...){

	return(pemaxsimobj(x,nprint))
}



"print.emaxsimBobj" <-
function(x,nprint=min(length(x$y),20),...){

	return(pemaxsimobj(x,nprint))
}


"pemaxsimobj" <-
function(x,nprint, ...){
	if(length(nprint)==1)nprint<-c(1,nprint)
	tmp<-cbind(1:length(x$y),x$dose,x$y)[nprint[1]:nprint[2],]
	colnames(tmp)<-c("ID","dose","y")
	rownames(tmp)<-rep("",nrow(tmp))
	cat(paste("IDs printed are ",nprint[1],":",nprint[2]," out of ",
        length(x$y),"\n",sep=""))

	print(tmp,quote= FALSE)

	return(invisible(tmp))
}

