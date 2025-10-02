"summary.emaxsimobj" <-
function(object,...){


    noFit<-( (as.character(object$modType)!= object$fitType) &
             (!object$negC) & (!object$bigC) )

	cat(paste("Original model converged:       ",!noFit,"\n",sep=""))
	cat(paste("ED50 estimate>ED50 upper limit: ",object$bigC,"\n",sep=""))
	cat(paste("ED50 estimate<ED50 lower limit: ",object$negC,"\n",sep=""))
	cat(paste("Type of  model fit:             ",object$fitType,"\n",sep=""))
	print(summary(object$fit))
	return(invisible())
}


"summary.emaxsimBobj" <-
function(object,...){
	
	dat<-cbind(dose=object$bfit$dose,y=object$bfit$y)
	print(apply(dat,2,summary))

	cat("\n")
	
	print(object$bfit$estanfit)
	return(invisible())
	
}

