"plot.emaxsimB" <-
function(x,id=x$idmax,plotDif= TRUE,...)
{
	### produces qq plot for z-statistics for primary comparison of
	### idmax dose with placebo
	
	if(plotDif){
		fitdifv<-x$fitdifv
		fitdifP<-x$predpop-x$predpop[,1]
    }else{
       fitdifv<-x$fitpredv
       fitdifP<-x$predpop
    }

    doselev<-x$genObj$genP$doselev
    xx<-qnorm(ppoints(nrow(fitdifv)))
    if(plotDif)y<-(fitdifv[, id] - fitdifP[,id])/x$sedifv[, id]
    if(!plotDif)y<-(fitdifv[, id] - fitdifP[,id])/x$sepredv[, id]
    ord<-order(y)
    y<-y[ord]
    plot(xx,y,
         type="n",
         main=paste("QQ plot of population-Z for dose=",
         					 doselev[id],sep=""),
    		 xlab="Theoretical Quantiles",ylab="Sample Quantiles",... )
    if(plotDif)mtext('Difference with placebo')
    points(xx,y)
    abline(0, 1)
    cat("Population-Z is (posterior mean - pop)/(posterior sd)\n")	
    return(invisible())
}

