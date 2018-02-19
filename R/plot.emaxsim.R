"plot.emaxsim" <-
function(x,id=x$idmax,plotDif= TRUE,...)
{
	### produces qq plot for z-statistics for primary comparison of
	### idmax dose with placebo

    if(plotDif){
       fitdifv<-x$fitpredv-x$fitpredv[,1]
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
    altord<-x$fitType[ord]
    modChar<-as.character(x$modType)
    plot(xx,y,
         type="n",
         main=paste("QQ plot of population-Z for dose=",
	    doselev[id],sep=""),
                    xlab="Theoretical Quantiles",ylab="Sample Quantiles",... )
    points(xx[altord!=modChar],y[altord!=modChar],col="red")
    points(xx[altord==modChar],y[altord==modChar])
    if(plotDif)mtext('Difference with placebo')
	abline(0, 1)
        cat("Results based on alternative model fits are displayed in red\n")
	return(invisible())
}

