DRDensityPlot<-function(x,qL,qH,qlevL=c(0.025,0.05,0.10,0.25),
						xlim,ylim,xlab='x',ylab='y'){

    ### check input consistency
    if(!is.matrix(qL))stop("qL must be a matrix")
    if(!is.matrix(qH))stop("qH must be a matrix") 
    if(nrow(qL)!=length(x))stop("rows of qL must correspond to x values")
    if(nrow(qH)!=length(x))stop("rows of qH must correspond to x values")
    if(length(qlevL)!=ncol(qL))stop("columns of qL must correpond to qlevL")
    if(length(qlevL)!=ncol(qH))stop("columns of qH must correpond to 1-qlevL")
		
		### remove row/col names that create warning messages
		rownames(qL)<-NULL
		colnames(qL)<-NULL
		rownames(qH)<-NULL
		colnames(qH)<-NULL
		
    ng<-length(qlevL)
    gscale<-c(1:ng)*floor(70/ng)
    gscale<-paste("gray",gscale,sep="")
    gscale<-rev(gscale)

    if(missing(ylim)){
		ylim=range(as.vector(cbind(qL,qH)))
		ylim<-c(ylim[1]-0.05*(ylim[2]-ylim[1]),ylim[2]+0.05*(ylim[2]-ylim[1]))
	}
	if(missing(xlim))xlim=range(x)

	dfplot<-NULL
    for(icol in c(1:ng)){
		dfplot<-rbind(dfplot,cbind(x,qL[,icol],qH[,icol],1-2*qlevL[icol]))
    }
	dfplot<-data.frame(dfplot)
	names(dfplot)<-c('x','ql','qh','clev')
	dfplot$clev<-factor(dfplot$clev,levels=sort(unique(dfplot$clev),decreasing=TRUE))
	
	### assign variables to avoid error messages in package check
	ql<-NULL
	qh<-NULL
	clev<-NULL

	ggp<-ggplot(data=dfplot,aes(x=x,ymin=ql,ymax=qh,fill=clev))+geom_ribbon()
	ggp<-ggp+coord_cartesian(xlim=xlim,ylim=ylim) + xlab(xlab) + ylab(ylab)
	ggp<-ggp+theme_bw()+scale_fill_grey(start = .9, end = 0.5)
	ggp<-ggp+labs(fill='Level')

	print(ggp)

    return(invisible(ggp))
}



