"plot.emaxsimobj" <-
  function(x,xlim,xat=NULL,ylim,xlab,ylab,plotDif=FALSE,
           plotResid=FALSE,clev=0.9,plotPop=c('m','3','4'),negC= FALSE,logScale=FALSE, 
  				 predict=TRUE,plot=TRUE,...)
  {
    binary<-x$binary
    doselev<-sort(unique(x$dose))
    est4<-x$est4
    est3<-x$est3
    
    cimul<-abs(qnorm((1-clev)/2))
    n<-table(x$dose)
    
    plotPop<-match.arg(plotPop)
   	if((plotPop=='3')&&(length(x$pop)!=3))stop("Invalid populaton Emax parameters")
   	if((plotPop=='4')&&(length(x$pop)!=4))stop("Invalid populaton Emax parameters")   
    ngrid<-101
    
    if(!plotDif){
      dm<-x$dm
      fitpred<-x$fitpred
      predpop<-x$predpop
      sepred<-x$sepred
      resid<-dm-fitpred
      
      if(!binary){
        sem<-x$residSD/sqrt(n)
        seP<-sqrt(sepred^2+sem^2)
        low<-fitpred-cimul*sepred
        high<-fitpred+cimul*sepred
        lowP<-fitpred-cimul*seP
        highP<-fitpred+cimul*seP
      }else{
        semL<-1/sqrt(fitpred*(1-fitpred)*n)   ##logit scale
        sepredL<-sepred/(fitpred*(1-fitpred))
        sePL<-sqrt(sepredL^2+semL^2)
        low<-plogis(qlogis(fitpred)-cimul*sepredL)
        high<-plogis(qlogis(fitpred)+cimul*sepredL)
        lowP<-plogis(qlogis(fitpred)-cimul*sePL)
        highP<-plogis(qlogis(fitpred)+cimul*sePL)
      }
      
    }else{
      if(!binary){
        sem<-sqrt((x$residSD/sqrt(n))^2 + (x$residSD/sqrt(n[1]))^2 )
      }else{
        sem<-sqrt(x$fitpred*(1-x$fitpred)/n + x$fitpred[1]*(1-x$fitpred[1])/n[1])
      }
      seP<-sqrt(x$sedif^2+sem^2)
      seP[1]<-0.0
      
      dm<-x$dm-x$dm[1]
      fitpred<-x$fitpred-x$fitpred[1]
      predpop<-x$predpop-x$predpop[1]
      sepred<-x$sedif
      resid<-dm-fitpred
      
      low<-fitpred-cimul*sepred
      high<-fitpred+cimul*sepred
      lowP<-fitpred-cimul*seP
      highP<-fitpred+cimul*seP
    }
    
    ### set plot limits and labels
    if(missing(xlim))xlim<-c(0-.05*max(doselev),1.05*max(doselev))
    if(!is.null(xat)){
      if(min(xat)<xlim[1] | max(xat)>xlim[2]) stop('Tickmark locations are outside X limit')
    }    
    if(missing(xlab))xlab<-"Dose"
    if(missing(ylim)){  
      if(!plotResid){
        ylim<-c(min(dm,lowP,predpop), 
                max(dm,highP,predpop))
      }else if(plotResid)ylim<-range(resid)
      ylim<-c(ylim[1]-.05*(ylim[2]-ylim[1]),ylim[2]+.05*(ylim[2]-ylim[1]))
    }
    if(plotDif && missing(ylab)){
      ylab<-'Diff with PBO'
    }else if(plotResid && missing(ylab)){   
      ylab<-'Residuals'
    }else if(missing(ylab))ylab<-"Y" 
    
    
    if(plotResid & !logScale){
      gp2<-ggplot(data.frame(doselev,resid),aes(x=doselev,y=resid))+
        geom_point(shape=8,color='red',size=4)
      gp2<-gp2+xlab(xlab)+ylab(ylab)
      gp2<-gp2+coord_cartesian(xlim=xlim,ylim=ylim)
      gp2<-gp2+geom_hline(yintercept=0,linetype=2)
      #return(invisible())
    }else if(plotResid & logScale){
      x0 <- doselev       
      if(sum(x0==0)){
        xtemp <- doselev[2]^2/doselev[3]
        doselevlog <- doselev
        doselevlog[doselev==0] <- xtemp 
        xlimlog <- c(log(xtemp)-(log(xlim[2])-log(xtemp))/10, 
                     log(xlim[2])+(log(xlim[2])-log(xtemp))/10)        
      }else{
        xtemp <- doselev[1]
        doselevlog <- doselev
        xlimlog <- c(log(xtemp)-(log(xlim[2])-log(xtemp))/10, 
                     log(xlim[2])+(log(xlim[2])-log(xtemp))/10)
      }
      gp2<-ggplot(data.frame(doselevlog,resid),aes(x=log(doselevlog),y=resid))+
        geom_point(shape=8,color='red',size=4)
      gp2<-gp2+xlab(xlab)+ylab(ylab)
      gp2<-gp2+coord_cartesian(xlim=xlimlog,ylim=ylim)
      gp2<-gp2+geom_hline(yintercept=0,linetype=2)  
      if(is.null(xat))  gp2 <- gp2 + scale_x_continuous(breaks=log(doselevlog),
                                                            labels=doselev)       
    }
    
    werrbar<-min(diff(sort(unique(doselev))))*(0.4)
    
    if(!plotResid & !logScale){
      
      df1<-data.frame(low,high,as.numeric(lowP),as.numeric(highP),doselev,dm)
      names(df1)<-c('low','high','lowP','highP','doselev','dm')
      gp2<-ggplot(df1,aes(x=doselev,y=dm))
      gp2<-gp2+geom_point(shape=8,color='red',size=4)
      if(predict)gp2<-gp2+geom_errorbar(aes(ymin=lowP,ymax=highP),width=werrbar,size=1.1,color='grey')
      gp2<-gp2+geom_errorbar(aes(ymin=low,ymax=high),width=0,size=1.1,color='black')
      gp2<-gp2+xlab(xlab)+ylab(ylab)
      gp2<-gp2+coord_cartesian(xlim=xlim,ylim=ylim)
      
      ### compute predictions on grid
      dgrid<-seq(doselev[1],doselev[length(doselev)],length=ngrid)
      if(!plotDif){
        fitgrid<-predict(x,dgrid)$fitpred
      }else fitgrid<-predict(x,dgrid)$fitdif
      
      df2<-data.frame(dgrid,fitgrid)
      names(df2)<-c('dgrid','fitgrid')
      gp2<-gp2+geom_line(data=df2,aes(x=dgrid,y=fitgrid),
                         col='black')
      
      
      if(plotPop=='m' || is.null(x$pop)){
        dfg<-data.frame(doselev,predpop)
        names(dfg)<-c('doselev','predpop')
        gp2<-gp2+geom_line(data=dfg,aes(x=doselev,y=predpop),
                           col='black',linetype=2)
      }else {
        pm<-emaxfun(dgrid,x$pop)
        pm0<-emaxfun(0,x$pop)
        if(binary){
          pm<-plogis(pm)
          pm0<-plogis(pm0)
        }
        if(plotDif)pm<-pm-pm0
        dfg<-data.frame(dgrid,pm)
        names(dfg)<-c('dgrid','pm')
        gp2<-gp2+geom_line(data=dfg,aes(x=dgrid,y=pm),
                           col='black',linetype=2)
      }
      
      if((x$modType==3 || x$modType==4) && !plotDif){
        if(x$modType==3)est<-est3 else est<-est4
        if(isTRUE(negC) && isTRUE(x$negC)){
          pm<-emaxfun(dgrid,est)
          if(binary)pm<-plogis(pm)
          dfg<-data.frame(dgrid,pm)
          names(dfg)<-c('dgrid','pm')
          gp2<-gp2+geom_line(data=dfg,aes(x=dgrid,y=pm),
                             col='black')
        }
      }
    }else if(!plotResid & logScale){
      
      x0 <- doselev 
      
      if(sum(x0==0)){
        ##emaxsim function requires placebo
        ##no need to include without placebo case  
        
        xtemp <- doselev[2]^2/doselev[3]
        doselevlog <- doselev
        doselevlog[doselev==0] <- xtemp
        xlimlog <- c(log(xtemp)-(log(xlim[2])-log(xtemp))/10, 
                     log(xlim[2])+(log(xlim[2])-log(xtemp))/10)
        werrbarlog <- min(diff(sort(unique(log(doselevlog)))))*(0.4)
        ##do not include dose=0 for plotdiff        
        if(plotDif){
        	mindose <- min(doselev[doselev>0])
        	index <- which(doselevlog>=mindose)
        	doselev <- doselev[index]
        	doselevlog <- doselevlog[index]
        	predpop <- predpop[index]
        	low <- low[index]
        	high <- high[index]
        	lowP <- lowP[index]
        	highP <- highP[index]
        	dm <- dm[index]
        	xlimlog <- c(log(mindose)-(log(xlim[2])-log(mindose))/10, 
        							 log(xlim[2])+(log(xlim[2])-log(mindose))/10)  
        	xat <- xat[index]
        }
        
        
        df1<-data.frame(low,high,as.numeric(lowP),as.numeric(highP),doselevlog,dm)
        names(df1)<-c('low','high','lowP','highP','doselevlog','dm')
        gp2<-ggplot(data=df1,aes(x=log(doselevlog),y=dm))
        gp2<-gp2+geom_point(shape=8,color='red',size=4)        
        if(predict)gp2<-gp2+geom_errorbar(aes(ymin=lowP,ymax=highP),width=werrbarlog,size=1.1,color='grey')
        gp2<-gp2+geom_errorbar(aes(ymin=low,ymax=high),width=0,size=1.1,color='black')
        gp2<-gp2+xlab(xlab)+ylab(ylab)      
        gp2<-gp2+coord_cartesian(xlim=xlimlog,ylim=ylim) 
        
        ### compute predictions on grid
        
        dgrid<-seq(doselev[1],doselev[length(doselev)],length=ngrid)
        dgrid <- c(dgrid, doselev)
        dgrid <- sort(unique(dgrid))
        
        if(!plotDif){
          fitgrid<-predict(x,dgrid)$fitpred
        }else fitgrid<-predict(x,dgrid)$fitdif
        
        dgridlog <- dgrid 
        dgridlog[dgrid==0] <- xtemp
        
        ##do not include dose=0 for plotdiff        
        if(plotDif){
        	mindose <- min(doselev[doselev>0])
        	index <- which(doselevlog>=mindose)
          indexgrid <- which(dgridlog>=mindose)   
          dgridlog <- dgridlog[indexgrid]
          fitgrid <- fitgrid[indexgrid]
          df2<-data.frame(dgridlog,fitgrid)
          names(df2)<-c('dgridlog','fitgrid')
          gp2<-gp2+geom_line(data=df2,aes(x=log(dgridlog),y=fitgrid),
          									 col='black')            
        }else if(!plotDif){
        	df2<-data.frame(dgridlog,fitgrid)
        	names(df2)<-c('dgridlog','fitgrid')
        	df21 <- subset(df2, df2$dgridlog <= doselevlog[2] & df2$dgridlog >= doselevlog[1])
        	#df22 <- subset(df2, df2$dgridlog > doselevlog[2])
        	df22 <- subset(df2, df2$dgridlog >= doselevlog[2])
        	gp2<-gp2+geom_line(data=df21,aes(x=log(dgridlog),y=fitgrid),
        										 col='blue')
        	gp2<-gp2+geom_line(data=df22,aes(x=log(dgridlog),y=fitgrid),
        										 col='black')   
        	
        }
        
        if(plotPop=='m' || is.null(x$pop)){
          dfg<-data.frame(doselevlog,predpop)
          names(dfg)<-c('doselevlog','predpop')
          if(plotDif){
          	gp2<-gp2+geom_line(data=dfg,aes(x=log(doselevlog),y=predpop),
          										 col='black',linetype=2)
          }else if (!plotDif){
          	dfg1 <- subset(dfg, dfg$doselevlog <= doselevlog[2] & dfg$doselevlog >= doselevlog[1])
          	dfg2 <- subset(dfg, dfg$doselevlog >= doselevlog[2])
          	gp2<-gp2+geom_line(data=dfg1,aes(x=log(doselevlog),y=predpop),
          										 col='blue',linetype=2)          
          	gp2<-gp2+geom_line(data=dfg2,aes(x=log(doselevlog),y=predpop),
          										 col='black',linetype=2)
          }

        }else if(plotPop%in%c('3','4')){
          if((plotPop=='3')&&(length(x$pop)!=3))stop("Invalid populaton Emax parameters")
          if((plotPop=='4')&&(length(x$pop)!=4))stop("Invalid populaton Emax parameters")
          pm<-emaxfun(dgrid,x$pop)
          pm0<-emaxfun(0,x$pop)
          if(binary){
            pm<-plogis(pm)
            pm0<-plogis(pm0)
          }
          ##do not include dose=0 for plotdiff           
          if(plotDif){
          	pm<-pm-pm0
          	mindose <- min(doselev[doselev>0])          	
          	indexgrid <- which(dgridlog>=mindose)          	
          	pm <- pm[indexgrid]
          	dgridveclog <- dgridveclog[indexgrid]            	
          }
          dfg<-data.frame(dgridlog,pm)
          names(dfg)<-c('dgridlog','pm')
          dfg1 <- subset(dfg, dfg$dgridlog <= doselevlog[2] & dfg$dgridlog >= doselevlog[1])
          dfg2 <- subset(dfg, dfg$dgridlog >= doselevlog[2])
          gp2<-gp2+geom_line(data=dfg1,aes(x=log(dgridlog),y=pm),
                             col='blue')          
          gp2<-gp2+geom_line(data=dfg2,aes(x=log(dgridlog),y=pm),
                             col='black', linetype=2)          
        }
        
        if((x$modType==3 || x$modType==4) && !plotDif){
          if(x$modType==3)est<-est3 else est<-est4
          if(isTRUE(negC) && isTRUE(x$negC)){
            pm<-emaxfun(dgrid,est)
            if(binary)pm<-plogis(pm)
            dfg<-data.frame(dgridlog,pm)
            names(dfg)<-c('dgridlog','pm')
            dfg1 <- subset(dfg, dfg$dgridlog <= doselevlog[2] & dfg$dgridlog >= doselevlog[1])
            dfg2 <- subset(dfg, dfg$dgridlog >= doselevlog[2])
            gp2<-gp2+geom_line(data=dfg1,aes(x=log(dgridlog),y=pm),
                               col='blue')          
            gp2<-gp2+geom_line(data=dfg2,aes(x=log(dgridlog),y=pm),
                               col='black', linetype=2)              
          }
        }
        if(is.null(xat)) gp2<-gp2 + scale_x_continuous(breaks=log(doselevlog),
                                      labels=doselev)       
      }
      
    }
    
    
    ## remove the vertical grid lines
    gp2 <- gp2+ theme(panel.grid.major.x=element_blank(),
                      panel.grid.minor.x=element_blank(),
                      panel.grid.major.y=element_line(size=0.1))   
    
    if(!is.null(xat)){
      if(!logScale)      gp2 <- gp2 + scale_x_continuous(breaks=xat, 
                                                        labels =xat)
      if(logScale){
        xatbench <- xat
        xat[xat==0] <- doselev[2]^2/doselev[3]
        gp2 <- gp2 + scale_x_continuous(breaks=log(xat), 
                                            labels =xatbench)           
      }         
    }    
    
    gp2<-gp2+ theme_bw()
    if(plot){
      print(gp2)
      if(!plotResid) cat("Note:  Dashed curve is population, solid curve is estimated\n")
      if(!plotResid & !plotDif & logScale) cat("Note:  Blue curve is the approximate curve in log scale\n")
    }
    return(invisible(gp2))
  }

"plot.emaxsimBobj" <-
function(x, clev=0.9, plotDif=FALSE, plotPop=c('m','3','4'), 
				 logScale=FALSE, plotResid=FALSE, plot=TRUE, ... )
{
	if(!missing(logScale) || !missing(plotResid))stop('logScale and plotResid options not currently implemented')
	out<-plot(x$bfit, clev=clev, plotDif=plotDif, plot=FALSE, ... )  	
	lplot<-out$lplot
	
  plotPop<-match.arg(plotPop)
 	if((plotPop=='3')&&(length(x$pop)!=3))stop("Invalid populaton Emax parameters")
 	if((plotPop=='4')&&(length(x$pop)!=4))stop("Invalid populaton Emax parameters")
  ngrid<-101

	pop<-x$pop
	if(!is.null(pop) && plotPop!='m'){
		dgrid<-seq(0,1.1*max(x$bfit$dose),length=ngrid)
		if(x$binary){
			popg<-plogis(emaxfun(dgrid,x$pop))
		}else  popg<-emaxfun(dgrid,x$pop)
	}else{
		dgrid<-sort(unique(x$dose))
		popg<-x$predpop
	}
	
	if(plotDif){
		popg<-popg-popg[1]
		if(logScale){
			popg<-popg[-1]
			dgrid<-dgrid[-1]
		}
	}
	
	lplot<-lplot+geom_line(data=data.frame(dgrid=dgrid,popg=popg),
                            aes(x=dgrid,y=popg),color='black',
												 		linetype='dashed',size=1.1) 	
	
	print(lplot)
  if(!plotResid) cat("Note:  Dashed curve is population, solid curve is estimated\n")
	return(invisible(lplot))
}


