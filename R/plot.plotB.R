"plot.plotB"<-function(x,
                     plotDif= FALSE,plotMed= FALSE,plotResid=FALSE,
                     predict= TRUE,logScale=FALSE,
                     xlim,xat=NULL,ylim,xlab,ylab,labac='Act Comp',shapeac=8,colac='red',  
                     symbolLabel='Group',symbolShape=8,symbolColor='red',symbolSize=4,...)
{
  
  activeControl<-x$activeControl
  
  if(plotResid && plotMed)warning(paste("When plotResid=TRUE, ",
                                        "the plotMed option is ignored",sep=""))
  if(activeControl & plotDif)message('Dose response less active control is plotted')
  if(plotResid & plotDif)stop('Both plotResid and plotDif cannot be requested')
  
  dgrid<-x$dgrid
  dose<-x$dose
  clev<-x$clev
  nolegend<-x$nolegend
  doselev <- sort(unique(dose))
  
  cadj <- abs(qnorm((1 - clev)/2))
  low.clev <- (1 - clev)/2
  up.clev <- clev + (1 - clev)/2
  
  ### assign variables to avoid error messages in package check
  yac<-NULL
  aclowLpred<-NULL
  acupLpred<-NULL
  aclowL<-NULL
  acupL<-NULL
  
  dmax<-max(dose,dgrid)
  dmin<-min(dose,dgrid)
  if(missing(xlim))xlim <- c(min(0,dmin-0.1*dmax), 1.1*dmax) 
  if(missing(xlab))xlab <- "Dose" 
  if(!is.null(xat)){
    if(min(xat)<xlim[1] | max(xat)>xlim[2]) stop('Tickmark locations are outside X limit')
  }      
  symbol<-x$symbol
  if(activeControl)levels(symbol)<-c(levels(symbol),labac)
  nsymlev<-length(levels(symbol))
  if(activeControl){
    symbolShape<-c(symbolShape,shapeac)
    symbolColor<-c(symbolColor,colac)
  }
  nshape<-length(symbolShape)
  if(nshape!=nsymlev){
    if(nshape==1){
      symbolShape<-rep(symbolShape,nsymlev)
    }else if(nshape<nsymlev)stop('symbol and symbolShape incompatible')
  }
  ncolor<-length(symbolColor)
  if(ncolor!=nsymlev){
    if(ncolor==1){
      symbolColor<-rep(symbolColor,nsymlev)
    }else if(ncolor<nsymlev)stop('symbol and symbolColor are incompatible')
  }
  
  
  ##########################################################
  ### resid plots
  ##########################################################
  if(plotResid & !logScale){
    ym<-x$pairwise[,'ym']
    pred<-x$modelABS[,'Bmean']
    resid<-ym-pred
    if(missing(ylab))ylab<-'Residuals'
    maxr<-max(resid)
    minr<-min(resid)
    rdif<-maxr-minr
    if(missing(ylim))ylim<-c(minr-.1*rdif,maxr+.1*rdif)
    gp2<-ggplot(data.frame(dose,resid,symbol),aes(x=dose,y=resid))
    gp2<-gp2+geom_hline(yintercept=0,linetype=2)
    gp2<-gp2+geom_point(aes(shape=symbol,color=symbol),size=symbolSize)
    gp2<-gp2+scale_color_manual(name=symbolLabel,values=symbolColor)
    gp2<-gp2+scale_shape_manual(name=symbolLabel,values=symbolShape)   
    gp2<-gp2+ylab(ylab)+xlab(xlab)+theme_bw()
    gp2<-gp2+coord_cartesian(xlim=xlim,ylim=ylim)
    
    if(nolegend)gp2<-gp2 + theme(legend.position = "none")
  
  }else if(plotResid & logScale){
    ym<-x$pairwise[,'ym']
    pred<-x$modelABS[,'Bmean']
    resid<-ym-pred
    if(missing(ylab))ylab<-'Residuals'
    maxr<-max(resid)
    minr<-min(resid)
    rdif<-maxr-minr
    if(missing(ylim))ylim<-c(minr-.1*rdif,maxr+.1*rdif)
    x0 <- dose
    if(sum(x0==0)){
      xtemp <- doselev[2]^2/doselev[3]
      doselog <- dose 
      doselog[dose==0] <- xtemp 
      xlimlog <- c(log(xtemp)-(log(xlim[2])-log(xtemp))/10, 
                   log(xlim[2])+(log(xlim[2])-log(xtemp))/10)        
    }else{
      xtemp <- doselev[1]
      doselog <- dose
      xlimlog <- c(log(xtemp)-(log(xlim[2])-log(xtemp))/10, 
                   log(xlim[2])+(log(xlim[2])-log(xtemp))/10)
    }
    
    gp2<-ggplot(data.frame(doselog,resid,symbol),aes(x=log(doselog),y=resid))

    gp2<-gp2+geom_hline(yintercept=0,linetype=2)
    gp2<-gp2+geom_point(aes(shape=symbol,color=symbol),size=symbolSize)
    gp2<-gp2+scale_color_manual(name=symbolLabel,values=symbolColor)
    gp2<-gp2+scale_shape_manual(name=symbolLabel,values=symbolShape)   
    gp2<-gp2+ylab(ylab)+xlab(xlab)+theme_bw()
    gp2<-gp2+coord_cartesian(xlim=xlimlog,ylim=ylim)
    if(is.null(xat))   gp2<-gp2 + scale_x_continuous(breaks=log(sort(unique(doselog))),
                                  labels=sort(unique(dose)))     
    if(nolegend)gp2<-gp2 + theme(legend.position = "none")

    
  }
 ##########################################################
 ### end resid plotting
 ##########################################################
  
  if(plotDif){
    plotupP<-x$modelDIF[,"upLdifpred"]
    plotlowP<-x$modelDIF[,"lowLdifpred"]
    ylimHP<-max(x$modelDIF[,"upLdifpred"],x$modelDIFG[,"upLdifGpred"])
    ylimLP<-min(x$modelDIF[,"lowLdifpred"],x$modelDIFG[,"lowLdifGpred"])
    
    plotup<-x$modelDIF[,"upLdif"]
    plotlow<-x$modelDIF[,"lowLdif"]
    ylimH<-max(x$modelDIF[,"upLdif"],x$modelDIFG[,"upLdifG"],ylimHP)
    ylimL<-min(x$modelDIF[,"lowLdif"],x$modelDIFG[,"lowLdifG"],ylimLP)
    
    ploty<-x$pairwise[,"ymdif"]
    if(plotMed){
      plotmg<-x$modelDIFG[,"BdifmedG"]
      plotm<-x$modelDIF[,"Bdifmed"]
    }else{
      plotmg<-x$modelDIFG[,"BdifmeanG"]
      plotm<-x$modelDIF[,"Bdifmean"]
    }
  }else{
    plotupP<-x$modelABS[,"upLpred"]
    plotlowP<-x$modelABS[,"lowLpred"]
    ylimHP<-max(x$modelABS[,"upLpred"],x$modelABSG[,"upLGpred"])
    ylimLP<-min(x$modelABS[,"lowLpred"],x$modelABSG[,"lowLGpred"])
    
    plotup<-x$modelABS[,"upL"]
    plotlow<-x$modelABS[,"lowL"]
    ylimH<-max(x$modelABS[,"upL"],x$modelABSG[,"upLG"],ylimHP)
    ylimL<-min(x$modelABS[,"lowL"],x$modelABSG[,"lowLG"],ylimLP)
    
    ploty<-x$pairwise[,"ym"]
    if(plotMed){
      plotmg<-x$modelABSG[,"BmedG"]
    }else{
      plotmg<-x$modelABSG[,"BmeanG"]
    }
  }
  
  if(missing(ylim))ylim<-c(ylimL-0.05*(ylimH-ylimL),ylimH+0.05*(ylimH-ylimL)) 
  if(missing(ylab))ylab <- ifelse(plotDif,"Diff with Placebo","Y") 
  
  if(!plotResid & !logScale){
    werrbar<-min(diff(sort(unique(dose))))*(0.4)
    
    if(activeControl){
      pdose <- length(doselev)
      yac=rep(x$AC['yac'], pdose)
      aclowL=rep(x$AC['aclowL'], pdose)
      acupL=rep(x$AC['acupL'], pdose)
      aclowLpred=rep(x$AC['aclowLpred'], pdose)
      acupLpred=rep(x$AC['acupLpred'], pdose)     
      gp2<-ggplot(data.frame(doselev,plotlow,plotup,plotlowP,plotupP,
                             yac, aclowL,acupL,aclowLpred,acupLpred))
    }else{
      gp2<-ggplot(data.frame(doselev,plotlow,plotup,plotlowP,plotupP))     
    }
    
   
    if(predict)gp2<-gp2+geom_errorbar(aes(x=doselev,ymin=plotlowP,ymax=plotupP),
                                      width=werrbar,size=1.1,color='grey')
    gp2<-gp2+geom_errorbar(aes(x=doselev,ymin=plotlow,ymax=plotup),width=0,size=1.1,color='black')
    gp2<-gp2+geom_line(data=data.frame(dgrid,plotmg),
                       aes(x=dgrid,y=plotmg),col='black',size=1.1)
    gp2<-gp2+ylab(ylab) + xlab(xlab) + theme_bw()
    gp2<-gp2+geom_point(data=data.frame(dose,ploty,symbol),aes(x=dose,y=ploty,shape=symbol,
                                                               color=symbol),size=symbolSize)
    
    
    if(activeControl & !plotDif)xlim[2] <- xlim[2]+max(diff(doselev))
    gp2<-gp2+coord_cartesian(xlim=xlim,ylim=ylim)
    
    if(nolegend)gp2<-gp2 + theme(legend.position = "none")
    
    if(activeControl & !plotDif){
      
      if(predict)gp2<-gp2+geom_errorbar(aes(x=max(doselev)+max(diff(doselev)),
                                            ymin=aclowLpred,ymax=acupLpred),
                                        width=0,size=1.1,color='grey')
      gp2<-gp2+geom_errorbar(aes(x=max(doselev)+max(diff(doselev)),
                                 ymin=aclowL,ymax=acupL),width=0,size=1.1,color='black')    
      
      gp2 <- gp2 + geom_point(aes(x=max(doselev)+max(diff(doselev)),y=yac),
                              size=symbolSize,color=colac,shape=shapeac)      
      
      if(is.null(xat))   gp2 <- gp2 + scale_x_continuous(breaks=c(doselev, max(doselev)+max(diff(doselev))),
                                      labels=c(doselev, labac)) 
      
      ### move legend in DR plot and eliminate plot border
      if(!nolegend)gp2<-gp2 + theme(legend.position = "left")
      gp2<-gp2+theme(axis.line = element_line(color = 'black'))
    }
    gp2<-gp2+scale_color_manual(name=symbolLabel,values=symbolColor,drop=FALSE)
    gp2<-gp2+scale_shape_manual(name=symbolLabel,values=symbolShape,drop=FALSE)
  }
  else if(!plotResid & logScale){
    x0 <- dose
    if(sum(x0==0)){
      
      xtemp <- doselev[2]^2/doselev[3]
      doselog <- dose
      dgridlog <- dgrid
      doselevlog <- doselev 
      doselog[dose==0] <- dgridlog[dgrid==0] <- doselevlog[doselev==0] <- xtemp
      xlimlog <- c(log(xtemp)-(log(xlim[2])-log(xtemp))/10, 
                   log(xlim[2])+(log(xlim[2])-log(xtemp))/10)
      werrbarlog<-min(diff(sort(unique(log(doselog)))))*(0.4)
      data0=data.frame(dgridlog,plotmg)
      data1 <- subset(data0, data0$dgridlog < doselevlog[2] & data0$dgridlog >= xtemp)
      data2 <- subset(data0, data0$dgridlog >= doselevlog[2])    
    }else{
      
      xtemp <- doselev[1]^2/doselev[2]
      doselog <- dose
      dgridlog <- dgrid
      doselevlog <- doselev 
      xlimlog <- c(log(xtemp)-(log(xlim[2])-log(xtemp))/10, 
                   log(xlim[2])+(log(xlim[2])-log(xtemp))/10)
      werrbarlog <- min(diff(sort(unique(log(doselog)))))*(0.4)
      data0=data.frame(dgridlog,plotmg)
      data1 <- subset(data0, data0$dgridlog < doselevlog[1] & data0$dgridlog >= xtemp)
      data2 <- subset(data0, data0$dgridlog >= doselevlog[1])      
    }  
    
    if(activeControl){
      pdose <- length(doselev)
      yac=rep(x$AC['yac'], pdose)
      aclowL=rep(x$AC['aclowL'], pdose)
      acupL=rep(x$AC['acupL'], pdose)
      aclowLpred=rep(x$AC['aclowLpred'], pdose)
      acupLpred=rep(x$AC['acupLpred'], pdose)      
      gp2<-ggplot(data.frame(doselevlog,plotlow,plotup,plotlowP,plotupP, 
                             yac, aclowL,acupL,aclowLpred,acupLpred))
    }else{
      gp2<-ggplot(data.frame(doselevlog,plotlow,plotup,plotlowP,plotupP))
    }
   
    if(predict)gp2<-gp2+geom_errorbar(aes(x=log(doselevlog),ymin=plotlowP,ymax=plotupP),
                                      width=werrbarlog,size=1.1,color='grey')
    gp2<-gp2+geom_errorbar(aes(x=log(doselevlog),ymin=plotlow,ymax=plotup),width=0,size=1.1,color='black')
    gp2<-gp2+geom_line(data=data0,
                       aes(x=log(dgridlog),y=plotmg),col='black',size=1.1, linetype="dashed")
    gp2<-gp2+geom_line(data=data2,
                       aes(x=log(dgridlog),y=plotmg),col='black',size=1.1, linetype="solid")    
    
    gp2<-gp2+ylab(ylab) + xlab(xlab) + theme_bw()
    gp2<-gp2+geom_point(data=data.frame(doselog,ploty,symbol),aes(x=log(doselog),y=ploty,shape=symbol,
                                                                  color=symbol),size=symbolSize)
    
    if(activeControl & !plotDif)xlimlog[2] <- xlimlog[2]+max(diff(log(doselevlog)))
     gp2<-gp2+coord_cartesian(xlim=xlimlog,ylim=ylim)
    
    if(nolegend)gp2<-gp2 + theme(legend.position = "none")
    if(!activeControl & is.null(xat)) gp2 <- gp2 + scale_x_continuous(breaks=log(doselevlog),
                                                         labels=doselev)    
    if(activeControl & !plotDif){

      if(predict)gp2<-gp2+geom_errorbar(aes(x=log(max(doselevlog))+max(diff(log(doselevlog))),
                                            ymin=aclowLpred,ymax=acupLpred),
                                        width=0,size=1.1,color='grey')
      gp2<-gp2+geom_errorbar(aes(x=log(max(doselevlog))+max(diff(log(doselevlog))),
                                 ymin=aclowL,ymax=acupL),width=0,size=1.1,color='black')    
      if(is.null(xat))   gp2 <- gp2 + scale_x_continuous(breaks=c(log(doselevlog), log(max(doselevlog))+max(diff(log(doselevlog)))),
                                      labels=c(doselev, labac)) 
      gp2 <- gp2 + geom_point(aes(x=log(max(doselevlog))+max(diff(log(doselevlog))),y=yac),
                              size=symbolSize,color=colac,shape=shapeac)      
      
      ### move legend in DR plot and eliminate plot border
      if(!nolegend)gp2<-gp2 + theme(legend.position = "left")
      gp2<-gp2+theme(axis.line = element_line(color = 'black'))
    }    
    gp2<-gp2+scale_color_manual(name=symbolLabel,values=symbolColor,drop=FALSE)
    gp2<-gp2+scale_shape_manual(name=symbolLabel,values=symbolShape,drop=FALSE)    
  }
  
  ## remove the vertical grid lines
  gp2 <- gp2+ theme(panel.grid.major.x=element_blank(),
                    panel.grid.minor.x=element_blank(),
                    panel.grid.major.y=element_line(size=0.1))  
  
  if(!is.null(xat)){
    if(activeControl){
    if(!logScale)      gp2 <- gp2 + scale_x_continuous(breaks=c(xat, max(doselev)+max(diff(doselev))),
                                                  labels=c(xat, labac)) 
    if(logScale){
      xatbench <- xat
      xat[xat==0] <- doselev[2]^2/doselev[3]
      gp2 <- gp2 + scale_x_continuous(breaks=c(log(xat), log(max(doselevlog))+max(diff(log(doselevlog)))),
                                                    labels=c(xatbench, labac))           
    }  
    }else{
      if(!logScale)      gp2 <- gp2 + scale_x_continuous(breaks=xat,
                                                    labels=xat)
      if(logScale){
        xatbench <- xat
        xat[xat==0] <- doselev[2]^2/doselev[3]
        gp2 <- gp2 + scale_x_continuous(breaks=c(log(xat)),
                                        labels=xatbench)           
      } 
    }
  }    
  
  return(gp2)
}


