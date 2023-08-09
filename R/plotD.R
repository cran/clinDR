"plotD"<-function (y, dose, baseline, se = TRUE, line = TRUE, 
                   meansOnly=FALSE,sem=NULL,clev = 0.9, 
                   xlab='Dose',ylab='Response', logScale=FALSE) 
{
  cadj <- abs(qnorm((1 - clev)/2))
  if(meansOnly && se && missing(sem))stop('sem must be specified with means only and se=TRUE')
  if(meansOnly && !missing(baseline))stop('Covariate adjustment and means only cannot be specified')
  if(meansOnly){
    dord<-order(dose)
    ym<-y[dord]
    doselev<-dose[dord]
  }else{
    doselev <- sort(unique(dose))
    if(missing(baseline)){
      ym <- tapply(y, dose, mean, na.rm = TRUE)
      sem <- sqrt(tapply(y, dose, var, na.rm = TRUE)/tapply(!is.na(y), 
                                                            dose, sum, na.rm = TRUE))
    }else{
      anova.fit<-lm(y~baseline+factor(dose))
      predat<-data.frame(dose=doselev,baseline=rep(mean(baseline),length(doselev)))
      pred.anova<-predict(anova.fit,predat,se.fit= TRUE)
      ym<-pred.anova$fit
      sem<-pred.anova$se.fit
    }
  }
  
  if(is.null(sem))sem<-numeric(length(ym))
  chvec<-ym + cadj * sem
  clvec<-ym - cadj * sem
  
  if(!logScale){
    ggp<-ggplot(data.frame(ym,doselev,chvec,clvec),aes(x=doselev,y=ym))
    ggp<-ggp + geom_point(size=4,shape=8,color='red')  
    if (se == TRUE) 
      ggp<-ggp + geom_errorbar(aes(ymax=chvec,ymin=clvec))
    if (line == TRUE) 
      ggp<-ggp + geom_line()
    
  }else{
    x0 <- doselev
    if(sum(x0==0)){
      doselevlog <- doselev      
      doselevlog[doselev==0] <- doselev[2]^2/doselev[3]
      ggp<-ggplot(data.frame(ym,doselevlog,chvec,clvec),aes(x=log(doselevlog),y=ym))   
      ggp<-ggp + scale_x_continuous(breaks=log(doselevlog),
                                    labels=c(format(round(doselev, 0), nsmall=0)))       
      ggp<-ggp + geom_point(size=4,shape=8,color='red')  
      if (se == TRUE) 
        ggp<-ggp + geom_errorbar(aes(ymax=chvec,ymin=clvec))
      if (line == TRUE){
        data <- data.frame(ym,doselevlog,chvec,clvec)
        #data1 <- subset(data, data$doselevlog < doselevlog[2] & data$doselevlog >= doselevlog[1])
        data2 <- subset(data, data$doselevlog >= doselevlog[2])
        ggp<-ggp+geom_line(data=data, aes(x=log(doselevlog),ym),color='black',size=0.6, linetype="dashed") 
        ggp<-ggp+geom_line(data=data2, aes(x=log(doselevlog),ym),color='black',size=0.6, linetype="solid") 
        
      } 
      
    }else{
      doselevlog <- doselev
      
      ggp<-ggplot(data.frame(ym,doselevlog,chvec,clvec),aes(x=log(doselevlog),y=ym))   
      ggp<-ggp + scale_x_continuous(breaks=log(doselevlog),
                                    labels=doselev)  
      ggp<-ggp + geom_point(size=4,shape=8,color='red')  
      if (se == TRUE) 
        ggp<-ggp + geom_errorbar(aes(ymax=chvec,ymin=clvec))
      if (line == TRUE) 
        ggp<-ggp + geom_line()
    } 

  
    
  }
  

  ggp <- ggp +ylab(ylab) + xlab(xlab) + ggplot2::theme_bw()   
  ## remove the vertical grid lines
  ggp <- ggp+ ggplot2::theme(panel.grid.major.x=element_blank(),
                    panel.grid.minor.x=element_blank(),
                    panel.grid.major.y=element_line(size=0.1))  
  
  return(list(ggp = ggp, means = ym, se = sem))
}


