plotBdensity<-function(dgrid,parm,
    modelFun=emaxfun,
    qlevL=c(0.025,0.05,0.10,0.25),
    plotDif= FALSE, 
    logit= FALSE, ...) {

    qlevH<- 1-qlevL

    pbo<-modelFun(0,parm)
    if(logit== TRUE)pbo<-plogis(pbo)

    qL<-NULL
    qH<-NULL
    for (dlev in dgrid){
        dr<-modelFun(dlev,parm)
        if(logit== TRUE)dr<-plogis(dr)
        if(plotDif== TRUE)dr<-dr-pbo
        qL<-rbind(qL,quantile(dr,probs=qlevL))
        qH<-rbind(qH,quantile(dr,probs=qlevH))
    }

    DRDensityPlot(dgrid,qL,qH,qlevL=qlevL, ...)

    return(list(qL=qL, qH=qH))
    }


