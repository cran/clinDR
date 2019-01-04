startEmax<-function(y,dose,baseline,count=rep(1,length(y)),
            modType=3,binary=FALSE,
            lbED50=doselev[2]/10,ubED50=max(doselev),
            lbLambda=0.5,ubLambda=5){

    if(binary && !(y%in%c(0,1)))stop('y must be 0/1 for binary data')
    if(! modType%in%c(3,4))stop("modType must be 3 or 4")

    ### require at least 4 points because ends of designs
    ### are unstably estimated and discarded
    if((length(unique(dose))<4) & (modType==4))stop(paste("At least 4 dose group means", 
        " are required to obtain starting values"))
		
		nbase<-0
		if(!missing(baseline)){
			if(is.vector(baseline)){
				if(length(baseline)!=length(y))stop('baseline and y differ in length')
				baseline<-matrix(baseline,ncol=1)
			}else if (is.matrix(baseline)){
				nbase<-ncol(baseline)
				if(nrow(baseline)!=length(y))stop('baseline and y differ in length')
			}
		}
		
		parm<-modType+nbase

    ### compute dose group means (with possible baseline adjustment)
    doselev <- sort(unique(dose))
    ndose<-length(doselev)

    dord<-order(dose)
    dose<-dose[dord]
    y<-y[dord]
    count<-count[dord]
    if(!missing(baseline)){
    	baseline<-baseline[dord,,drop=FALSE]
      colnames(baseline)<- paste("baseline",1:nbase,sep="")
    }

    if (missing(baseline)) {
        dm<-numeric(ndose)
        for(i in 1:ndose){ 
            dm[i]<-weighted.mean(y[dose==doselev[i]],
                                 count[dose==doselev[i]],na.rm=TRUE)
        }
    }else {
        if(!binary){
            anova.fit <- lm(y ~ factor(dose) + baseline,weights=count) 
        }else{
            anova.fit<-glm(y ~ factor(dose) + baseline,family=binomial(), 
                           weights=count,glm.control(maxit = 100))
        }
        b<- coef(anova.fit)[(1+ndose):(ndose+nbase)]
        for(i in 1:nbase)names(b)[i]<- paste("b",i,sep='')
        bdat<-matrix(rep(apply(baseline,2,mean),
        														ndose),ncol=nbase,byrow=TRUE)
        predat <- data.frame(dose = doselev, baseline=I(bdat)) 
        pred.anova <- predict(anova.fit, newdata=predat, se.fit = TRUE, type='response')
        dm <- pred.anova$fit
    }

    dord<-order(dm) 
    if(dord[1]<dord[ndose]){
      if(!binary){
          e0<-dm[dord[1]]
          emax<-dm[dord[ndose]]-e0
      }else{
          e0<-qlogis(min(max(dm[dord[1]],0.001),0.999))
          emax<-qlogis(min(max(dm[dord[ndose]],0.001),0.999))-e0
      }
    }else{
      if(!binary){
          e0<-dm[dord[ndose]]
          emax<-dm[dord[1]]-e0
      }else{
          e0<-qlogis(min(max(dm[dord[ndose]],0.001),0.999))
          emax<-qlogis(min(max(dm[dord[1]],0.001),0.999))-e0
      }
    }

    if(!binary){
        pY<-(dm[2:(length(dm)-1)] - e0 + 0.0005 * emax)/(1.001 * emax)
    }else{
        pY<-dm[2:(length(dm)-1)]
        pY<-ifelse(pY>.999,.999,pY)
        pY<-ifelse(pY<.001,.001,pY)
    }
    logitY <- log(pY/(1-pY))

    ### log-linear model for logit of the proportion of the mean effect
    if(modType==3){
       ed50 <- mean(exp(-logitY + log(doselev[2:(length(dm)-1)])))
       if(ed50<lbED50)ed50<-lbED50
       if(ed50>ubED50)ed50<-ubED50
       Sparm <- c(log(ed50), emax, e0)
       names(Sparm)<-c("led50","emax","e0")
    }else{
       x<-log(doselev[2:(length(dm)-1)])
       lmfit<-lm(logitY~x,weights=count[2:(length(dm)-1)])
       lambda<-coef(lmfit)[2]
       ed50<- exp( -coef(lmfit)[1]/lambda  )
       if(ed50<lbED50)ed50<-lbED50
       if(ed50>ubED50)ed50<-ubED50
       if(lambda<lbLambda)lambda<-lbLambda
       if(lambda>ubLambda)lambda<-ubLambda
       Sparm<- c(log(ed50), lambda, emax, e0)
       names(Sparm)<-c("led50","lambda","emax","e0")
    }

    if(!missing(baseline))Sparm<-c(Sparm,b)

    return(Sparm)
 }

