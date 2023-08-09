"plot.fitEmax" <-
  function(x,int=0,plotResid=FALSE,clev=0.9,
           predict=TRUE,plotci=TRUE,plotDif=FALSE,
           xlab='Dose',
  				 ylab=ifelse(plotResid,'Residuals',ifelse(plotDif,
                 'Difference With Placebo','Response')),
  				 ncol=NULL,
           symbol=NULL,symbolLabel='Group',symbolShape=8,symbolColor='red',symbolSize=4,
           bwidth=NULL, xlim=NULL, xat=NULL, ylim=NULL, logScale=FALSE, ngrid=200,
           plot=TRUE, ...){
  	
	outplot<-combplot(x=x,bayes=FALSE,int=int,plotResid=plotResid,
					clev=clev, predict=predict,plotci=plotci,plotDif=plotDif,
          xlab=xlab,ylab=ylab,
					ncol=ncol,
          symbol=symbol,symbolLabel=symbolLabel,symbolShape=symbolShape,
					symbolColor=symbolColor,symbolSize=symbolSize,
          bwidth=bwidth, xlim=xlim, xat=xat, ylim=ylim, logScale=logScale, ngrid=ngrid,
          plot=plot, ...)
	
	return(invisible(outplot))  	 	
}
  	
"plot.fitEmaxB" <-
  function(x,int=0,plotResid=FALSE,clev=0.9,
           predict=TRUE,plotci=TRUE,plotDif=FALSE,
           xlab='Dose',
   				 ylab=ifelse(plotResid,'Residuals',ifelse(plotDif,
                 'Difference With Placebo','Response')), 				 
  				 ncol=NULL,
           symbol=NULL,symbolLabel='Group',symbolShape=8,symbolColor='red',symbolSize=4,
           bwidth=NULL, xlim=NULL, xat=NULL, ylim=NULL, logScale=FALSE, ngrid=200,
           plot=TRUE, ...){
  	
	outplot<-combplot(x=x,bayes=TRUE,int=int,plotResid=plotResid,
					clev=clev, predict=predict,plotci=plotci,plotDif=plotDif,
          xlab=xlab,ylab=ylab,
					ncol=ncol,
          symbol=symbol,symbolLabel=symbolLabel,symbolShape=symbolShape,
					symbolColor=symbolColor,symbolSize=symbolSize,
          bwidth=bwidth, 
					xlim=xlim,xat=xat,ylim=ylim,
          logScale=logScale, ngrid=ngrid,
          plot=plot, ...)
	return(invisible(outplot))  	 	
}
 
  	
 combplot <-  	
  function(x,bayes,int=0,plotResid=FALSE,clev=0.9,
	       predict=TRUE,plotci=TRUE,plotDif=FALSE,
	       xlab='Dose',
 				 ylab=ifelse(plotResid,'Residuals',ifelse(plotDif,
							'Difference With Placebo','Response')), 				  			 
 				 ncol=NULL,
	       symbol=NULL,symbolLabel='Group',symbolShape=8,symbolColor='red',
 				 symbolSize=4,
	       bwidth=NULL,
				 xlim=NULL, xat=NULL,
				 ylim=NULL,
	       logScale=FALSE,
				 ngrid=ngrid,
	       plot=TRUE, ...){
	
    #### plot of a continuous emax model (in x) and dose group means
    ####
  	if(plotDif && is.null(ylab)){
  		ylab<-'Y Diff'
  	}else if(plotResid && is.null(ylab)){   
  		ylab<-'Residuals'
  	}else if(is.null(ylab))ylab<-"Y" 
	 	clevup<- 0.5+clev/2
		clevlow<-0.5-clev/2   
		cadj<- abs(qnorm((1-clev)/2))	
    nolegend<-is.null(symbol)
    
    y<-x$y
    dose<-x$dose
    prot<-as.numeric(x$prot)
    protlab<-levels(x$prot)
    count<-x$count
    binary<-x$binary
    
    if(bayes){
    	dimFit<-x$dimFit
    	vcest<-x$vcest
    	if(binary && dimFit)y<-plogis(y)
    }
    
    if(!bayes){
    	if(!binary){
    		residSD<-sigma(x)
    	}else residSD<-NA
  	}
    
    if(isTRUE(plotDif & plotResid))warning('plotDif is ignored when plotResid is TRUE')
    
    if(is.null(symbol)){symbol<-factor(rep(1,length(y)))
    }else symbol<-factor(symbol)  ### make sure symbol is in factor format
    
    if(length(symbol)!=length(y))stop('Symbol variable has invalid length')
    
    if(length(symbol)!=length(y))stop('Symbol variable has invalid length')
    
    
    
    if(int>0){
      nprot<-int
      nstart<-int
      doselev<-sort(unique(dose[prot==int]))
    }else{
      nprot<-max(prot)
      nstart<-1
      doselev<-sort(unique(dose))
    }
    ### global dose range
    dmax<-max(doselev)*1.1
    dmin<-min(c(0,doselev))-0.1*dmax
    if(is.null(xlim))xlim<-c(dmin,dmax)
    dgrid<-seq(0,xlim[2],length=ngrid)
    dgrid <- c(dgrid, doselev)
    dgrid <- sort(unique(dgrid)) 
    ngridsup<-length(dgrid)
    
    
    
    
    if(!is.null(xat)){
      if(min(xat)<xlim[1] | max(xat)>xlim[2]) stop('Tickmark locations are outside X limit')
    }
    
    ### dose level results
    protD<-NULL
    doseLevVec<-NULL
    fitvec<-NULL
    cilvec<-NULL
    cihvec<-NULL
    pilvec<-NULL
    pihvec<-NULL
    
    ### grid level results
    protG<-NULL    
    predvecG<-NULL
    dgridvec<-NULL
    
    ### dose/symbol level results
    protDS<-NULL
    ymvecDS<-NULL
    dosevecDS<-NULL
    predvalDS<-NULL
    symDS<-NULL
    
    for(k in nstart:nprot){
      
      ### protocol specific samples
      yS<-y[prot==k]
      doseS<-dose[prot==k]
      countS<-count[prot==k]
      protS<-prot[prot==k]
      symS<-levels(symbol)[symbol[prot==k]]
      doselev<-sort(unique(doseS))
      ndose<-length(doselev)
      symlev<-sort(unique(symS))
      nsym<-length(symlev)
      ym <- NULL
      dvec <- NULL
      symvec <- NULL
      nds <- 0
      ### dose/symbol group means and se's
      for(i in 1:ndose){
        for(j in 1:nsym){
          ind<-((doseS==doselev[i]) & (symS==symlev[j]))
          ng0<-sum(ind)
          if(ng0>0){
            ymhold<-weighted.mean(yS[ind],countS[ind])
            if(isTRUE(plotDif & i==1))ym0<-ymhold
            if(plotDif)ymhold<-ymhold-ym0
            ym<-c(ym,ymhold)
            dvec<-c(dvec,doselev[i])
            symvec<-c(symvec,symlev[j])
            nds<-nds+1
          }
        }
      }
      
			if(bayes){      
	      ### at each unique dose level for protocol
	      seout<-predict(x,doselev,int=k,clev=clev)  
	      if(!plotDif){
	        fitp<-seout$predMed
	        cil<-seout$lb     ### ci bounds 
	        cih<-seout$ub
	      }else{
	        fitp<-seout$fitdifMed
	        cil<-seout$lbdif     ### ci bounds 
	        cih<-seout$ubdif
	      }
	      
	      ### predictive intervals 
	      nsub<-apply(tapply(countS,list(doseS,symS),sum),1,min,na.rm=TRUE)
	      pmean<-seout$simResp   
	      sigsim<-seout$sigsim
	      nsim<-nrow(pmean)
	      pil<-numeric(ndose)
	      pih<-numeric(ndose)
	      if(dimFit){
	      	ppred<-matrix(numeric(nsim*ndose),ncol=ndose)
	      	if(binary){
	      		ppred<-plogis(qlogis(pmean)+rmvnorm(nsim,rep(0,dimFit),vcest))
	      	}else ppred<-pmean+rmvnorm(nsim,rep(0,dimFit),vcest)
	      	
	      	for(i in 1:ndose){
			      if(i==1)pp0<-ppred[,1]
			      if(plotDif)phold<-ppred[,i]-pp0	else phold<-ppred[,i]
			      pil[i]<-quantile(phold,probs=clevlow,na.rm=TRUE)
			      pih[i]<-quantile(phold,probs=clevup,na.rm=TRUE)      			
	      	}
	      }else{
		      for(i in 1:ndose){
			      if(!binary){ ppred<-rnorm(nsim,pmean[,i],sigsim/sqrt(nsub[i]))
			      }else{ 
			      	ppred<-rbinom(nsim,nsub[i],pmean[,i])/nsub[i]
			      }
			      if(i==1)pp0<-ppred
			      if(plotDif)ppred<-ppred-pp0	
			      pil[i]<-quantile(ppred,probs=clevlow,na.rm=TRUE)
			      pih[i]<-quantile(ppred,probs=clevup,na.rm=TRUE)
		      }
	      }
	      
	      ### predicted values on grid
	      predvals<-predict(x,dgrid,int=k,clev=clev)$predMed
	      ### predicted values for each unique dose/symbol group combination
	      predvalSym<-predict(x,dvec,int=k,clev=clev)$predMed
	      if(plotDif){
	        predvals<-predvals-predvals[1]
	        predvalSym<-predvalSym-predvalSym[1]
	      }
	      
			}else{
				
				### at each unique dose level for protocol
				seout<-predict(x,doselev,int=k,clev=clev)  
				if(!plotDif){
					fitp<-seout$pred
					se<-seout$se
					cil<-seout$lb     ### ci bounds 
					cih<-seout$ub
				}else{
					ppred<-seout$pred   ### save for use in predictive se below
					fitp<-seout$fitdif
					se<-seout$sedif
					cil<-seout$lbdif     ### ci bounds 
					cih<-seout$ubdif
				}
				
				### predictive intervals 
				nsub<-apply(tapply(countS,list(doseS,symS),sum),1,min,na.rm=TRUE)
				if(!binary){
					if(!plotDif){
						sepred<-sqrt( se^2 + (residSD^2)/nsub )
						pil<-fitp-cadj*sepred
						pih<-fitp+cadj*sepred
					}else{
						sepred<-sqrt( se^2 + (residSD^2)/nsub + (residSD^2)/nsub[1] )
						sepred[1]<-0                  ### pbo-pbo
						pil<-fitp-cadj*sepred
						pih<-fitp+cadj*sepred
					}
				}else{
					if(!plotDif){
						sepred<-sqrt( (se/(fitp*(1-fitp)))^2 + 1/(nsub*fitp*(1-fitp)) ) 
						pil<-plogis(qlogis(fitp)-cadj*sepred)
						pih<-plogis(qlogis(fitp)+cadj*sepred)
					}else{
						### normal approximation to binomial for differences
						sepred<-sqrt(se^2+ppred*(1-ppred)/nsub + ppred[1]*(1-ppred[1])/nsub[1])
						sepred[1]<-0   #pbo-pbo
						pil<-fitp-cadj*sepred
						pih<-fitp+cadj*sepred
					}
				}
				
				### predicted values on grid
				predvals<-predict(x,dgrid,int=k,clev=clev)$pred
				### predicted values for each unique dose/symbol group combination
				predvalSym<-predict(x,dvec,int=k,clev=clev)$pred
				if(plotDif){
					predvals<-predvals-predvals[1]
					predvalSym<-predvalSym-predvalSym[1]
				}
				
			}
      
      #### accumulate results for protocol
      ### dose level results
      
      protD<-c(protD,rep(k,ndose))
      doseLevVec<-c(doseLevVec,doselev)
      fitvec<-c(fitvec,fitp)
      cilvec<-c(cilvec,cil)
      cihvec<-c(cihvec,cih)
      pilvec<-c(pilvec,pil)
      pihvec<-c(pihvec,pih)
      
      ### grid level
      protG<-c(protG,rep(k,ngridsup))
      predvecG<-c(predvecG,predvals)
      dgridvec<-c(dgridvec,dgrid)
      
      ### dose/symbol group level
      protDS<-c(protDS,rep(k,nds))
      ymvecDS<-c(ymvecDS,ym)
      dosevecDS<-c(dosevecDS,dvec)
      symDS<-c(symDS,symvec)
      predvalDS<-c(predvalDS,predvalSym)
    } 
    residuals<-ymvecDS-predvalDS
    
    symDS<-factor(symDS)  
    nlevsym<-length(levels(symDS))
    nshape<-length(symbolShape)
    if(nshape!=nlevsym && nshape==1)symbolShape<-rep(symbolShape,nlevsym)
    ncolor<-length(symbolColor)
    if(ncolor!=nlevsym && ncolor==1)symbolColor<-rep(symbolColor,nlevsym)
    
    if(int)protlab<-protlab[int]
    
    if(plotResid & !logScale){
      lplot<-ggplot(data.frame(dosevecDS,residuals,symDS,protDS=factor(protDS,labels=protlab)),
                    aes(x=dosevecDS,y=residuals))
      lplot<-lplot+geom_point(aes(shape=symDS,color=symDS),size=symbolSize)        
      lplot<-lplot+scale_color_manual(name=symbolLabel,values=symbolColor)
      lplot<-lplot+scale_shape_manual(name=symbolLabel,values=symbolShape)
      lplot<-lplot + geom_hline(yintercept=0,linetype=2)
      if(length(unique(prot))>1){
      	if(is.null(ncol)) ncol<-min(nprot,3)
      	lplot<-lplot+facet_wrap(~protDS,ncol=ncol) 
      }
      lplot<-lplot+ylab(ylab) + xlab(xlab) + ggplot2::theme_bw()   
      if(!is.null(ylim)){lplot<-lplot+coord_cartesian(xlim=xlim,ylim=ylim)
      }else lplot<-lplot+coord_cartesian(xlim=xlim)
    }
    else if(plotResid & logScale){
      x0 <- dosevecDS
      if(sum(x0==0)){
        xtemp <- sort(unique(dosevecDS))[2]^2/sort(unique(dosevecDS))[3]
        dosevecDSlog <- dosevecDS 
        dosevecDSlog[dosevecDS==0] <- xtemp 
        xlimlog <- c(log(xtemp)-(log(xlim[2])-log(xtemp))/10, 
                     log(xlim[2])+(log(xlim[2])-log(xtemp))/10)        
      }else{
        xtemp <- sort(unique(dosevecDS))[1]
        dosevecDSlog <- dosevecDS
        xlimlog <- c(log(xtemp)-(log(xlim[2])-log(xtemp))/10, 
                     log(xlim[2])+(log(xlim[2])-log(xtemp))/10)
      }
      lplot<-ggplot(data.frame(dosevecDSlog,residuals,symDS,protDS=factor(protDS,labels=protlab)),
                    aes(x=log(dosevecDSlog),y=residuals))
      lplot<-lplot+geom_point(aes(shape=symDS,color=symDS),size=symbolSize)        
      lplot<-lplot+scale_color_manual(name=symbolLabel,values=symbolColor)
      lplot<-lplot+scale_shape_manual(name=symbolLabel,values=symbolShape)
      lplot<-lplot + geom_hline(yintercept=0,linetype=2)
      if(is.null(xat))  lplot <- lplot + scale_x_continuous(breaks=log(sort(unique(dosevecDSlog))),
                                          labels=sort(unique(dosevecDS))) 
      if(length(unique(prot))>1){
      	if(is.null(ncol))ncol<-min(nprot,3)
      	lplot<-lplot+facet_wrap(~protDS,ncol=ncol) 
      }
      if(!is.null(ylim)){lplot<-lplot+coord_cartesian(xlim=xlimlog,ylim=ylim)
      }else lplot<-lplot+coord_cartesian(xlim=xlimlog)
    } 
    
    else if(!plotResid & !logScale){
      
      if(is.null(bwidth)){werrbar<-min(diff(sort(unique(dosevecDS))))*(0.4)
      }else werrbar<-bwidth
      lplot<-ggplot(data.frame(cilvec,cihvec,pilvec,pihvec,
                               doseLevVec,protD=factor(protD,labels=protlab)))
      lplot<-lplot+scale_color_manual(name=symbolLabel,values=symbolColor)
      lplot<-lplot+scale_shape_manual(name=symbolLabel,values=symbolShape)
      if(length(unique(prot))>1){
      	if(is.null(ncol))ncol<-min(nprot,5)
      	lplot<-lplot+facet_wrap(~protD,ncol=ncol) 
      }
      if(predict)lplot<-lplot+geom_errorbar(aes(x=doseLevVec,ymax=pihvec,ymin=pilvec),size=1.1,
                                            color='grey',width=werrbar)
      if(plotci)lplot<-lplot+geom_errorbar(aes(x=doseLevVec,ymax=cihvec,ymin=cilvec),width=0,size=1.1,color='black')
      lplot<-lplot+geom_line(data=data.frame(dgridvec=dgridvec,predvecG=predvecG,
                                             protD=factor(protG,labels=protlab)),
                             aes(x=dgridvec,y=predvecG),color='black',size=1.1, ...)   
      lplot<-lplot+geom_point(data=data.frame(ymvecDS,dosevecDS,symDS,protD=factor(protDS,labels=protlab)),
                              aes(x=dosevecDS,y=ymvecDS,shape=symDS,color=symDS),
                              size=symbolSize)
      lplot<-lplot+ylab(ylab) + xlab(xlab) + ggplot2::theme_bw()   
      if(!is.null(ylim)){lplot<-lplot+coord_cartesian(xlim=xlim,ylim=ylim)
      }else lplot<-lplot+coord_cartesian(xlim=xlim)
      
    }
    else if(!plotResid & logScale){
      x0 <- doseLevVec 
      if(sum(x0==0)){
        
        xtemp <- sort(unique(doseLevVec))[2]^2/sort(unique(doseLevVec))[3]
        doseLevVeclog <- doseLevVec
        dgridveclog <- dgridvec
        dosevecDSlog <- dosevecDS 
        doseLevVeclog[doseLevVec==0] <- dgridveclog[dgridvec==0] <- dosevecDSlog[dosevecDS==0] <- xtemp
        bench_doseLevVeclog <-  sort(unique(doseLevVeclog))
        xlimlog <- c(log(xtemp)-(log(xlim[2])-log(xtemp))/10, 
                     log(xlim[2])+(log(xlim[2])-log(xtemp))/10)
        werrbarlog <- min(diff(sort(unique(log(dosevecDSlog)))))*(0.4)
        
        ##do not include dose=0 for plotdiff
        if(plotDif){
          mindose <- min(dosevecDS[dosevecDS>0])
          index <- which(doseLevVeclog>=mindose)
          indexgrid <- which(dgridveclog>=mindose)
          doseLevVec <- doseLevVec[index]
          doseLevVeclog <- doseLevVeclog[index]
          dgridveclog <-dgridveclog[indexgrid] 
          predvecG <- predvecG[indexgrid]   
          protG <- protG[indexgrid]          
          dosevecDSlog <- dosevecDSlog[index]
          cilvec <- cilvec[index]
          cihvec <- cihvec[index] 
          pilvec <- pilvec[index]
          pihvec <- pihvec[index]
          ymvecDS <- ymvecDS[index]
          protD <- protD[index]
          protDS <- protDS[index] 
          symDS <- symDS[index]
          fitvec <- fitvec[index]     
          xlimlog <- c(log(mindose)-(log(xlim[2])-log(mindose))/10, 
                       log(xlim[2])+(log(xlim[2])-log(mindose))/10)     
        }
        
        lplot<-ggplot(data.frame(cilvec,cihvec,pilvec,pihvec,
                                 doseLevVec, doseLevVeclog, protD=factor(protD,labels=protlab)))
        lplot<-lplot+scale_color_manual(name=symbolLabel,values=symbolColor)
        lplot<-lplot+scale_shape_manual(name=symbolLabel,values=symbolShape)
        if(length(unique(prot))>1){
      	if(is.null(ncol))ncol<-min(nprot,5)
        	lplot<-lplot+facet_wrap(~protD,ncol=ncol)
        }
        
        
        if(predict)lplot<-lplot+geom_errorbar(aes(x=log(doseLevVeclog),ymax=pihvec,ymin=pilvec),size=1.1,
                                              color='grey',width=werrbarlog)
        
        if(plotci)lplot<-lplot+geom_errorbar(aes(x=log(doseLevVeclog),ymax=cihvec,ymin=cilvec),width=0,size=1.1,
                                             color='black')
        
        data=data.frame(dgridveclog,predvecG=predvecG,protD=factor(protG,labels=protlab))
        data1 <- subset(data, data$dgridveclog < bench_doseLevVeclog [2] & data$dgridveclog >= xtemp)
        data2 <- subset(data, data$dgridveclog >= bench_doseLevVeclog [2])
        lplot<-lplot+geom_line(data=data1, aes(x=log(dgridveclog),y=predvecG),color='black',size=1.1, linetype="dashed") 
        lplot<-lplot+geom_line(data=data2, aes(x=log(dgridveclog),y=predvecG),color='black',size=1.1, linetype="solid") 
        
        
        
        lplot<-lplot+geom_point(data=data.frame(ymvecDS,dosevecDSlog,symDS,protD=factor(protDS,labels=protlab)),
                                aes(x=log(dosevecDSlog),y=ymvecDS,shape=symDS,color=symDS),
                                size=symbolSize)
        if(is.null(xat)) lplot <- lplot + scale_x_continuous(breaks=log(sort(unique(doseLevVeclog))),
                                            labels=sort(unique(doseLevVec)))
        if(!is.null(ylim)){lplot<-lplot+coord_cartesian(xlim=xlimlog,ylim=ylim)
        }else lplot<-lplot+coord_cartesian(xlim=xlimlog)
      }else{
        
        xtemp <- sort(unique(doseLevVec))[1]^2/sort(unique(doseLevVec))[2]
        doseLevVeclog <- doseLevVec
        dgridveclog <- dgridvec
        dosevecDSlog <- dosevecDS 
        dgridveclog[dgridvec==0] <- xtemp
        bench_doseLevVeclog <-  sort(unique(doseLevVeclog))        
        xlimlog <- c(log(xtemp)-(log(xlim[2])-log(xtemp))/10, 
                     log(xlim[2])+(log(xlim[2])-log(xtemp))/10)
        werrbarlog <- min(diff(sort(unique(log(dosevecDSlog)))))*(0.4)
        
        ## if plotDif is true, the plot starts from the minimum non-zero dose
        ## as dose 0 always corresponds to Dif=0
        if(plotDif){
          mindose <- min(dosevecDS[dosevecDS>0])
          index <- which(doseLevVeclog>=mindose)
          indexgrid <- which(dgridveclog>=mindose)
          doseLevVec <- doseLevVec[index]        
          doseLevVeclog <- doseLevVeclog[index]
          dgridveclog <-dgridveclog[indexgrid] 
          predvecG <- predvecG[indexgrid] 
          protG <- protG[indexgrid]          
          dosevecDSlog <- dosevecDSlog[index]
          cilvec <- cilvec[index]
          cihvec <- cihvec[index] 
          pilvec <- pilvec[index]
          pihvec <- pihvec[index]
          ymvecDS <- ymvecDS[index]
          protD <- protD[index]
          protDS <- protDS[index]
          symDS <- symDS[index]  
          fitvec <- fitvec[index]
          xlimlog <- c(log(mindose)-(log(xlim[2])-log(mindose))/10, 
                       log(xlim[2])+(log(xlim[2])-log(mindose))/10)          
        }        
        lplot<-ggplot(data.frame(cilvec,cihvec,pilvec,pihvec,
                                 doseLevVec,doseLevVeclog,protD=factor(protD,labels=protlab)))
        lplot<-lplot+scale_color_manual(name=symbolLabel,values=symbolColor)
        lplot<-lplot+scale_shape_manual(name=symbolLabel,values=symbolShape)
        if(length(unique(prot))>1){
      		if(is.null(ncol))ncol<-min(nprot,5)
        	lplot<-lplot+facet_wrap(~protD,ncol=ncol) 
        }
        if(predict)lplot<-lplot+geom_errorbar(aes(x=log(doseLevVeclog),ymax=pihvec,ymin=pilvec),size=1.1,
                                              color='grey',width=werrbarlog)
        if(plotci)lplot<-lplot+geom_errorbar(aes(x=log(doseLevVeclog),ymax=cihvec,ymin=cilvec),width=0,size=1.1,color='black')
        
        data=data.frame(dgridveclog,predvecG=predvecG,protD=factor(protG,labels=protlab))
        data1 <- subset(data, data$dgridveclog < bench_doseLevVeclog[1] & data$dgridveclog >= xtemp)
        data2 <- subset(data, data$dgridveclog >= bench_doseLevVeclog[1])
        
        lplot<-lplot+geom_line(data=data1, aes(x=log(dgridveclog),y=predvecG),color='black',size=1.1, linetype="dashed") 
        lplot<-lplot+geom_line(data=data2, aes(x=log(dgridveclog),y=predvecG),color='black',size=1.1, linetype="solid") 
        
        lplot<-lplot+geom_point(data=data.frame(ymvecDS,dosevecDSlog,symDS,protD=factor(protDS,labels=protlab)),
                                aes(x=log(dosevecDSlog),y=ymvecDS,shape=symDS,color=symDS),
                                size=symbolSize)
        if(is.null(xat)) lplot <- lplot + scale_x_continuous(breaks=log(sort(unique(doseLevVeclog))),
                                            labels=sort(unique(doseLevVec))) 
        if(!is.null(ylim)){lplot<-lplot+coord_cartesian(xlim=xlimlog,ylim=ylim)
        }else lplot<-lplot+coord_cartesian(xlim=xlimlog)
        
      }
    }
    
    lplot <- lplot +ylab(ylab) + xlab(xlab) + ggplot2::theme_bw()   
    ## remove the vertical grid lines
    lplot <- lplot+ ggplot2::theme(panel.grid.major.x=element_blank(),
    											panel.grid.minor.x=element_blank(),
    											panel.grid.major.y=element_line(size=0.1))   
    
    if(nolegend)lplot<-lplot + ggplot2::theme(legend.position = "none")

    
    if(!is.null(xat)){
      if(!logScale)      lplot <- lplot + scale_x_continuous(breaks=xat, 
                                          labels =xat)
      if(logScale){
        xatbench <- xat
        xat[xat==0] <- doselev[2]^2/doselev[3]
        lplot <- lplot + scale_x_continuous(breaks=log(xat), 
                                            labels =xatbench)           
      }         
    }
    
    if(plot){
      print(lplot)
    }
    
    return(invisible(list(lplot=lplot,plotdata=cbind(prot=protD,dose=doseLevVec,fit=fitvec,
                                                     cil=cilvec,cih=cihvec,
                                                     pil=pilvec,pih=pihvec))))
  }



