"targetD" <- function(fit,target,modType=4,binary=FALSE){ 

   ### re-use code for other input types if emaxalt 
   ### and fitType is "3" or "4" 
   simIn<-(class(fit)=='emaxalt' || class(fit)=='emaxsimobj')
	if(simIn){
		fitType<-fit$fitType
		binary<-fit$binary
	}

   if(simIn && (fitType=='3' | fitType=='4')){ 
       fitType<-fitType
       modType<- 3*(fitType=='3')+4*(fitType=='4')
	   if(modType==3){fit<-list(fit$est3,fit$vc)
	   }else fit<-list(fit$est4,fit$vc)
       simIn<-FALSE
   }

   ### process all but 2-parameter model fits
   if(!simIn){
       ### coef=(ed50,lambda,emax,e0) or (ed50,emax,e0)
       if(class(fit)=='nls'){
          parm<-coef(fit)[1:modType]
          vfcov <- vcov(fit)[1:modType,1:modType]
       }else{
          parm<-fit[[1]][1:modType]
          vfcov<-fit[[2]][1:modType,1:modType]
       }
       
        if(modType==4){
           e0<-parm[4]
           ed50<-parm[1]
           lambda<-parm[2]
           emax<-parm[3]
        }else{
           e0<-parm[3]
           ed50<-parm[1]
           lambda<-1
           emax<-parm[2]
        }

        ed50<-exp(ed50)
		### binary risk difference converted to logistic difference
		if(binary){
			p0<-plogis(e0)
			Q<-qlogis(target+p0)-e0 
		}else Q<-target 


        if ( sign(Q)!=sign(emax) ) {
            warning("The sign of target must match the sign of emax") 
            return(c(targetDose=NA,seTD=NA))          
            }

        if ( abs(Q)>=abs(emax) ) {
            warning("The target cannot exceed emax") 
            return(c(targetDose=NA,seTD=NA))          
            }

		### intermediate computations
		trat<-Q/(emax-Q)
		tratl<-trat^(1/lambda)

        ### mle for dose meeting the target
        td<- ed50*tratl 
        names(td)<-'target dose'

        ### compute derivatives for se of td
        dered<-td
        deremax<- -td/(lambda*(emax-Q))
        derlam<- -td*log(trat)/lambda^2
		### e0
		if(!binary){
			dere0<-0
		}else{
			dQ<-p0*(1-p0)/((target+p0)*(1-(target+p0))) - 1
			dere0<-td*dQ/lambda
			dere0<-dere0*emax/(Q*(emax-Q))
		}
        
        if(modType==4){    
           L<-c(dered,derlam,deremax,dere0)
        }else { L<-c(dered,deremax,dere0)}
        setd<- sqrt(t(L)%*%vfcov%*%L)
    }else{  ### 2-parameter fit
        e0<-fit$estA[1]
        beta<-fit$estA[2]
        vfcov<-matrix(fit$vc,ncol=2)

		if(fitType=='E')f0<-1 else f0<-0

		if(!binary){
			Q<-target+beta*f0
			dQb<-f0
			dQe<-0
		}else{
			peb<-plogis(e0+beta*f0)
			tpeb<-target+peb
			Q<-qlogis(tpeb)-e0
			dQb<-f0*peb*(1-peb)/( tpeb*(1-tpeb))
			dQe<-peb*(1-peb)/( tpeb*(1-tpeb)) -1
		}

        if ( sign(Q)!=sign(beta) ) {
            warning("The sign of target must match the sign of the change in the estimated curve") 
            return(c(targetDose=NA,seTD=NA))          
        }
		if(fitType=='E' && (Q/beta<=0)){
            warning("The target cannot be achieved from the estimated exponential curve") 
            return(c(targetDose=NA,seTD=NA))
        }

		if(fitType=='L'){
			fder<-1
			td<-Q/beta
		}else if(fitType=='LL'){
			fder<-exp(-Q/beta)
			td<-exp(Q/beta)-1
		}else{
            c0<- as.numeric(names(fit$dm))[length(fit$dm)]
			fder<-Q/(c0*beta)
			td<-c0*log(Q/beta)
		}

		tdb<-(dQb/beta-Q/beta^2)/fder
		tde<-(dQe/beta)/fder


		der<-c(tde,tdb)
		setd<-sqrt(max(0,t(der)%*%vfcov%*%der))
	}
    out<-c(td,setd)
    names(out)<-c('targetDose','seTD')
    return(out)  
}


