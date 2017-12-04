"emaxfun" <-
function(dose,parm){
### sigmoid emax
### parm=led50,(lambda),emax,e0 

    if(is.matrix(parm)){
       if(isTRUE(all.equal(ncol(parm),3))){ftype<-3
       }else if (isTRUE(all.equal(ncol(parm),4))){ftype<-4
       }else stop('Invalid parm')
    }else if(is.vector(parm)){
       if(isTRUE(all.equal(length(parm),3))){ftype<-3
       }else if (isTRUE(all.equal(length(parm),4))){ftype<-4
       }else stop('Invalid parm')
       parm<-matrix(parm,nrow=1)
    }else stop('Invalid parm')

    dose<-matrix(dose,ncol=1)

    if(ftype==4)return(apply(X=dose,MARGIN=1,FUN=e4,led50 = parm[,1], 
        lambda = parm[,2], emax = parm[,3], e0 = parm[,4]))
    if(ftype==3)return(apply(X=dose,MARGIN=1,FUN=e3,led50 = parm[,1], 
        emax = parm[,2], e0 = parm[,3]))
}

e4<-function(dose,led50,lambda,emax,e0){
    
    return( e0+emax - emax/(1+((dose/exp(led50))^lambda)) )
}

e3<-function(dose,led50,emax,e0){
    
    return( e0+emax - emax/(1+(dose/exp(led50))) )
}

