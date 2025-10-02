targetBeta<-
function(minval,pminV,pmaxV,maxval=1,aInit=1,bInit=1,upB=1){
    #minval lowest dose in a design
    #pminV,pmaxV  target prob less than minval and maxval(=1)
    #aInit,bInit  parameters of the beta distribution

    minfun<-function(pvec,minval,pminV,pmaxV,upB){
        t1<-pbeta(minval/upB,pvec[1],pvec[2])-pminV
        t2<-pbeta(maxval/upB,pvec[1],pvec[2])-pmaxV
        return(t1^2+t2^2)
   }
    opt<-optim(c(aInit,bInit),minfun,method="Nelder-Mead",minval=minval,
          pminV=pminV,pmax=pmaxV,upB=upB)

    if(opt$convergence==0)
    {
        aparm <- opt$par[1]
        bparm <- opt$par[2]
        return(parm=c(aparm,bparm))
    }else{
        stop("No solution found") #when there is no solution
    }
}


