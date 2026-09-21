"nllogis" <-
function(parms,y,dose,prot=rep(1,length(y)),count=rep(1,length(y)),
				 xbase=NULL){
   ### if prot is specified, it must be sequential 1,2,3,...
		
	nprot<-length(unique(prot))
	nparm<-length(parms)
	
	if(is.null(xbase))nbase<-0 else nbase<-ncol(xbase)
	if(nbase)bparms<-parms[(1+nparm-nbase):nparm]
	nsub<-nparm-nprot-nbase

	
	### y must be coded 0/1
	if(any(y!=0 & y!=1))stop("y must be 0/1")
	
	b<- -parms[nsub]
	c<- 1/exp(parms[1])
	
	if(nsub==3){lambda<-parms[2]
	}else lambda<-1
	
	a<- parms[nsub+prot]-b
	
	evec<- a + b/(1+(c*dose)^lambda)
	if(nbase)evec<-evec+as.vector(xbase%*%bparms)
	
	ll1<-0
	ll0<-0
	
	if(length(y==1)>0){
	   esub<-evec[y==1]
	   countsub<-count[y==1]
	   ll1<- sum(countsub*plogis(esub,log.p= TRUE))
	}
	if(length(y==0)>0){
	   esub<-evec[y==0]
	   countsub<-count[y==0]
	   ll0<- sum(countsub*plogis(-esub,log.p= TRUE))
	}
	return(-(ll1+ll0))
}

