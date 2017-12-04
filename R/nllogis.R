"nllogis" <-
function(parms,y,dose,prot=rep(1,length(y)),count=rep(1,length(y))){
   ### if prot is specified, it must be sequential 1,2,3,...

   protid<-unique(prot)
   nprot<-length(protid)
   nsub<-length(parms)-nprot

   ### y must be coded 0/1
   if(any(y!=0 & y!=1))stop("y must be 0/1")

   b<- -parms[nsub]
   c<- 1/exp(parms[1])

   if(nsub==3){lambda<-parms[2]
   }else lambda<-1

   a<- parms[nsub+prot]-b
   evec<- a + b/(1+(c*dose)^lambda)
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

