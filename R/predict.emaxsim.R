"predict.emaxsim" <-
function(object,dose, dref=0, ...){

   vc <- object$vc
   estA <- object$estA
   est3 <- object$est3
   est4 <- object$est4
   fitType <- object$fitType
   maxE<- max(object$genObj$genP$doselev)
   binary<-object$binary
   
   nsim<-length(fitType)
   fitpredv<-matrix(numeric(nsim*length(dose)),nrow=nsim)
   fitdifv<- matrix(numeric(nsim*length(dose)),nrow=nsim)
   sepredv<- matrix(numeric(nsim*length(dose)),nrow=nsim)
   sedifv <- matrix(numeric(nsim*length(dose)),nrow=nsim)
   for(i in 1:nsim){
    if (fitType[i] == "4") {
        vcf <- matrix(vc[i, 1:16], ncol = 4)
    }
    else if (fitType[i] == "3") {
        vcf <- matrix(vc[i, 1:9], ncol = 3)
    }else vcf <- matrix(vc[i, 1:4], ncol = 2)

      out<-predictEmax(dose, estA[i, ], est3[i, ], est4[i, 
            ], vcf, fitType[i], dref, maxE, binary=binary)
      fitpredv[i,]<-out$fitpred
      fitdifv[i,]<- out$fitdif
      sepredv[i,]<- out$sepred
      sedifv[i,]<-  out$sedif
   }
   colnames(fitpredv)<-dose
   colnames(fitdifv)<-dose
   colnames(sepredv)<-dose
   colnames(sedifv)<-dose
   return(list(fitpredv=fitpredv,fitdifv=fitdifv,sepredv=sepredv,sedifv=sedifv))
}

