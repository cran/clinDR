"predict.emaxalt" <-
function(object,dose, dref=0, ...){

	binary<-object$binary

   ## select and dimension vc based on fitType
   if(object$fitType=='4'){vcf<-matrix(object$vc[1:16],ncol=4)
   }else if(object$fitType=='3'){vcf<-matrix(object$vc[1:9],ncol=3)
   }else vcf<-matrix(object$vc[1:4],ncol=2)

   ### extract max dose in design from names in mean vector
   maxE<-as.numeric(names(object$dm)[length(object$dm)])

   predE<-predictEmax(dose,object$estA,object$est3,object$est4,
                      vcf,object$fitType,dref,maxE,binary=binary)

   return(predE)
}

