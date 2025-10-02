"update.emaxsimobj" <-
function(object,new.parm,modType=object$modType,...){
    
    if(length(new.parm)!=modType)stop("modType does not match
parameter length")

	binary<-object$binary
     dose<-object$dose
     y<-object$y
     
     if(is.null(names(new.parm))){
        if(modType==4){names(new.parm)<-c("led50","lambda","emax","e0")
        }else names(new.parm)<-c("led50","emax","e0")
     }else if(!all(names(new.parm)%in%c("led50","emax","e0"))){
              stop("Invalid names for new.parm")
     }

     simout<-emaxalt(y,dose,modType,binary=binary,iparm=new.parm,
                     ed50cutoff=object$ed50cutoff,
                     ed50lowcutoff=object$ed50lowcutoff,switchMod=object$switchMod)

     object$init<-simout$Sparm
     object$fitType<-simout$fitType
     object$fitpred <-simout$fitpred
     object$sepred <-simout$sepred
     object$sedif <- simout$sedif
     object$est4<-simout$est4
     object$est3<-simout$est3
     object$estA<-simout$estA
     object$vc<-simout$vc
	 object$residSD<-simout$residSD
                
     return(object)
}

