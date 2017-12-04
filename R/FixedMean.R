"FixedMean" <-
function(n,doselev,meanlev,resSD,parm=NULL,binary=FALSE){
        dord<-order(doselev)
        doselev<-doselev[dord]
		dose<-rep(doselev,n[dord])

	if(binary){
		resSD<-NULL
	}else if(missing(resSD))stop('resSD must be specified for normal data')


###     genP will be passed to genFun
###     genP must contain n,resSD,doselev,dose + any other inputs 
###     used by genfun
        genP<-list(n=n,resSD=resSD,doselev=doselev,dose=dose,
        					 parm=parm,binary=binary,meanlev=meanlev)

		genFun <- function(obj){
			binary<-obj$binary
      meanR<-rep(obj$meanlev,obj$n)
			if(binary){
				y<-rbinom(sum(obj$n),1,meanR)
			}else y<-rnorm(sum(obj$n), mean = meanR, sd = obj$resSD)
###
###         genFun must return population mean level for each dose
###                group, parameters defining the DR curve
###                (application specific), and the response
###                data, in the order listed
			return(list(meanlev=obj$meanlev,parm=obj$parm,resSD=obj$resSD,y=y))
		}
		return( list(genP=genP,genFun=genFun)  )
}

