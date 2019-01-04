"[.emaxsim" <-
function(x,i,...){
  if(exists('.Random.seed'))save.seed<-.Random.seed
  save.rng<-RNGkind()[1]
  RNGkind("L'Ecuyer-CMRG")
  .Random.seed<<-as.integer(x$rseed[i,])
    
	binary<-x$binary

  ### re-create dose
  doselev<-x$genObj$genP$doselev
  n<-x$genObj$genP$n
  dose<-rep(doselev,n)


  ### regenerate data set
  gendat <- x$genObj$genFun(x$genObj$genP)
  predpop<-gendat[['meanlev']]
  y<-gendat[['y']]

  ### convert vc to appropriate dimensional matrix
  if(x$fitType[i]=="4"){
        vc<-matrix(x$vc[i,],ncol=4)
  }else if(x$fitType[i]=="3"){
        vc<-matrix(x$vc[i,1:9],ncol=3)
  } else{ 
        n2<-ncol(x$estA)
        vc<-matrix(x$vc[i,c(1:(n2^2))],ncol=n2)
  }

  if(exists('save.seed')).Random.seed<<-save.seed   ### re-assign seed to original sequence order
  RNGkind(save.rng)
 
     return(structure(list(y=y, dose=dose, pop=x$pop[i,],
     											popSD=x$popSD[i],binary=binary,
            init=x$init[i,], est4=x$est4[i,],
            est3=x$est3[i,],estA=x$estA[i,], 
            vc=vc,residSD=x$residSD[i],bigC=x$bigC[i],
            negC=x$negC[i],modType=x$modType, 
            fitType=x$fitType[i], ed50cutoff=x$ed50cutoff, 
            ed50lowcutoff=x$ed50lowcutoff,switchMod=x$switchMod, 
            predpop=x$predpop[i,],
            dm=x$mv[i,],dsd=x$sdv[i,], fitpred = x$fitpredv[i,],
            sepred = x$sepredv[i,], sedif = x$sedifv[i,],                     
   	        pVal=x$pVal[i],selContrast=x$selContrast[i],
						idmax=x$idmax),class="emaxsimobj"))
}


"[.emaxsimB" <-
function(x,i,...){
  if(exists('.Random.seed'))save.seed<-.Random.seed
  save.rng<-RNGkind()[1]
  RNGkind("L'Ecuyer-CMRG")
  .Random.seed<<-as.integer(x$rseed[i,])
  modType<-x$modType
	binary<-x$binary
	msSat<-x$msSat[i]
	prior<-x$prior
	mcmc<-x$mcmc

  ### re-create dose
  doselev<-x$genObj$genP$doselev
  Ndose<-length(doselev)
  n<-x$genObj$genP$n
  dose<-rep(doselev,n)
  

  ### regenerate data set
  gendat <- x$genObj$genFun(x$genObj$genP)
  predpop<-gendat[['meanlev']]
  y<-gendat[['y']]
  
  ### format data as counts if binary
  if(binary){
  	cin<-tapply(y,dose,sum)
  	cin<-	c(cin,n-cin)
  	yin<-c(rep(1,Ndose),rep(0,Ndose))
  }else{
  	cin<-n
  	yin<-tapply(y,dose,mean)
  } 
  
 	### set up input dose variable
 	if(binary)din<-c(doselev,doselev) else din<-doselev
  
  bfit<-fitEmaxB(yin,din,prior,modType,count=cin,binary=binary,
   								 msSat=msSat,mcmc=mcmc,estan=NULL,diagnostics=FALSE)


  if(exists('save.seed')).Random.seed<<-save.seed   ### re-assign seed to original sequence order
  RNGkind(save.rng)
  
  return(structure(list(y=y, dose=dose, pop=x$pop[i,],
  											popSD=x$popSD[i],binary=binary,
					modType=modType, 
          predpop=x$predpop[i,],
          dm=x$mv[i,],dsd=x$sdv[i,], fitpred = x$fitpredv[i,],
          sepred = x$sepredv[i,], sedif = x$sedifv[i,],                     
					bfit=bfit,prior=prior,mcmc=mcmc,
	        pVal=x$pVal[i],selContrast=x$selContrast[i],
        	idmax=x$idmax),class="emaxsimBobj"))
}




