"emaxsimB" <-
function(nsim, genObj, prior, modType=3,binary=FALSE,seed=12357,
				 check=FALSE,nproc=parallel::detectCores(),
				 negEmax=FALSE,ed50contr=NULL, lambdacontr=NULL,testMods=NULL,
         idmax=length(doselev),mcmc=mcmc.control(),customCode=NULL,
				 customParms=NULL,
         description="")
	{
	####  assumes dose levels are sorted on input
	####  placebo (comparator) is in the first position
	####  primary test dose is in position idmax

  if(! modType %in% c(3,4))stop("modType must be 3 or 4")
  if(length(ed50contr)!=length(lambdacontr))stop('The number of ED50 and Lambda defining contrasts must be equal')
	
	if(exists('.Random.seed'))save.seed<-.Random.seed
	save.rng<-RNGkind()[1]
	
	on.exit( RNGkind(save.rng))
	on.exit(if(exists('save.seed')).Random.seed<<-save.seed,add=TRUE)
	
  ### extract design parameters from genObj
	n<-genObj$genP$n
	doselev<-genObj$genP$doselev ## doselev should be sorted/unique
  Ndose<-length(doselev)
  if(!binary)ddf<-sum(n)-Ndose else ddf<-Inf
	nfrac<-n/sum(n)              ##allocate fractional obs for stability
	nfrac2<-0.5*nfrac
	dose<-genObj$genP$dose
	
	estan<-selEstan(modType=modType,binary=binary)

	 ### set up emax-model contrasts for null hypothesis test
	contMat<-NULL
  if(missing(testMods)){
    if(is.null(ed50contr)){
        if(Ndose<=4){
            ed50contr<-c((doselev[1]+doselev[2])/2,
                                  (doselev[Ndose-1]+doselev[Ndose])/2)       
            lambdacontr<-rep(1,2)
        }else{
            ed50contr<-c((doselev[1]+doselev[2])/2,median(doselev),
                                  (doselev[Ndose-1]+doselev[Ndose])/2)
            lambdacontr<-rep(1,3)
        }
    }
    parmscontr<-cbind(ed50contr,lambdacontr)
    testMods<-Mods(sigEmax=parmscontr,doses=doselev, placEff = 0, 
                   maxEff=1-2*(negEmax))
  }
	if(!binary) contMat <-optContr(testMods,w=n) 

  ### simulation result holders
	if(modType==3){nparm<-3
    }else {nparm<-4}

	### posterior intervals to compute and store
 	llev<-c(0.025,0.05,0.1)
	ulev<-c(0.975,0.95,0.9)
      

  ### set up independent stream of random numbers for
  ### each simulation iteration.
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  rseed<-matrix(integer(nsim*7),ncol=7)
  rseed[1,]<-as.integer(.Random.seed)
  for(i in 2:nsim){
  	rseed[i,]<-nextRNGStream(rseed[i-1,])
  }
  
  if(isTRUE(.Platform$OS.type=='unix') && missing(nproc))
  	stop('nproc must be specified on multi-user unix machines')
 
	if(check)nproc<-1 else{
		if(nproc>detectCores())stop("The number of processes requested exceeds the number of processors available.")
			if(nproc>nsim){
				warning(paste('The number of processors requested exceeded the number of simulations\n', 
								'This iikely a mistake.  nproc set to 1'))
				nproc<-1
			}
	}
  
  if(mcmc$chains>1 && !check)stop('The number of chains should be 1 except when testing.')

  ### set up indices for computing consecutive blocks
  ### of simulations
  nblock<-as.integer(trunc(nsim/nproc))
  nleft<-nsim%%nproc
  indmat<-matrix(integer(nproc*2),ncol=2)
  indmat[1,1]<-1
  indmat[1,2]<-nblock+1*(1<=nleft)
  if(nproc>=2){
	  for (i in 2:nproc){
	  	indmat[i,1]<-indmat[i-1,2]+1
	  	indmat[i,2]<-indmat[i,1]+nblock-1 + 1*(i<=nleft)
	  }
  }
  
	inlist<-list(indmat=indmat,rseed=rseed,Ndose=Ndose,dose=dose,
						ddf=ddf,doselev=doselev,nparm=nparm,modType=modType,
						binary=binary,genObj=genObj,testMods=testMods,
						contMat=contMat,negEmax=negEmax,
						check=check,estan=estan,prior=prior,mcmc=mcmc,
						customCode=customCode,customParms=customParms,
						nfrac=nfrac,nfrac2=nfrac2,n=n,
						ulev=ulev,llev=llev)

	if(nproc==1){
		simout<-simrepB(1,inlist)
		if(check)return(simout)
		simout<-list(simout)
	}else{
		cl<-makeCluster(nproc)
		registerDoParallel(cl)	
		simout<-foreach(j=1:nproc, .packages='clinDR') %dopar%{
			simrepB(j,inlist)
		}
		stopCluster(cl)
	}  

 	####################################
	### assign output to matrices/vectors

	predpop <- matrix(rep(NA, nsim * Ndose), ncol = Ndose)
	fitpredv <- matrix(rep(NA, nsim * Ndose), ncol = Ndose)
	sepredv <- matrix(rep(NA, nsim * Ndose), ncol = Ndose)
	sedifv <- matrix(rep(NA, nsim * Ndose), ncol = Ndose)
	mv <- matrix(rep(NA, nsim * Ndose), ncol = Ndose)
	msSat<-rep(NA,nsim)
  pVal<-rep(NA,nsim)
  gofP<-rep(NA,nsim)
  selContrast<-rep(NA,nsim)
  colnames(predpop)<-doselev
  colnames(fitpredv)<-doselev
  colnames(sepredv)<-doselev
  colnames(sedifv)<-doselev
  colnames(mv)<-doselev

	residSD<-numeric(nsim)
	est<-matrix( rep(NA,nsim*nparm),ncol=nparm )
  if(modType==3){
   colnames(est)<-c("led50","emax","e0")
  }else{
   colnames(est)<-c("led50","lambda","emax","e0")
  }
	if(!binary){
		sdv <- matrix(rep(NA, nsim * Ndose), ncol = Ndose)
	}else sdv<-NULL

	nlev<-length(llev)
	nd1<-Ndose-1
	lb<-array(numeric(nd1*nsim*nlev),dim=c(nd1,nsim,nlev))
	ub<-array(numeric(nd1*nsim*nlev),dim=c(nd1,nsim,nlev))
	dimnames(lb)<-list(doselev[2:Ndose],1:nsim,llev)
	dimnames(ub)<-list(doselev[2:Ndose],1:nsim,ulev)
	
  pop<-NULL
  popSD<-NULL
	if(is.null(customCode))customOut<-NULL
		else customOut<-vector("list",nsim)

	for(j in 1:nproc){
		predpop[indmat[j,1]:indmat[j,2],]<-simout[[j]]$predpop	
		fitpredv[indmat[j,1]:indmat[j,2],]<-simout[[j]]$fitpredv
		sepredv[indmat[j,1]:indmat[j,2],]<-simout[[j]]$sepredv
		sedifv[indmat[j,1]:indmat[j,2],]<-simout[[j]]$sedifv
		mv[indmat[j,1]:indmat[j,2],]<-simout[[j]]$mv
		sdv[indmat[j,1]:indmat[j,2],]<-simout[[j]]$sdv
		msSat[indmat[j,1]:indmat[j,2]]<-simout[[j]]$msSat
		pVal[indmat[j,1]:indmat[j,2]]<-simout[[j]]$pVal
		gofP[indmat[j,1]:indmat[j,2]]<-simout[[j]]$gofP
		selContrast[indmat[j,1]:indmat[j,2]]<-simout[[j]]$selContrast
		residSD[indmat[j,1]:indmat[j,2]]<-simout[[j]]$residSD
		est[indmat[j,1]:indmat[j,2],]<-simout[[j]]$est
		lb[,indmat[j,1]:indmat[j,2],]<-simout[[j]]$lb
		ub[,indmat[j,1]:indmat[j,2],]<-simout[[j]]$ub
		customOut[indmat[j,1]:indmat[j,2]]<-simout[[j]]$customOut
		pop<-rbind(pop,simout[[j]]$pop)
		popSD<-c(popSD,simout[[j]]$popSD)		
	}
  
	return(structure(list(description=description,
				binary=binary,modType=modType,genObj=genObj, 
        pop=pop,popSD=popSD,mcmc=mcmc,prior=prior,
				est=est,residSD=residSD,
        pVal=pVal,selContrast=selContrast,
				testMods=testMods,gofP=gofP,
        negEmax=negEmax,
        predpop=predpop,        
        mv = mv, sdv = sdv, msSat=msSat, fitpredv = fitpredv,
        sepredv = sepredv, sedifv = sedifv, 
				lb=lb,ub=ub,
        rseed=rseed, idmax=idmax,customOut=customOut
        ),class="emaxsimB") )

}

simrepB<-function(j,inlist)
	{

	indmat<-inlist$indmat
	rseed<-inlist$rseed
	Ndose<-inlist$Ndose
	dose<-inlist$dose
	ddf<-inlist$ddf
	doselev<-inlist$doselev
	nparm<-inlist$nparm
	modType<-inlist$modType
	binary<-inlist$binary
	genObj<-inlist$genObj
	testMods<-inlist$testMods
	contMat<-inlist$contMat
	negEmax<-inlist$negEmax
	check<-inlist$check
	estan<-inlist$estan
	prior<-inlist$prior
	mcmc<-inlist$mcmc
	customCode<-inlist$customCode
	customParms<-inlist$customParms
	nfrac<-inlist$nfrac
	nfrac2<-inlist$nfrac2
	n<-inlist$n
	llev<-inlist$llev
	ulev<-inlist$ulev
	
 	nrep<-indmat[j,2]-indmat[j,1]+1 
 	
 	### set up input dose variable
 	if(binary)din<-c(doselev,doselev) else din<-doselev
 	if(binary)sigsim<-1  ### default placeholder
 	
	predpop <- matrix(rep(NA, nrep * Ndose), ncol = Ndose)
	fitpredv <- matrix(rep(NA, nrep * Ndose), ncol = Ndose)
	sepredv <- matrix(rep(NA, nrep * Ndose), ncol = Ndose)
	sedifv <- matrix(rep(NA, nrep * Ndose), ncol = Ndose)
	mv <- matrix(rep(NA, nrep * Ndose), ncol = Ndose)
	msSat<-rep(NA,nrep)
  pVal<-rep(NA,nrep)
  selContrast<-rep(NA,nrep)
  colnames(predpop)<-doselev
  colnames(fitpredv)<-doselev
  colnames(sepredv)<-doselev
  colnames(sedifv)<-doselev
  colnames(mv)<-doselev

	residSD<-numeric(nrep)
	gofP<-numeric(nrep)
	est<-matrix( rep(NA,nrep*nparm),ncol=nparm )
  if(modType==3){
   colnames(est)<-c("led50","emax","e0")
  }else{
   colnames(est)<-c("led50","lambda","emax","e0")
  }
	if(!binary){
		sdv <- matrix(rep(NA, nrep * Ndose), ncol = Ndose)
	}else sdv<-NULL

	nlev<-length(llev)
	nd1<-Ndose-1
	lb<-array(numeric(nd1*nrep*nlev),dim=c(nd1,nrep,nlev))
	ub<-array(numeric(nd1*nrep*nlev),dim=c(nd1,nrep,nlev))
	dimnames(lb)<-list(doselev[2:Ndose],1:nrep,llev)
	dimnames(ub)<-list(doselev[2:Ndose],1:nrep,ulev)
	
  pop<-NULL
  popSD<-NULL
  
  if(is.null(customCode))customOut<-NULL
  	 else customOut<-vector("list",nrep)
  
  if(!negEmax)trend<-'positive' else trend<-'negative'

	for(i in 1:nrep) {
		ioff<-i+indmat[j,1]-1

    ### generate simulated data set by calling
    ### the function in genObj
		.Random.seed<<-rseed[ioff,]

    ### simulate design and data
		gendat <- genObj$genFun(genObj$genP)
		predpop[i,]<-gendat[['meanlev']]
    pop<-rbind(pop,gendat[['parm']])
    popSD<-c(popSD,gendat[['resSD']])
		y<-gendat[['y']]

		#### set up normal approximation to logit(phat) if
		#### binary data
		####
		if(binary){
			phat<-tapply(y,dose,sum)
			phat<-(phat+nfrac2)/(n+nfrac)            ###shrink to .5 for stability 
			yhat<-qlogis(phat)
			V<-diag(1/(phat*(1-phat)*n))
			contMat <-optContr(testMods,S=V)
		}else{
			yhat<-tapply(y,dose,mean)
			ms<-summary(lm(y~factor(dose)))$sigma
			msSat[i]<-ms^2
			V<-diag(msSat[i]/n)
		}

    ### compute p-value for preliminary test
    holdC<-MCTtest(doselev, yhat, contMat = contMat,
						  S=V,df=ddf,type='general')$tStat
    holdP<-attr(holdC,'pVal')
    orderh<-order(holdC)
		ncontr<-length(orderh)
    pVal[i]<-holdP[orderh[ncontr]]
    selContrast[i]<-orderh[ncontr]
    
    ### format data as counts if binary
    if(binary){
    	cin<-tapply(y,dose,sum)
    	cin<-	c(cin,n-cin)
    	yin<-c(rep(1,Ndose),rep(0,Ndose))
    }else{
    	cin<-n
    	yin<-yhat
    } 
    bfit<-fitEmaxB(yin,din,prior,modType,count=cin,binary=binary,
    							 msSat=msSat[i],mcmc=mcmc,estan=estan,
    							 diagnostics=check,nproc=1)
    
    ### return mcmc for preliminary checking
    if(check){
    	if(binary){
	    	return(list(estanfit=bfit$estanfit,parms=coef(bfit),pVal=pVal[i],dose=dose,y=y))
    	}else{
	    	return(list(estanfit=bfit$estanfit,parms=coef(bfit),residSD=sigma(bfit),pVal=pVal[i],dose=dose,y=y))
    	}
    }

    ### assign simulation output
    mv[i,  ]  <- tapply(y,dose,mean)
		if(!binary){
			sdv[i,  ] <- sqrt(tapply(y,dose,var))
		}
    
    simout<-predict(bfit,doselev)
    fitpredv[i,]<-simout$pred
    sepredv[i,  ] <-simout$se
    sedifv[i,  ] <-simout$sedif
    simdif<-simout$simResp[,2:Ndose]-simout$simResp[,1]
    for(k in 1:nlev){
    	qtile<-apply(simdif,2,quantile,c(llev[k],ulev[k]))
    	lb[,i,k]<-qtile[1,]
    	ub[,i,k]<-qtile[2,]
    }
    
    ### extract generated parameters
  	parms<-coef(bfit)
		if(!binary){
			sigsim<- sigma(bfit)
			residSD[i]<-median(sigsim)
		}
    est[i,]<- apply(parms,2,median)
    
    ### compute gof test
    gofP[i]<-checkMonoEmax(yin,din,parms,sigsim^2,cin,
    											 trend=trend,logit=binary)
    
    ### execute customized code
    if(!is.null(customCode)){
    	
    	if(!is.null(customParms)){
	    	if(binary)customOut[[i]]<-customCode(parms,pVal[i],dose,y,
	    																			 customParms=customParms )	
	    	else 	customOut[[i]]<-customCode(parms,sigsim,pVal[i],dose,y,
	    																	 customParms=customParms )	
    	}else{
	    	if(binary)customOut[[i]]<-customCode(parms,pVal[i],dose,y)
	    	else 	customOut[[i]]<-customCode(parms,sigsim,pVal[i],dose,y)
    	}
    }
	}

	return(list(pop=pop,popSD=popSD,  
        est=est,residSD=residSD,
        pVal=pVal,selContrast=selContrast,
        gofP=gofP,predpop=predpop,        
        mv = mv, sdv = sdv, msSat=msSat, fitpredv = fitpredv,
        sepredv = sepredv, sedifv = sedifv, lb=lb,ub=ub,
				customOut=customOut
        ))
}

