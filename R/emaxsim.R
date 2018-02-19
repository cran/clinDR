"emaxsim" <-
function(nsim, genObj, modType=3,binary=FALSE,
				seed=12357,nproc=parallel::detectCores(),negEmax=FALSE,
				ed50contr=NULL, lambdacontr=NULL,testMods=NULL,
				idmax=length(doselev),iparm=rep(NA,modType), 
				ed50cutoff=2.5*max(doselev),
				ed50lowcutoff=doselev[2]/1000,switchMod= TRUE,truncLambda=6,
				description="")
	{
	####  assumes dose levels are sorted on input
	####  placebo (comparator) is in the first position
	####  primary test dose is in position idmax

	if(! modType %in% c(3,4))stop("modType must be 3 or 4")
	if(length(ed50contr)!=length(lambdacontr))stop('The number of ED50 and Lambda defining contrasts must be equal')
	
	if(exists('.Random.seed'))save.seed<-.Random.seed
	save.rng<-RNGkind()[1]
	

	### extract design parameters from genObj
	n<-genObj$genP$n
	doselev<-genObj$genP$doselev ## doselev should be sorted/unique
	 Ndose<-length(doselev)
	nfrac<-n/sum(n)              ##allocate fractional obs for stability
	nfrac2<-0.5*nfrac
	dose<-genObj$genP$dose

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
	nparmA<-2

	### set up independent stream of random numbers for
	### each simulation iteration.
	RNGkind("L'Ecuyer-CMRG")
	set.seed(seed)
	rseed<-matrix(integer(nsim*7),ncol=7)
	rseed[1,]<-as.integer(.Random.seed)
	for(i in 2:nsim){
		rseed[i,]<-nextRNGStream(rseed[i-1,])
	}
 
	### set up indices for computing consecutive blocks
	### of simulations
	if(isTRUE(.Platform$OS.type=='unix') && missing(nproc))
		stop('nproc must be specified on multi-user unix machines')

	if(nproc>detectCores())stop("The number of processes requested exceeds the number of processors available.")
	if(nproc>nsim){
	warning(paste('The number of processors requested exceeded the number of simulations\n', 
					'This iikely a mistake.  nproc set to 1'))
	nproc<-1
	}
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
						doselev=doselev,nparm=nparm,nparmA=nparmA,modType=modType,
						iparm=iparm,binary=binary,genObj=genObj,testMods=testMods,
						contMat=contMat,negEmax=negEmax,
						ed50cutoff=ed50cutoff,
						ed50lowcutoff=ed50lowcutoff,
						switchMod=switchMod,
						truncLambda=truncLambda,nfrac=nfrac,nfrac2=nfrac2,n=n)


	if(nproc==1){
		simout<-list(simrep(1,inlist))
	}else{
		cl<-makeCluster(nproc)
		registerDoParallel(cl)	
		simout<-foreach(j=1:nproc, .packages='clinDR') %dopar%{
			simrep(j,inlist)
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
	pVal<-rep(NA,nsim)
	selContrast<-rep(NA,nsim)
	colnames(predpop)<-doselev
	colnames(fitpredv)<-doselev
	colnames(sepredv)<-doselev
	colnames(sedifv)<-doselev
	colnames(mv)<-doselev
	
	init<-matrix( rep(NA,nsim*nparm),ncol=nparm )
	est3<-matrix( rep(NA,nsim*3),ncol=3 )
	est4<-matrix( rep(NA,nsim*4),ncol=4 )
	estA<-matrix( rep(NA,nsim*nparmA),ncol=nparmA )
	vc<-matrix(rep(NA,nsim*nparm^2),nrow=nsim)
	if(modType==3){
		colnames(init)<-c("led50","emax","e0")
		colnames(est3)<-c("led50","emax","e0")
	}else{
		colnames(init)<-c("led50","lambda","emax","e0")
		colnames(est4)<-c("led50","lambda","emax","e0")
	}
	if(!binary){
		sdv <- matrix(rep(NA, nsim * Ndose), ncol = Ndose)
	}else sdv<-NULL
	negC <- rep(F,nsim)
	bigC <- rep(F,nsim)
	
	fitType<-rep("",nsim)
	if(!binary)residSD<-numeric(nsim) else residSD<-NULL
	pop<-NULL
	popSD<-NULL
	
	for(j in 1:nproc){
		predpop[indmat[j,1]:indmat[j,2],]<-simout[[j]]$predpop	
		fitpredv[indmat[j,1]:indmat[j,2],]<-simout[[j]]$fitpredv
		sepredv[indmat[j,1]:indmat[j,2],]<-simout[[j]]$sepredv
		sedifv[indmat[j,1]:indmat[j,2],]<-simout[[j]]$sedifv
		mv[indmat[j,1]:indmat[j,2],]<-simout[[j]]$mv
		sdv[indmat[j,1]:indmat[j,2],]<-simout[[j]]$sdv
		pVal[indmat[j,1]:indmat[j,2]]<-simout[[j]]$pVal
		selContrast[indmat[j,1]:indmat[j,2]]<-simout[[j]]$selContrast
		init[indmat[j,1]:indmat[j,2],]<-simout[[j]]$init
		est3[indmat[j,1]:indmat[j,2],]<-simout[[j]]$est3
		est4[indmat[j,1]:indmat[j,2],]<-simout[[j]]$est4
		estA[indmat[j,1]:indmat[j,2],]<-simout[[j]]$estA
		vc[indmat[j,1]:indmat[j,2],]<-simout[[j]]$vc
		if(!binary)residSD[indmat[j,1]:indmat[j,2]]<-simout[[j]]$residSD
		negC[indmat[j,1]:indmat[j,2]]<-simout[[j]]$negC
		bigC[indmat[j,1]:indmat[j,2]]<-simout[[j]]$bigC
		fitType[indmat[j,1]:indmat[j,2]]<-simout[[j]]$fitType
		pop<-rbind(pop,simout[[j]]$pop)
		popSD<-c(popSD,simout[[j]]$popSD)
	}
	
	RNGkind(save.rng)
	if(exists('save.seed')).Random.seed<<-save.seed

  return(structure(list(description=description,
  											binary=binary,modType=modType,genObj=genObj, 
  											pop=pop,popSD=popSD,init=init,  
  											ed50contr=ed50contr, lambdacontr=lambdacontr,testMods=testMods,			
  											est4=est4,est3=est3,estA=estA, vc=vc,residSD=residSD,
  											fitType=fitType,pVal=pVal,selContrast=selContrast,
  											negEmax=negEmax,
  											ed50cutoff=ed50cutoff,ed50lowcutoff=ed50lowcutoff,switchMod=switchMod, 
  											negC=negC,
  											bigC=bigC,predpop=predpop,        
  											mv = mv, sdv = sdv, fitpredv = fitpredv,
  											sepredv = sepredv, sedifv = sedifv, 
  											rseed=rseed, idmax=idmax
											  ),class="emaxsim") )

}


simrep<-function(j,inlist)
	{

	indmat<-inlist$indmat
	rseed<-inlist$rseed
	Ndose<-inlist$Ndose
	dose<-inlist$dose
	doselev<-inlist$doselev
	nparm<-inlist$nparm
	nparmA<-inlist$nparmA
	modType<-inlist$modType
	iparm<-inlist$iparm
	binary<-inlist$binary
	genObj<-inlist$genObj
	testMods<-inlist$testMods
	contMat<-inlist$contMat
	negEmax<-inlist$negEmax
	ed50cutoff<-inlist$ed50cutoff	
	ed50lowcutoff<-inlist$ed50lowcutoff
	switchMod<-inlist$switchMod
	truncLambda<-inlist$truncLambda
	n<-inlist$n
	nfrac<-inlist$nfrac
	nfrac2<-inlist$nfrac2
	
 	nrep<-indmat[j,2]-indmat[j,1]+1 
 	
	predpop <- matrix(rep(NA, nrep * Ndose), ncol = Ndose)
	fitpredv <- matrix(rep(NA, nrep * Ndose), ncol = Ndose)
	sepredv <- matrix(rep(NA, nrep * Ndose), ncol = Ndose)
	sedifv <- matrix(rep(NA, nrep * Ndose), ncol = Ndose)
	mv <- matrix(rep(NA, nrep * Ndose), ncol = Ndose)
	pVal<-rep(NA,nrep)
	selContrast<-rep(NA,nrep)
	colnames(predpop)<-doselev
	colnames(fitpredv)<-doselev
	colnames(sepredv)<-doselev
	colnames(sedifv)<-doselev
	colnames(mv)<-doselev

	init<-matrix( rep(NA,nrep*nparm),ncol=nparm )
	est3<-matrix( rep(NA,nrep*3),ncol=3 )
	est4<-matrix( rep(NA,nrep*4),ncol=4 )
	estA<-matrix( rep(NA,nrep*nparmA),ncol=nparmA )
	vc<-matrix(rep(NA,nrep*nparm^2),nrow=nrep)
	if(modType==3){
		colnames(init)<-c("led50","emax","e0")
		colnames(est3)<-c("led50","emax","e0")
	}else{
		colnames(init)<-c("led50","lambda","emax","e0")
		colnames(est4)<-c("led50","lambda","emax","e0")
	}
	if(!binary){
		sdv <- matrix(rep(NA, nrep * Ndose), ncol = Ndose)
	}else sdv<-NULL
	negC <- rep(F,nrep)
	bigC <- rep(F,nrep)

	fitType<-rep("",nrep)
	fit<-vector('list',nrep)
	if(!binary)residSD<-numeric(nrep) else residSD<-NULL
	pop<-NULL
	popSD<-NULL

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
			lphat<-qlogis(phat)
			V<-diag(1/(phat*(1-phat)*n))
			contMat <-optContr(testMods,S=V)
		}

		### compute p-value for preliminary test
		if(binary){holdC<-MCTtest(doselev, lphat, contMat = contMat,
							S=V,df=Inf,type='general')$tStat
		}else holdC<-MCTtest(dose, y, contMat = contMat)$tStat

		holdP<-attr(holdC,'pVal')
		orderh<-order(holdC)
		ncontr<-length(orderh)
    pVal[i]<-holdP[orderh[ncontr]]
    selContrast[i]<-orderh[ncontr]
 	
		### main estimation code
		simout<-emaxalt(y,dose,modType,binary=binary,iparm,ed50cutoff,
						ed50lowcutoff,switchMod,truncLambda=truncLambda)

		### assign simulation output
		mv[i, ]  <- simout$dm
		if(!binary){
			sdv[i, ] <- simout$dsd
			residSD[i]<-simout$residSD
		}
		if(modType==3){init[i,1:3]<-simout$Sparm
		}else {init[i,]<-simout$Sparm}
		negC[i]<-simout$negC
		bigC[i]<-simout$bigC
		fitType[i]<-simout$fitType 
		if(fitType[i]=="4"){vc[i,]<-simout$vc
		}else if(fitType[i]=="3"){vc[i,1:9]<-simout$vc
		}else {vc[i,1:(nparmA^2)]<-simout$vc}
		fitpredv[i,]<-simout$fitpred
		sepredv[i,  ] <-simout$sepred
		sedifv[i,  ] <-simout$sedif
		est4[i,]<-simout$est4
		est3[i,]<-simout$est3
		estA[i,]<-simout$estA
	}
  
	return(list(pop=pop, popSD=popSD, init=init, 
				est4=est4,est3=est3,estA=estA, vc=vc,residSD=residSD,
				fitType=fitType,pVal=pVal,selContrast=selContrast,
				negC=negC,bigC=bigC,predpop=predpop,        
				mv = mv, sdv = sdv, fitpredv = fitpredv,
				sepredv = sepredv, sedifv = sedifv 
				))
}

