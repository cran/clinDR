"fitEmaxB"<-
function(y,dose,prior,modType=4,prot=rep(1,length(y)),count=rep(1,length(y)),
		 xbase=NULL,binary=FALSE,msSat=NULL,vcest=NULL,pboAdj=FALSE,
		 mcmc=mcmc.control(),estan=NULL,diagnostics=TRUE,
		 nproc=getOption("mc.cores", 1L)){

	if(nproc>parallel::detectCores())nproc<-parallel::detectCores()
	tol<-sqrt(.Machine$double.eps)
	protorig<-factor(prot)             #ensure that prot is a factor
	prot<-as.numeric(protorig)         #convert to 1,2,3,....
  nprot<-max(prot)
	dlev<-sort(unique(dose))
	sigmoid<-1*(modType==4)
	intercept<-1*(!pboAdj)
	
	if(!missing(vcest)){
		if(!is.matrix(vcest) || nrow(vcest)!=ncol(vcest) )
			stop('vcest must be a square matrix')
		dimFit<-nrow(vcest)	
		if(!missing(xbase))stop('only one of vcest and xbase can be specified')
		if(!isSymmetric(unname(vcest)))stop('vcest is not symmetric')
		npd<-any(eigen(vcest,symmetric=TRUE)$values<=0)
		if(npd)stop('vcest must be a positive-definite variance-covariance matrix')
		if(!missing(count))stop('only one of vcest and count can be specified')
		if(!missing(msSat))stop('only one of vcest and msSat can be specified')
	}else{
		dimFit<-0
		vcest<-matrix(1,nrow=1,ncol=1)
	}

	### check input consistency
	if(isTRUE(pboAdj && binary))stop('PBO adjustment not available with binary data')
	
	if(isTRUE( binary!=prior$binary ))stop('Binary specification in prior and model do not match')

	### check for missing data
	if(any(is.na(y)))stop('Missing y data should be imputed or deleted')

	lengthvec<-c(length(y),length(dose),length(prot),length(count))
	if(any(lengthvec)==0)stop('y, dose, prot, or count has length 0')
	if(any(abs(lengthvec-lengthvec[1])>tol))stop('Length of y,dose,prot, and count must be equal')
	
	nbase<-0
	if(!is.null(xbase)){
		if(is.vector(xbase)){
			nbase<-1
		}else if(!is.matrix(xbase)){
			stop('xbase must be a vector or matrix')
		}else nbase<-ncol(xbase)
	}
	
	if(nbase!=length(prior$basemu))stop('The covariate dimensions do not match the prior distribution')
	
	if(nbase>0){
		if(!all(abs(count-1)<tol))stop('Covariate adjustment cannot be requested when counts!=1')
		if(pboAdj)stop('Covariate adjustment cannot be combined with PBO adjustment')
		if(nbase>1 && !all.equal(dim(xbase),c(length(y),nbase))){
			stop('Dimensions of xbase are invalid')
		}
		if(nbase==1 && length(xbase)!=length(y))stop('Length of xbase does not match y')
		for(i in 1:nbase){
			xm<-tapply(xbase[,i],prot,mean)
			if(!all(abs(xm)<tol))stop('Xbase must be centered about the protocol means')
		}
	}
	

	patDat<-isTRUE(all(count==1) && !dimFit)
	

	### make local copy of control variables
	if(inherits(prior,'emaxPrior'))localParm<-TRUE else localParm<-FALSE
	if(!localParm && dimFit)stop('vcest available only with emaxPrior input')

	p50<-prior$p50
	
	if(localParm){
		epmu<-prior$epmu
		epsca<-prior$epsca
		difTargetmu<-prior$difTargetmu
		difTargetsca<-prior$difTargetsca
		dTarget<-prior$dTarget
		effDF<-prior$effDF
		if(length(effDF)==2){
			e0DF<-effDF[1]
			diftDF<-effDF[2]
		}else{
			e0DF<-effDF
			diftDF<-effDF			
		}
		parmDF<-prior$parmDF
		loged50mu<-prior$loged50mu
		loged50sca<-prior$loged50sca
		loglammu<-prior$loglammu
		loglamsca<-prior$loglamsca
		parmCor<-prior$parmCor
		### truncation of extreme ed50,lambda
		lowled50<-prior$lowled50
		highled50<-prior$highled50
		lowllam<-prior$lowllam
		highllam<-prior$highllam
	}else{
		epmu<-prior$epmu
		epsd<-prior$epsd
		emaxmu<-prior$emaxmu
		emaxsd<-prior$emaxsd
		edDF<-prior$edDF
		led50mu<-prior$led50mu
		led50sca<-prior$led50sca
		lama<-prior$lama
		lamb<-prior$lamb
		lamsca<-prior$lamsca
	}
	
	if(!binary && !dimFit){
		sigmalow<-prior$sigmalow
		sigmaup<-prior$sigmaup
	}else { # placeholder values ignored by stan but required for input
		sigmalow<-0.01
		sigmaup<- 1.0
	}		
	if(!is.null(prior$basemu)){
		basemu<-array(prior$basemu,dim=nbase)
		basevar<-prior$basevar
	}else{
		xbase<-matrix(1,nrow=0,ncol=0)
		basevar<-matrix(1,nrow=1,ncol=1)
		basemu<-array(1,dim=0)
	}
	
	chains<-mcmc$chains
	if(nproc>chains)nproc<-chains
	thin<-mcmc$thin
	warmup<-mcmc$warmup
	iter<-mcmc$iter
	propInit<-mcmc$propInit
	seed<-mcmc$seed
	adapt_delta<-mcmc$adapt_delta
	

	## locate compiled stan models
	if(!is.null(estan) && !inherits(estan,'stanmodel'))stop('estan model invalid')	

	if(isTRUE(grep("64",Sys.getenv("R_ARCH"))>0)){
		emod<-'comp64'
	}else emod<-'comp32'
	if(localParm)emod<-file.path(emod,'mrmodel.rds') else emod<-file.path(emod,'basemodel.rds')

	#############################################################################
	#### near-duplicate code for localParm and old code
	###  when deprecated, just remove everything within else block
	###
	if(localParm){
		### create initial values
		if(chains==1){ 
			e0init<-epmu
			difTargetinit<-difTargetmu
			loged50init<-loged50mu
			loglaminit<-loglammu
			if(nbase>0){
				binit<-matrix(numeric(nbase),nrow=nbase)
				binit[,1]<-basemu
			}
			if(!binary){
				mids<-(sigmaup+sigmalow)/2
				siginit<-mids 	
			}
		}else{
			low<-epmu-propInit*epsca
			high<-epmu+propInit*epsca	
			e0init<-seq(low,high,length=chains)		
			
			low<-difTargetmu-propInit*difTargetsca
			high<-difTargetmu+propInit*difTargetsca	
			difTargetinit<-seq(low,high,length=chains)				
			
			low<- -propInit*loged50sca
			high<- propInit*loged50sca
			loged50init<-seq(low,high,length=chains)				
			
			low<- -propInit*loglamsca
			high<- propInit*loglamsca
			loglaminit<-seq(low,high,length=chains)				
				
			if(nbase>0){
				blow<-basemu-propInit*sqrt(diag(basevar))
				bhigh<-basemu+propInit*sqrt(diag(basevar))
				binit<-matrix(numeric(nbase*chains),nrow=nbase)
				for(i in 1:nbase){
					binit[i,]<-seq(from=blow[i],to=bhigh[i],length=chains)
				}
			}
			if(!binary){
				mids<-(sigmaup+sigmalow)/2
				lows<-max(sigmalow,mids/2)
				highs<-mids
				siginit<-seq(lows,highs,length=chains)
			}
		}
		
		###############################################################################
		#### binary data for fitted model input 
		if(binary || dimFit){
			gp<-0  ## placeholder not used in stan
			cont<-0
			
			if(binary && !dimFit){
				### y must be coded 0/1
				if(any(y!=0 & y!=1))stop("y must be 0/1")
				
				### data summaries for input to stan
				
				protv<-prot
				if(!(nbase>0)){
					xind<-unique(cbind(dose,prot))		
					N<-nrow(xind)
					dv<-numeric(N)
					yvb<-numeric(N)
					nvb<-numeric(N)
					
					protv<-numeric(N)
					for(i in 1:N){
						nvb[i]<-sum(count[dose==xind[i,1] & prot==xind[i,2]])	
						yvb[i]<-sum(count[dose==xind[i,1] & prot==xind[i,2] & y==1])	
						dv[i]<-xind[i,1]
						protv[i]<-xind[i,2]
					}
				}else{
					N<-length(y); yvb<-y;  dv<-dose; nvb<-count 
				}
				yv<-numeric(N)  ## placeholder for stan	
				nv<-numeric(N)  ## 
			}else{  #fitted model
				N<-length(dose)
				### data summaries for input to stan
				yv<-y;  dv<-dose; nv<-count; protv<-prot
				### placeholder for binary
				yvb<-numeric(N); nvb<-numeric(N)
			}
			
			df2<-0;  ## placeholder for stan
			ssy<-0   ##

			inits<- vector("list", chains)
			if(nbase==0 && !pboAdj){
				if(modType==4){
					for(j in 1:chains){
						inits[[j]]<-
							list(e0=array(rep(e0init[j],nprot),dim=nprot),difTarget=difTargetinit[j],
									 loged50=loged50init[j],llam=loglaminit[j])				
					}	
					parameters<-c('led50','lambda','emax','e0','difTarget','loglambda')
				}else{
					for(j in 1:chains){
						inits[[j]]<-
							list(e0=array(rep(e0init[j],nprot),dim=nprot),difTarget=difTargetinit[j],
									 loged50=loged50init[j],llam=loglaminit[j])
					}	
					parameters<-c('led50','emax','e0','difTarget')
				}
			}else if(!pboAdj){
				if(modType==4){
					for(j in 1:chains){
						inits[[j]]<-
							list(e0=array(rep(e0init[j],nprot),dim=nprot),difTarget=difTargetinit[j],
									 loged50=loged50init[j],llam=loglaminit[j],
									 bslope=array(binit[,j],dim=nbase))				
					}	
					parameters<-c('led50','lambda','emax','e0','bslope','difTarget','loglambda')
				}else{
					for(j in 1:chains){
						inits[[j]]<-
							list(e0=array(rep(e0init[j],nprot),dim=nprot),difTarget=difTargetinit[j],
									 loged50=loged50init[j],llam=loglaminit[j],
									 bslope=array(binit[,j],dim=nbase))				
					}	
					parameters<-c('led50','emax','e0','bslope','difTarget')
				}			
			} else{  ### pbo excluded, applies with vcest only
				if(modType==4){
					for(j in 1:chains){
						inits[[j]]<-
							list(difTarget=difTargetinit[j],
									 loged50=loged50init[j],llam=loglaminit[j])				
					}					
					parameters<-c('led50','lambda','emax','difTarget','loglambda')
				}else{
					for(j in 1:chains){
						inits[[j]]<-
							list(difTarget=difTargetinit[j],
									 loged50=loged50init[j],llam=loglaminit[j])				
					}								
					parameters<-c('led50','emax','difTarget')
				}
			}			
			
			###############################################################################
			#### continuous data
		}else{
			cont<-1
			N<-length(dose)
			
			### data summaries for input to stan
			yv<-y;  dv<-dose; nv<-count; protv<-prot
			### placeholder for binary
			yvb<-numeric(N); nvb<-numeric(N)
			
			gp<-1-1*(isTRUE(patDat | is.null(msSat)) )
			if(isTRUE(is.null(msSat) && !patDat))warning('Grouped data and msSat not specified.  SD based on low DF')
			if(!gp){
				df2<-0
				ssy<-0
			}else{
				df2<-(sum(count)-N)/2
				ssy<-2*df2*msSat
			}
			
			inits<- vector("list", chains)
			if(!pboAdj){
				if(nbase==0){
					if(modType==4){
						for(j in 1:chains){
							inits[[j]]<-
								list(e0=array(rep(e0init[j],nprot),dim=nprot),difTarget=difTargetinit[j],
										 loged50=loged50init[j],llam=loglaminit[j],
										 sigma=array(siginit[j],dim=1))
						}	
						parameters<-c('led50','lambda','emax','e0','sigma','difTarget','loglambda')
					}else{
						for(j in 1:chains){
							inits[[j]]<-
								list(e0=array(rep(e0init[j],nprot),dim=nprot),max=difTargetinit[j],
										 loged50=loged50init[j],llam=loglaminit[j], 
										 sigma=array(siginit[j],dim=1) )				
						}	
						parameters<-c('led50','emax','e0','sigma','difTarget')
					}
				}else{
					if(modType==4){
						for(j in 1:chains){
							inits[[j]]<-
								list(e0=array(rep(e0init[j],nprot),dim=nprot),difTarget=difTargetinit[j],
										 loged50=loged50init[j],llam=loglaminit[j], 
										 bslope=array(binit[,j],dim=nbase),
										 sigma=array(siginit[j],dim=1))
						}	
						parameters<-c('led50','lambda','emax','e0','bslope','sigma','difTarget','loglambda')
					}else{
						for(j in 1:chains){
							inits[[j]]<-
								list(e0=array(rep(e0init[j],nprot),dim=nprot),difTarget=difTargetinit[j],
										 loged50=loged50init[j],llam=loglaminit[j], 
										 bslope=array(binit[,j],dim=nbase),
										 sigma=array(siginit[j],dim=1) )				
						}	
						parameters<-c('led50','emax','e0','bslope','sigma','difTarget')
					}			
				}
			}else{  ### pbo excluded
				if(modType==4){
					for(j in 1:chains){
						inits[[j]]<-
							list(difTarget=difTargetinit[j],
									loged50=loged50init[j],llam=loglaminit[j], 
									 sigma=array(siginit[j],dim=1))				
					}					
					parameters<-c('led50','lambda','emax','sigma','difTarget','loglambda')
				}else{
					for(j in 1:chains){
						inits[[j]]<-
							list(difTarget=difTargetinit[j],
									loged50=loged50init[j],llam=loglaminit[j], 
									 sigma=array(siginit[j],dim=1))				
					}								
					parameters<-c('led50','emax','sigma','difTarget')
				}
			}	
		}
		
		
		indata<-c('N','nprot','protv','cont','sigmoid','gp',
							'intercept','nbase',
							'yv','nv','yvb','nvb','dv','xbase','df2','ssy',
							'epmu','epsca','difTargetmu','difTargetsca','dTarget','sigmalow','sigmaup','p50',
							'loged50mu','loged50sca','e0DF','diftDF','parmDF','loglammu','loglamsca',
							'parmCor','basemu','basevar','lowled50','highled50','lowllam','highllam',
							'dimFit','vcest')		
	}else{
		### create initial values
		if(chains==1){ 
			e0init<-epmu
			emaxinit<-emaxmu
			led50init<-0
			if(modType==4)laminit<-1/lamsca
			if(nbase>0){
				binit<-matrix(numeric(nbase),nrow=nbase)
				binit[,1]<-basemu
			}
			if(!binary){
				mids<-(sigmaup+sigmalow)/2
				siginit<-mids 	
			}
		}else{
			low<-epmu-propInit*epsd
			high<-epmu+propInit*epsd	
			e0init<-seq(low,high,length=chains)		
			
			low<-emaxmu-propInit*emaxsd
			high<-emaxmu+propInit*emaxsd	
			emaxinit<-seq(low,high,length=chains)				
			
			low<- -led50mu
			high<-led50mu+propInit*led50sca
			led50init<-seq(low,high,length=chains)				
			
			if(modType==4){
				laminit<-seq(0.5/lamsca,2/lamsca,length=chains)
			}
			if(nbase>0){
				blow<-basemu-propInit*sqrt(diag(basevar))
				bhigh<-basemu+propInit*sqrt(diag(basevar))
				binit<-matrix(numeric(nbase*chains),nrow=nbase)
				for(i in 1:nbase){
					binit[i,]<-seq(from=blow[i],to=bhigh[i],length=chains)
				}
			}
			if(!binary){
				mids<-(sigmaup+sigmalow)/2
				lows<-max(sigmalow,mids/2)
				highs<-mids
				siginit<-seq(lows,highs,length=chains)
			}
		}
		
		###############################################################################
		#### binary data
		if(binary){
			gp<-0  ## placeholder not used in stan
			cont<-0
			### y must be coded 0/1
			if(any(y!=0 & y!=1))stop("y must be 0/1")
	
			
			### data summaries for input to stan
	
			protv<-prot
			if(!(nbase>0)){
				xind<-unique(cbind(dose,prot))		
				N<-nrow(xind)
				dv<-numeric(N)
				yvb<-numeric(N)
				nvb<-numeric(N)
	
				protv<-numeric(N)
				for(i in 1:N){
					nvb[i]<-sum(count[dose==xind[i,1] & prot==xind[i,2]])	
					yvb[i]<-sum(count[dose==xind[i,1] & prot==xind[i,2] & y==1])	
					dv[i]<-xind[i,1]
					protv[i]<-xind[i,2]
				}
			}else{
				N<-length(y); yvb<-y;  dv<-dose; nvb<-count 
			}
			df2<-0;  ## placeholder for stan
			ssy<-0   ##
			yv<-numeric(N)  ## 
			nv<-numeric(N)  ## placeholder for stan	
			
			inits<- vector("list", chains)
			if(nbase==0){
				if(modType==4){
					for(j in 1:chains){
					inits[[j]]<-
						list(e0=array(rep(e0init[j],nprot),dim=nprot),emax=emaxinit[j],
								 ed50t=led50init[j],lamt=array(laminit[j],dim=1))				
					}	
					parameters<-c('led50','lambda','emax','e0')
				}else{
					for(j in 1:chains){
					inits[[j]]<-
						list(e0=array(rep(e0init[j],nprot),dim=nprot),emax=emaxinit[j],
								 ed50t=led50init[j])				
					}	
					parameters<-c('led50','emax','e0')
				}
			}else{
				if(modType==4){
					for(j in 1:chains){
					inits[[j]]<-
						list(e0=array(rep(e0init[j],nprot),dim=nprot),emax=emaxinit[j],
								 ed50t=led50init[j],lamt=array(laminit[j],dim=1),
								 bslope=array(binit[,j],dim=nbase))				
					}	
					parameters<-c('led50','lambda','emax','e0','bslope')
				}else{
					for(j in 1:chains){
					inits[[j]]<-
						list(e0=array(rep(e0init[j],nprot),dim=nprot),emax=emaxinit[j],
								 ed50t=led50init[j],bslope=array(binit[,j],dim=nbase))				
					}	
					parameters<-c('led50','emax','e0','bslope')
				}			
			}
			
		###############################################################################
		#### continuous data
		}else{
			cont<-1
			N<-length(dose)
			
			### data summaries for input to stan
			yv<-y;  dv<-dose; nv<-count; protv<-prot
			### placeholder for binary
			yvb<-numeric(N); nvb<-numeric(N)
	
			gp<-1-1*(isTRUE(patDat | is.null(msSat)) )
			if(isTRUE(is.null(msSat) && !patDat))warning('Grouped data and msSat not specified.  SD based on low DF')
			if(!gp){
				df2<-0
				ssy<-0
			}else{
				df2<-(sum(count)-N)/2
				ssy<-2*df2*msSat
			}
			
			inits<- vector("list", chains)
			lamt<-numeric(1)
			sigma<-numeric(1)
			if(!pboAdj){
				if(nbase==0){
					if(modType==4){
						for(j in 1:chains){
							inits[[j]]<-
								list(e0=array(rep(e0init[j],nprot),dim=nprot),emax=emaxinit[j],
										 ed50t=led50init[j],lamt=array(laminit[j],dim=1),
										 sigma=array(siginit[j],dim=1))
						}	
						parameters<-c('led50','lambda','emax','e0','sigma')
					}else{
						for(j in 1:chains){
							inits[[j]]<-
								list(e0=array(rep(e0init[j],nprot),dim=nprot),emax=emaxinit[j],
										 ed50t=led50init[j], sigma=array(siginit[j],dim=1) )				
						}	
						parameters<-c('led50','emax','e0','sigma')
					}
				}else{
						if(modType==4){
						for(j in 1:chains){
							inits[[j]]<-
								list(e0=array(rep(e0init[j],nprot),dim=nprot),emax=emaxinit[j],
										 ed50t=led50init[j],lamt=array(laminit[j],dim=1),
										 bslope=array(binit[,j],dim=nbase),
										 sigma=array(siginit[j],dim=1))
						}	
						parameters<-c('led50','lambda','emax','e0','bslope','sigma')
					}else{
						for(j in 1:chains){
							inits[[j]]<-
								list(e0=array(rep(e0init[j],nprot),dim=nprot),emax=emaxinit[j],
										 ed50t=led50init[j],bslope=array(binit[,j],dim=nbase),
										 sigma=array(siginit[j],dim=1) )				
						}	
						parameters<-c('led50','emax','e0','bslope','sigma')
					}			
				}
			}else{  ### pbo excluded
					if(modType==4){
						for(j in 1:chains){
							inits[[j]]<-
								list(emax=emaxinit[j],
										 ed50t=led50init[j],lamt=array(laminit[j],dim=1),sigma=array(siginit[j],dim=1))				
					}					
					parameters<-c('led50','lambda','emax','sigma')
				}else{
					for(j in 1:chains){
						inits[[j]]<-
							list(emax=emaxinit[j],
									 ed50t=led50init[j],sigma=array(siginit[j],dim=1))				
					}								
					parameters<-c('led50','emax','sigma')
				}
			}	
		}
	
	
		indata<-c('N','nprot','protv','cont','sigmoid','gp',
							'intercept','nbase',
							'yv','nv','yvb','nvb','dv','xbase','df2','ssy',
					'epmu','epsd','emaxmu','emaxsd','sigmalow','sigmaup','p50',
						'led50mu','led50sca','edDF','lama','lamb','lamsca',
					'basemu','basevar','dimFit','vcest')		
	}
#######################################################################
	### assign compiled model file
	if(is.null(estan)){
		emod<-system.file(package="clinDR", "models", emod)
		if(file.access(emod,mode=0)<0)stop(paste('The compiled rstan model',
												'could not be accessed.  You must',
												'run compileStanModels once before using',
												'the Bayesian functions'))	
		estan<-readRDS(emod)	
	}
	
	nproc<-as.integer(nproc)  ## required by rstan
	if(diagnostics){
		estanfit<-sampling(estan,data=indata,
										 warmup=warmup,iter=iter,thin=thin,seed=seed,
										 chains=chains,init=inits,pars=parameters,
										 control=list(adapt_delta=adapt_delta),
										 						 cores = nproc 
										 )		
	}else{
		capture.output(
		estanfit<-suppressWarnings(sampling(estan,data=indata,
											 warmup=warmup,iter=iter,thin=thin,seed=seed,
											 chains=chains,init=inits,pars=parameters,
											 control=list(adapt_delta=adapt_delta),
											 verbose=FALSE,open_progress=FALSE,refresh=-1,
											 cores = nproc ))
		,file=NULL)
	}
	
	
	fit<-list(
		estanfit=estanfit,
		y=y, dose=dose, prot=protorig,count=count,
		nbase=nbase,xbase=xbase,dimFit=dimFit,vcest=vcest,
		modType=modType,binary=binary,pboAdj=pboAdj,
		msSat=msSat,prior=prior,mcmc=mcmc,localParm=localParm		
	)
	class(fit)<-'fitEmaxB'

	return(fit)
	
}

prior.control<-function(epmu=NULL,epsd=NULL,emaxmu=NULL,emaxsd=NULL,p50=NULL,
												sigmalow=NULL,sigmaup=NULL,
												led50mu=0.79,led50sca=0.60,
												edDF=3,lama=3.03,lamb=18.15,
												lamsca=6,basemu=NULL,basevar=NULL,
												binary=FALSE)
{
	### defaults based on meta-analyses of dose response
	
	.Deprecated("emaxPrior",package="clinDR",
		msg=paste('This prior distribution and its implementation ',
							'in functions fitEmaxB/emaxsimB will be eliminated ',
							'in a future version of package clinDR. Specify the ',
							'prior/parameterization using emaxPrior.Control',sep=''))
	
	if(!binary){
		if( any( is.null(epmu),is.null(epsd),is.null(emaxmu),is.null(emaxsd),
					is.null(p50),is.null(sigmalow),is.null(sigmaup) ) 
			)stop('All parameters without default values must be specified')
	}else{ 
		if( any( is.null(epmu),is.null(epsd),is.null(emaxmu),is.null(emaxsd),
				is.null(p50) )
			)stop('All parameters without default values must be specified')
	}
	
	if(!is.null(basemu) | !is.null(basevar)){
		nbase<-length(basemu)
		if(!is.matrix(basevar))stop('basevar must be an nbbasexnbase variance-covariance matrix')
		if(any(dim(basevar)!=nbase))stop('basevar must be an nbbase x nbase variance-covariance matrix')
		if(!isSymmetric(unname(basevar)))stop('basevar is not symmetric')
		if(any(eigen(basevar,symmetric=TRUE)$values<=0))stop('basevar is not positive-definite')
	}

	list(epmu=epmu,
			 epsd=epsd,
			 emaxmu=emaxmu,
			 emaxsd=emaxsd,
			 p50=p50,
			 sigmalow=sigmalow,
			 sigmaup=sigmaup,
			 led50mu=led50mu,
			 led50sca=led50sca,
			 edDF=edDF,
			 lama=lama,
			 lamb=lamb,
			 lamsca=lamsca,
			 basemu=basemu,
			 basevar=basevar,
			 binary=binary)
}

emaxPrior.control<-function(epmu=NULL,epsca=NULL,
														difTargetmu=NULL,difTargetsca=NULL,dTarget=NULL,
														p50=NULL,sigmalow=NULL,sigmaup=NULL,
														effDF=parmDF,parmDF=5,
														loged50mu=0.0,loged50sca=1.73,
														loglammu=0.0,loglamsca=0.425,parmCor=-0.45,
														lowled50=log(0.001),highled50=log(1000),
														lowllam=log(0.3),highllam=log(4.0),
														basemu=NULL,basevar=NULL,
														binary=FALSE)
{
	### defaults based on meta-analyses of dose response
	
	if(!binary){
		if( any( is.null(epmu),is.null(epsca),is.null(difTargetmu),is.null(difTargetsca),
						 is.null(p50),is.null(sigmalow),is.null(sigmaup) ) 
		)stop('All parameters without default values must be specified')
	}else{ 
		if( any( is.null(epmu),is.null(epsca),is.null(difTargetmu),is.null(difTargetsca),
						 is.null(p50) )
		)stop('All parameters without default values must be specified')
	}
	
	if(!is.null(basemu) | !is.null(basevar)){
		nbase<-length(basemu)
		if(!is.matrix(basevar))stop('basevar must be an nbbasexnbase variance-covariance matrix')
		if(any(dim(basevar)!=nbase))stop('basevar must be an nbbase x nbase variance-covariance matrix')
		if(!isSymmetric(unname(basevar)))stop('basevar is not symmetric')
		if(any(eigen(basevar,symmetric=TRUE)$values<=0))stop('basevar is not positive-definite')
	}
	if(is.null(dTarget))stop('dTarget must be specified')
	if(!missing(effDF) && (length(effDF)>2))stop('effDF must have length <=2')
	
	default<-FALSE
	if(missing(loged50mu) && missing(loged50sca) &&
		 missing(loglammu) && missing(loglamsca))default<-TRUE 
	
	prior<-list(epmu=epmu,
							epsca=epsca,
							difTargetmu=difTargetmu,
							difTargetsca=difTargetsca,
							dTarget=dTarget,
							p50=p50,
							sigmalow=sigmalow,
							sigmaup=sigmaup,
							effDF=effDF,
							parmDF=parmDF,
							loged50mu=loged50mu,
							loged50sca=loged50sca,
							loglammu=loglammu,
							loglamsca=loglamsca,
							parmCor=parmCor,
							lowled50=lowled50,
							highled50=highled50,
							lowllam=lowllam,
							highllam=highllam,
							basemu=basemu,
							basevar=basevar,
							binary=binary,
							default=default)
	class(prior)<-'emaxPrior'
	return(prior)
}

mcmc.control<-function(chains=1,thin=1,warmup=1000,iter=3333*thin+warmup,
											 propInit=0.50,seed=12357,
											 adapt_delta=0.90)
{
	list(chains=chains,thin=thin,warmup=warmup,iter=iter,propInit=propInit,
			 seed=seed,adapt_delta=adapt_delta)	
}


