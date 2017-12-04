"fitEmaxB"<-
function(y,dose,prior,modType=3,prot=rep(1,length(y)),count=rep(1,length(y)),
		 binary=FALSE,msSat=NULL,pboAdj=FALSE,
		 mcmc=mcmc.control(),estan=NULL,diagnostics=TRUE){

	tol<-.Machine$double.eps
	protorig<-factor(prot)             #ensure that prot is a factor
	prot<-as.numeric(protorig)         #convert to 1,2,3,....
  nprot<-max(prot)
	dlev<-sort(unique(dose))

	### check input consistency
	if(isTRUE(pboAdj && binary))stop('PBO adjustment not available with binary data')
	
	if(isTRUE( binary==prior$binary ))stop('binary specification in prior and model do not match')

	### check for missing data
	if(any(is.na(y)))stop('Missing y data should be imputed or deleted')

	lengthvec<-c(length(y),length(dose),length(prot),length(count))
	if(any(abs(lengthvec-lengthvec[1])>tol))stop('Length of y,dose,prot, and count must be equal')

	patDat<-isTRUE(all(count==1))
	
	if(!is.null(estan) && class(estan)!='stanmodel')stop('estan model invalid')	
	

	if(isTRUE(grep("64",Sys.getenv("R_ARCH"))>0)){
		emod<-'comp64'
	}else emod<-'comp32'

	
	### make local copy of control variables
	epmu<-prior$epmu
	epsd<-prior$epsd
	emaxmu<-prior$emaxmu
	emaxsd<-prior$emaxsd

	if(!binary){
		sigmalow<-prior$sigmalow
		sigmaup<-prior$sigmaup
	}
	p50<-prior$p50
	led50mu<-prior$led50mu
	led50sca<-prior$led50sca
	edDF<-prior$edDF
	lama<-prior$lama
	lamb<-prior$lamb
	lamsca<-prior$lamsca
	
	chains<-mcmc$chains
	thin<-mcmc$thin
	warmup<-mcmc$warmup
	iter<-mcmc$iter
	propInit<-mcmc$propInit
	seed<-mcmc$seed
	adapt_delta<-mcmc$adapt_delta
	
	### create initial values
	if(chains==1){ 
		e0init<-epmu
		emaxinit<-emaxmu
		led50init<-0
		if(modType==4)laminit<-1/lamsca
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
	}
	
	if(!binary){
		mids<-(sigmaup+sigmalow)/2
		if(chains==1){
			siginit<-mids
		}else{
			lows<-max(sigmalow,mids/2)
			highs<-mids
			siginit<-seq(lows,highs,length=chains)
		}
	}

	###############################################################################
	#### binary data
	if(binary){
		### y must be coded 0/1
		if(any(y!=0 & y!=1))stop("y must be 0/1")

		emod<-file.path(emod,'estanbin')
		
		### data summaries for input to stan
		xind<-unique(cbind(dose,prot))		
		N<-nrow(xind)
		yv<-numeric(N)
		dv<-numeric(N)
		nv<-numeric(N)
		protv<-numeric(N)
		for(i in 1:N){
			nv[i]<-sum(count[dose==xind[i,1] & prot==xind[i,2]])	
			yv[i]<-sum(count[dose==xind[i,1] & prot==xind[i,2] & y==1])	
			dv[i]<-xind[i,1]
			protv[i]<-xind[i,2]
		}

		indata<-c('N','nprot','protv','yv','nv','dv','epmu','epsd','emaxmu','emaxsd','p50',
							'led50mu','led50sca','edDF','lama','lamb','lamsca')	
		
		inits<- vector("list", chains)
		if(modType==4){
			for(j in 1:chains){
			inits[[j]]<-
				list(e0=array(rep(e0init[j],nprot),dim=nprot),emax=emaxinit[j],
						 ed50t=led50init[j],lamt=laminit[j])				
			}	
			parameters<-c('led50','lambda','emax','e0')
		}else{
			inits<- vector("list", chains)
			for(j in 1:chains){
			inits[[j]]<-
				list(e0=array(rep(e0init[j],nprot),dim=nprot),emax=emaxinit[j],
						 ed50t=led50init[j])				
			}	
			parameters<-c('led50','emax','e0')
		}
		emod<-paste(emod,modType,sep='')
		
	###############################################################################
	#### continuous data
	}else{

		emod<-file.path(emod,'estancont')
		
		### data summaries for input to stan
		yv<-y;  dv<-dose;
		N<-length(yv)
		gp<-1-1*(isTRUE(patDat | is.null(msSat)) )
		if(isTRUE(is.null(msSat) && !patDat))warning('Grouped data and msSat not specified.  SD based on low DF')
		if(!gp){
			df2<-0
			ssy<-0
		}else{
			df2<-(sum(count)-N)/2
			ssy<-2*df2*msSat
		}
		indata<-c('N','nprot','prot','yv','count','dv','gp','df2','ssy','epmu','epsd','emaxmu','emaxsd','p50',
							'sigmalow','sigmaup',
							'led50mu','led50sca','edDF','lama','lamb','lamsca')	
		
		inits<- vector("list", chains)
		if(!pboAdj){
			if(modType==4){
				for(j in 1:chains){
					inits[[j]]<-
						list(e0=array(rep(e0init[j],nprot),dim=nprot),emax=emaxinit[j],
								 ed50t=led50init[j],lamt=laminit[j],sigma=siginit[j])				
				}	
				parameters<-c('led50','lambda','emax','e0','sigma')
			}else{
				for(j in 1:chains){
					inits[[j]]<-
						list(e0=array(rep(e0init[j],nprot),dim=nprot),emax=emaxinit[j],
								 ed50t=led50init[j],sigma=siginit[j])				
				}	
				parameters<-c('led50','emax','e0','sigma')
			}
			emod<-paste(emod,modType,sep='')
		}else{  ### pbo excluded
				if(modType==4){
					for(j in 1:chains){
						inits[[j]]<-
							list(emax=emaxinit[j],
									 ed50t=led50init[j],lamt=laminit[j],sigma=siginit[j])				
				}					
				parameters<-c('led50','lambda','emax','sigma')
			}else{
				for(j in 1:chains){
					inits[[j]]<-
						list(emax=emaxinit[j],
								 ed50t=led50init[j],sigma=siginit[j])				
				}								
				parameters<-c('led50','emax','sigma')
			}
			emod<-paste(emod,modType,'NoPbo',sep='')		
		}	
	
	


	}
	
	emod<-paste(emod,'.rds',sep='')
	
	### assign compiled model file
	if(is.null(estan)){
		emod<-system.file(package="clinDR", "models", emod)
		if(file.access(emod,mode=0)<0)stop(paste('The compiled rstan model',
												'could not be accessed.  You must',
												'run compileStanModels once before using',
												'the Bayesian functions'))	
		estan<-readRDS(emod)	
	}
	
	if(diagnostics){
		estanfit<-sampling(estan,data=indata,
										 warmup=warmup,iter=iter,thin=thin,seed=seed,
										 chains=chains,init=inits,pars=parameters,
										 control=list(adapt_delta=adapt_delta)
										 )		
	}else{
		capture.output(
		estanfit<-sampling(estan,data=indata,
											 warmup=warmup,iter=iter,thin=thin,seed=seed,
											 chains=chains,init=inits,pars=parameters,
											 control=list(adapt_delta=adapt_delta),
											 verbose=FALSE,open_progress=FALSE,refresh=-1
											 )
		,file=NULL)
	}
	
	
	fit<-list(
		estanfit=estanfit,
		y=y, dose=dose, prot=protorig,count=count,
		modType=modType,binary=binary,pboAdj=pboAdj,msSat=msSat,
		prior=prior,mcmc=mcmc		
	)
	class(fit)<-'fitEmaxB'

	return(fit)
	
}

prior.control<-function(epmu=NULL,epsd=NULL,emaxmu=NULL,emaxsd=NULL,p50=NULL,
												sigmalow=NULL,sigmaup=NULL,
												led50mu=0.79,led50sca=0.60,
												edDF=3,lama=3.03,lamb=18.15,
												lamsca=6,binary=FALSE)
{
	### defaults based on meta-analyses of dose response
	
	if(!binary){
		if( any( is.null(epmu),is.null(epsd),is.null(emaxmu),is.null(emaxsd),
					is.null(p50),is.null(sigmalow),is.null(sigmaup) ) 
			)stop('All parameters without default values must be specified')
	}else{ 
		if( any( is.null(epmu),is.null(epsd),is.null(emaxmu),is.null(emaxsd),
				is.null(p50) )
			)stop('All parameters without default values must be specified')
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
			 lamsca=lamsca)
}
mcmc.control<-function(chains=1,thin=1,warmup=500,iter=5000*thin,
											 propInit=0.25,seed=12357,
											 adapt_delta=0.8)
{
	list(chains=chains,thin=thin,warmup=warmup,iter=iter,propInit=propInit,
			 seed=seed,adapt_delta=adapt_delta)	
}
	