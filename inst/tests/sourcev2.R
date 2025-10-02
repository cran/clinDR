fitModel<-function(id,y,trt,visit,prmean0,prsd0,prmean,prsd,gparm1=3,gparm2=1.5,
									mcmc=mcmc.control()){
	## id must be 1,2,3...  without any skipped indices (a patient with
	##                    no observed resp must be deleted)
	## trt must be 1,2
	## visits must be numbered sequential 1,2,..., individual visits can be skipped
	##         but there must be at least 1 measurement for some patient at each visit
	## resp is 0/1
	
	## remove na
	indNA<-!is.na(y)
	id<-id[indNA]
	y<-y[indNA]
	trt<-trt[indNA]
	visit<-visit[indNA]
	
	### check and format inputs
	trtcheck<-sort(unique(trt))
	ntrt<-length(trtcheck)
	if(any(trtcheck!=1:ntrt))stop('trt must be sequentially numbered without skipping')
	
	if(!all(y %in% c(0,1)))stop('y must be 0/1')
	
	idcheck<-sort(unique(id))
	nsubj<-max(idcheck)
	if(any(idcheck!=1:nsubj))stop('id must be sequentially numbered without skipping')
	
	vcheck<-sort(unique(visit))
	nvisit<-max(vcheck)
	if(any(vcheck!=1:nvisit))stop('visits must be sequentially numbered')
	
	N<-length(id)
	ntrt1<-ntrt-1

	### stan fit
	indata<-c('N','nsubj','nvisit','ntrt','ntrt1','id','trt','visit','y',
						'prmean0','prsd0','prmean','prsd','gparm1','gparm2')
	
	parameters<-c('beta','sigma','theta')
	
	estanmod<-stan_model(file='imputeMIL2.stan',
											save_dso=TRUE,auto_write=FALSE,model_name='imputeMI')	
	stanfit<-sampling(estanmod,data=indata,
										chains=mcmc$chains,
										 warmup=mcmc$warmup,iter=mcmc$iter,thin=mcmc$thin,seed=mcmc$seed,
										 pars=parameters,
										 control=list(adapt_delta=mcmc$adapt_delta),
										 cores = mcmc$chains 
	)		
	#############################################################
	### convert generated parms to R matrices
	beta<-as.matrix(stanfit,pars='beta')
	beta<-array(as.vector(beta),dim=c(nrow(beta),max(visit),max(trt)))
	
	sigma<-as.vector(as.matrix(stanfit,pars='sigma'))
	theta<-as.matrix(stanfit,pars='theta')

	return(list(stanfit=stanfit,beta=beta,sigma=sigma,theta=theta))
}

inputmi<-function(id,trt,y,visit,trtsel=1:2,vsel=max(visit)){
	### data set with 2 selected trts and a single
	### visit on 1 record per patient is created
	### missing y MUST be represented by NA in the input data set
	
	indsub<- (trt%in%trtsel) & visit==vsel 
	id<-id[indsub]
	trt<-trt[indsub]
	y<-y[indsub]
	visit<-visit[indsub]
	m<-!is.na(y)
	return(list(id=id,trt=trt,y=y,m=m))
}

miprobs<-function(mdat,vsel,beta,sigma,nimp=100){
	## observed data and imputation probabilities conditional
	## on mcmc parameters for missing data
	if(nimp>dim(beta)[1])stop('too many imputations requested')
	nsubj<- length(mdat$id)

	probs<-matrix(numeric(nsubj*nimp),ncol=nimp)	
	probs[mdat$m==1,]<-mdat$y[mdat$m==1]
	
	missind<-which(mdat$m==0)	
	for(imp in c(1:nimp)){
		for(i in missind){
			probs[i,imp]<-plogis(beta[imp,vsel,mdat$trt[i]]+sigma[imp]*theta[imp,i])
		}
	}
	return(probs)
}

midat<-function(mprobs,trt,m,deltat=0,deltac=0,f=mnse, returnYimp=FALSE, ...){
	### multiple imputation results
	### mi probs output by miprobs
	### trt and m output by inputmi 
	### deltat, deltac for tipping point analyses
	###
	### default complete data SE using mnse with delmn=0
	### add delmn=xx to change the assumed dif in mnse
	### other se-computing functions can be supplied
	### they must have the same call and return forms
	### as mnse.  different additional parameters can
	### be supplied (in place of delmn)
	nsubj<-nrow(mprobs)		
	nimp<-ncol(mprobs)
	ptimp<-numeric(nimp)
	pcimp<-numeric(nimp)
	sedelimp<-numeric(nimp)
	yimp<-matrix(numeric(nsubj*nimp),ncol=nimp)
	
	yimp[m==1,]<-mprobs[m==1,1]  ## observed data
	missind<-which(m==0)	
	nmis<-length(missind)
	trtmis<-trt[m==0]
	if(nmis>0){
		doff<-rep(deltac,nmis)
		doff[trtmis==max(trt)]<-deltat
		for(imp in 1:nimp){
			genprob<-mprobs[missind,imp]-doff
			genprob[genprob<0]<-0
			genprob[genprob>1]<-1
			yimp[missind,imp]<-rbinom(nmis,1,genprob)
			impout<-f(yimp[,imp],trt, ...)		
			ptimp[imp]<-impout$phatt
			pcimp[imp]<-impout$phatc
			sedelimp[imp]<-impout$sedel
		}
	}else stop('no missing responses')
	
	difp<-ptimp-pcimp
	miest<-mean(difp)
	vb<-var(difp)          
	vw<-mean(sedelimp^2)
	mise<-sqrt( vw + (1+(1/nimp))*vb )  ### validated vs mitools package
	midf<-(nimp-1)*( 1 + vw/((1+1/nimp)*vb))^2
	if(!returnYimp){return(list(miest=miest,mise=mise,midf=midf))
	}else return(list(miest=miest,mise=mise,midf=midf,yimp=yimp))
}	
	
mnse<-function(y,trt,delmn=0){
	## mietenen se calculation
	nt<-sum(trt==2)
	nc<-sum(trt==1)
	phatt<-mean(y[trt==2])
	phatc<-mean(y[trt==1])
	
	theta<-nt/nc
	d<-phatc*delmn*(1-delmn)
	c<-delmn^2-delmn*(2*phatc+theta+1)+phatc+theta*phatt
	b<- -(1+theta+phatc+theta*phatt-delmn*(theta+2))
	a<-1+theta
	vu<-b^3/(3*a)^3 - b*c/(6*a^2) + d/(2*a)
	u<-sign(vu)*sqrt( b^2/(3*a)^2 - c/(3*a) )
	w<-(pi + acos(vu/u^3))/3
	
	pcr<-2*u*cos(w) - b/(3*a)
	ptr<- pcr+delmn
	
	sec<-sqrt( pcr*(1-pcr)/nc )
	setrt<-sqrt( ptr*(1-ptr)/nt )
	sedel<- sqrt( pcr*(1-pcr)/nc + ptr*(1-ptr)/nt )
	return(list(phatc=phatc,phatt=phatt,sedel=sedel))	
}


binNorm<-function(y,trt, ...){
	nt<-sum(trt==2)
	nc<-sum(trt==1)
	phatt<-mean(y[trt==2])
	phatc<-mean(y[trt==1])
	sedel<-sqrt( phatc*(1-phatc)/nc + phatt*(1-phatt)/nt )
return(list(phatc=phatc,phatt=phatt,sedel=sedel))	
}


