context('fiEmaxB with updated prior')

#### NOTE: processors set to 15 for final simulation

### stan parallel options
options(mc.cores = parallel::detectCores())


### continuous
set.seed(12357)

doselev<-c(0,5,25,50,100,350)
n<-c(78,81,81,81,77,80)
n1<-sum(n)
n2<-sum(n[1:4])

doselev<-c(doselev,doselev[1:4])
n<-c(n,n[1:4])

### population parameters for simulation
e0<-2.465375 
ed50<-67.481113 
emax<-15.127726
sdy<-8.0
pop<-c(log(ed50),emax,e0)    
dose<-rep(doselev,n)
meanlev<-emaxfun(dose,pop)  


y<-rnorm(n1+n2,meanlev,sdy)
prots<-c(rep(1,n1),rep(2,n2))

prior<-emaxPrior.control(0,30,0,30,dTarget=350,p50=50,sigmalow=.1,sigmaup=30,parmDF=5)
mcmc<-mcmc.control(chains=3,warmup=500,iter=3000,seed=53453,propInit=0.15,adapt_delta = .9)

testout<-suppressWarnings(fitEmaxB(y,dose,prior=prior,modType=4,prot=prots,
									mcmc=mcmc,diagnostics=FALSE,nproc=3))

parms<-coef(testout)
estimate<-apply(parms,2,mean)
se<-sqrt(diag(var(parms)))
z<-(estimate-c(pop[1],1,pop[2:3],pop[3]))/se

### check parameter estimates
test_that("model parameters agree within 2.5se",{
	expect_lt(as.numeric(max(abs(z))),2.5)
})

### check predictions
predout<-predict(testout,dosevec=c(20,80),int=2)

poppred<-emaxfun(c(20,80),pop[c(1:3)])
z<-(predout$pred-poppred)/predout$se
zdif<-(predout$fitdif-(poppred-e0))/predout$sedif


### check predictions 
test_that("predictions agree within 2.5se",{
	expect_lt(as.numeric(max(abs(z))),2.5)
	expect_lt(as.numeric(max(abs(zdif))),2.5)
})

test_that("check absolute levels",{
	expect_that(as.numeric(predout$pred),
							equals(poppred,tol=2*sdy/sqrt(70),scale=1))
	expect_that(as.numeric(predout$fitdif),
							equals((poppred-e0),tol=2*sdy/sqrt(70),scale=1))
})

##########################################################
### check aggregated, with larger n for better asymptotics

set.seed(12357)

doselev<-c(0,5,25,50,100,350)
n<-5*c(78,81,81,81,77,80)
n1<-sum(n)
n2<-sum(n[1:4])

doselev<-c(doselev,doselev[1:4])
n<-c(n,n[1:4])

### population parameters for simulation
e0<-2.465375 
ed50<-67.481113 
emax<-15.127726
sdy<-8.0
pop<-c(log(ed50),emax,e0)    
dose<-rep(doselev,n)
meanlev<-emaxfun(dose,pop)  


y<-rnorm(n1+n2,meanlev,sdy)
prots<-c(rep(1,n1),rep(2,n2))

ymean<-tapply(y,list(dose,prots),mean)
ymean<-ymean[!is.na(ymean)]
msSat<-(summary(lm(y~factor(dose)))$sigma)^2
protshort<-c(rep(1,6),rep(2,4))
nag<-table(dose,prots)
nag<-as.vector(nag[nag>0])


prior<-emaxPrior.control(0,30,0,30,350,50,0.1,30,parmDF=5)

mcmc<-mcmc.control(chains=3,warmup=500,iter=3000,seed=53453,propInit=0.15,adapt_delta = .9)

suppressWarnings(testout<-fitEmaxB(ymean,doselev,prior=prior,modType=4,
									prot=protshort,count=nag,msSat=msSat,
									mcmc=mcmc,diagnostics=FALSE,nproc=3))

parms<-coef(testout)
estimate<-apply(parms,2,mean)
se<-sqrt(diag(var(parms)))
z<-(estimate-c(pop[1],1,pop[2:3],pop[3]))/se

### check parameter estimates
test_that("model parameters agree within 2.5se",{
	expect_lt(as.numeric(max(abs(z))),2.5)
})

### check predictions
predout<-predict(testout,dosevec=c(20,80),int=2)

poppred<-emaxfun(c(20,80),pop[c(1:3)])
z<-(predout$pred-poppred)/predout$se
zdif<-(predout$fitdif-(poppred-e0))/predout$sedif


### check predictions 
test_that("predictions agree within 2.5se",{
	expect_lt(as.numeric(max(abs(z))),2.5)
	expect_lt(as.numeric(max(abs(zdif))),2.5)
})

test_that("check absolute levels",{
	expect_that(as.numeric(predout$pred),
							equals(poppred,tol=2*sdy/sqrt(70),scale=1))
	expect_that(as.numeric(predout$fitdif),
							equals((poppred-e0),tol=2*sdy/sqrt(70),scale=1))
})

###############################################################
###############################################################
#### repeat with 3 parm 
#### no intercept
####
set.seed(12357)

doselev<-c(0,5,25,50,100,350)
n<-10*c(78,81,81,81,77,80)
n1<-sum(n)
n2<-sum(n[1:4])

doselev<-c(doselev,doselev[1:4])
n<-c(n,n[1:4])

### population parameters for simulation
e0<-0
ed50<-67.481113 
emax<-15.127726
sdy<-8.0
pop<-c(log(ed50),emax,e0)    
dose<-rep(doselev,n)
meanlev<-emaxfun(dose,pop)  


y<-rnorm(n1+n2,meanlev,sdy)
prots<-c(rep(1,n1),rep(2,n2))

ysub<-y[dose!=0]
dsub<-dose[dose!=0]
protsub<-prots[dose!=0]

prior<-emaxPrior.control(0,30,0,30,350,50,0.1,30,parmDF=5)
mcmc<-mcmc.control(chains=3,warmup=500,iter=3000,seed=53453,propInit=0.15,adapt_delta = .9)


suppressWarnings(testout<-fitEmaxB(ysub,dsub,prior,modType=3,prot=protsub,pboAdj=TRUE,
									mcmc=mcmc,diagnostics=FALSE,nproc=3))

parms<-coef(testout)
estimate<-apply(parms,2,mean)
se<-sqrt(diag(var(parms)))
z<-(estimate-pop[1:2])/se

### check parameter estimates
test_that("pboadj model parameters agree within 2se",{
	expect_lt(as.numeric(max(abs(z))),2.0)
})

### check predictions
predout<-predict(testout,dosevec=c(20,80),int=2)

poppred<-emaxfun(c(20,80),pop)-emaxfun(0,pop)

z<-(predout$pred-poppred)/predout$se
zdif<-(predout$fitdif-(poppred-e0))/predout$sedif


test_that("pboadj predictions agree within 2.5se",{
	expect_lt(as.numeric(max(abs(z))),2.5)
	expect_lt(as.numeric(max(abs(zdif))),2.5)
})

test_that("pboadj check absolute levels",{
	expect_that(as.numeric(predout$pred),
							equals(poppred,tol=2*sdy/sqrt(700),scale=1))
	expect_that(as.numeric(predout$fitdif),
							equals((poppred-e0),tol=2*sdy/sqrt(70),scale=1))
})


###############################################################
#### no intercept with 4 parm, non-zero reference group
#### large sample sizes, grouped

set.seed(12357)

doselev<-c(0,5,25,50,100,350)
n<-100*c(78,81,81,81,77,80)
n1<-sum(n)
n2<-sum(n[1:4])

doselev<-c(doselev,doselev[1:4])
n<-c(n,n[1:4])

### population parameters for simulation
e0<-0
ed50<-67.481113 
emax<-15.127726
sdy<-8.0
pop<-c(log(ed50),emax,e0)    
dose<-rep(doselev,n)
meanlev<-emaxfun(dose,pop)  


y<-rnorm(n1+n2,meanlev,sdy)
prots<-c(rep(1,n1),rep(2,n2))

ysub<-y[dose!=0]
dsub<-dose[dose!=0]
protsub<-prots[dose!=0]

prior<-emaxPrior.control(0,30,0,30,350,50,0.1,30,parmDF=5)
mcmc<-mcmc.control(chains=3,warmup=500,iter=3000,seed=53453,propInit=0.15,adapt_delta = .9)



ymean<-tapply(ysub,list(dsub,protsub),mean)
ymean<-ymean[!is.na(ymean)]
msSat<-(summary(lm(y~factor(dose)))$sigma)^2
protshort<-c(rep(1,5),rep(2,3))
nag<-table(dsub,protsub)
nag<-as.vector(nag[nag>0])
dlevsub<-doselev[doselev!=0]

suppressWarnings(testout<-fitEmaxB(ymean,dlevsub,prior,modType=4,prot=protshort,
									count=nag,pboAdj=TRUE,msSat=msSat,
									mcmc=mcmc,diagnostics=FALSE,nproc=3))

parms<-coef(testout)
estimate<-apply(parms,2,mean)
se<-sqrt(diag(var(parms)))
z<-(estimate-c(pop[1],1,pop[2]))/se

### check parameter estimates
test_that("pboadj4 model parameters agree within 2se",{
	expect_lt(as.numeric(max(abs(z))),2.0)
})

### check predictions
predout<-predict(testout,dosevec=c(20,80),int=1,dref=50)

poppred<-emaxfun(c(20,80),c(pop[1],1,pop[2],0))
popref<-emaxfun(50,c(pop[1],1,pop[2],0))

z<-(predout$pred-poppred)/predout$se
zdif<-(predout$fitdif-(poppred-popref))/predout$sedif

### check predictions 
test_that("pboadj4 predictions agree within 2.5se",{
	expect_lt(as.numeric(max(abs(z))),2.5)
	expect_lt(as.numeric(max(abs(zdif))),2.5)
})

test_that("pboadj4 check absolute levels",{
	expect_that(as.numeric(predout$pred),
							equals(poppred,tol=2*sdy/sqrt(700),scale=1))
	expect_that(as.numeric(predout$fitdif),
							equals((poppred-popref),tol=2*sdy/sqrt(70),scale=1))
})


#############################################################
#############################################################
### 4 parm model and grouped data, one protocol
### and replicated measurement per dose/protocol condition
runsim<-function(j,seed,nsim){
	set.seed(seed[j])
	doselev<-c(0,5,25,50,100,350)
	n<-50*c(78,81,81,81,77,80)
	n1<-sum(n)
	n2<-sum(n[1:4])
	
	doselev<-c(doselev,doselev[1:4])
	n<-c(n,n[1:4])
	
	### population parameters for simulation
	e0<-2.465375 
	ed50<-67.481113 
	emax<-15.127726
	sdy<-1
	pop<-c(log(ed50),emax,e0)    
	dose<-rep(doselev,n)
	meanlev<-emaxfun(dose,pop)  
	poppred<-emaxfun(c(20,80),pop)
	popref<-emaxfun(50,pop)
	
	
	modtype<-4
	if(modtype==4){pparm<-c(pop[1],1,pop[2:3])
	}else pparm<-pop
	z<-matrix(numeric(modtype*nsim),ncol=modtype)
	zabs<-matrix(numeric(nsim*2),ncol=2)
	zdif<-matrix(numeric(nsim*2),ncol=2)
	prior<-emaxPrior.control(0,30,0,30,350,50,0.1,30,parmDF=5)
	mcmc<-mcmc.control(chains=1,warmup=500,iter=5000,seed=53453,propInit=0.15,adapt_delta = .9)
	estan<-selEstan('mrmodel')
	
	for(i in 1:nsim){
		y<-rnorm(n1+n2,meanlev,sdy)
		prots<-c(rep(1,n1),rep(2,n2))
		### by dose
		ysum<-tapply(y,dose,mean)
		nsum<-as.numeric(table(dose))
		msSat<-tapply(y,dose,var)
		msSat<-sum((nsum-1)*msSat)/(sum(nsum)-length(nsum))
		### by dose/prots
		ysum<-c(tapply(y[prots==1],dose[prots==1],mean),
						tapply(y[prots==2],dose[prots==2],mean))
		
		suppressWarnings(testout<-fitEmaxB(ysum,doselev,prior=prior,count=n,
											 modType=modtype,
											msSat=msSat,mcmc=mcmc,estan=estan,
											diagnostics=FALSE,nproc = 1))
		
		parms<-coef(testout)
		estimate<-apply(parms,2,mean)
		se<-sqrt(diag(var(parms)))
		z[i,]<-(estimate-pparm)/se
		predout<-predict(testout,dosevec=c(20,80),int=1,dref=50)
		zabs[i,]<-(predout$pred-poppred)/predout$se
		zdif[i,]<-(predout$fitdif-(poppred-popref))/predout$sedif
	}
	return(list(z=z,zabs=zabs,zdif=zdif))
}

nsim<-67
### set up independent stream of random numbers for
### each simulation iteration.
RNGkind("L'Ecuyer-CMRG")
set.seed(12357)
seed<-matrix(integer(nsim*7),ncol=7)
seed[1,]<-as.integer(.Random.seed)
for(i in 2:nsim){
 seed[i,]<-nextRNGStream(seed[i-1,])
}
 
cl<-makeCluster(nprocdef)
registerDoParallel(cl)	
outsim<-foreach(j=1:nprocdef, .packages=c('clinDR')) %dopar%{
	runsim(j,seed,nsim)
}
stopCluster(cl)
RNGkind("default")

z<-NULL
zabs<-NULL
zdif<-NULL
for(i in 1:nprocdef){
	z<-rbind(z,outsim[[i]]$z)
	zabs<-rbind(zabs,outsim[[i]]$zabs)
	zdif<-rbind(zdif,outsim[[i]]$zdif)
}

nsim<-nsim*nprocdef
### check parameter estimates
test_that("grouped data model parameters agree within 3se",{
	expect_lt(as.numeric(mean(apply(abs(z)>2,2,mean))),0.075)
})

test_that("predictions agree within 3se",{
	expect_that(0.05,
							equals(as.numeric(mean(apply(abs(zabs)>1.96,2,mean))),
										 tolerance=3*sqrt(.05*.95/nsim),scale=1))
})

test_that("dif predictions agree within 2se",{
	expect_lt( as.numeric(max(apply(abs(zdif)>1.96,2,mean))),0.075 )
})

#########################################################
#########################################################
#### continuous with vcest specified
####

#################
### base case with missing values
### fully saturated first-stage
#############
###
set.seed(12357)
nrep<-10
nd<-6
nv<-4
sig<-sqrt(10/4)

doselev<-c(0,1,2,4,8,16)
dose<-sort(rep(doselev,nrep*nv))

id<-sort(rep(1:(nrep*nd),nv))
vis<-rep(1:nv,nd*nrep)

led50<-log(3)
emax<-5
e0<-0
popparm<-c(led50,emax,e0)
misprop<-0.15  # mcar

modmean<-emaxfun(dose,popparm)+vis
poppred<-tapply(modmean,list(dose,vis),mean)[,nv]
poppred<-tapply(modmean,list(dose,vis),mean)[,nv]
popcov<-matrix(rep(0.25,nv^2),ncol=nv)
popcov[1,2]<-0.75
popcov[2,1]<-0.75
diag(popcov)<-1
popcov<-popcov*sig^2

y<-modmean+as.vector(t(rmvnorm(nd*nrep,rep(0,nv),popcov)))
dosefac<-factor(dose)
visfac<-factor(vis)

### impose missing data
ntot<-nd*nrep*nv
misid<-sample(1:ntot,round(ntot*misprop),replace = FALSE)
y<-y[-misid]
dosefac<-dosefac[-misid]
visfac<-visfac[-misid]
vis<-vis[-misid]
id<-id[-misid]

modfit<-gls(y ~ dosefac*visfac-1, 
						correlation = corSymm(form = ~ vis | id),
						weights = varIdent(form = ~ 1 | vis))

preddat<-data.frame(dosefac=factor(doselev),
										visfac=factor(rep(nv,nd),levels=c(1:4)))
predvals<-predict(modfit,preddat)

L<-model.matrix(~ dosefac*visfac-1,preddat)
vcpred<-L%*%tcrossprod(vcov(modfit),L)

###############################################
## bayes model fit
prior<-emaxPrior.control(epmu=5,epsca=20,
							difTargetmu=0,difTargetsca = 20,
							dTarget = 16,p50=4,
							sigmalow=0.001,sigmaup = 1)
mcmc=mcmc.control(chains=3)

bfitout<-fitEmaxB(y=predvals,dose=doselev,prior=prior,modType=3,
									vcest=vcpred,mcmc=mcmc,nproc=3)
#plot(bfitout)
bpred<-predict(bfitout,dosevec=doselev)
best<-bpred$predMed
bse<-bpred$se

test_that("z-stat dose estimated using vcest: 3-parm",{
expect_lt(max(abs((best-poppred)/bse)),
						2.5)
})

### repeat with 4-parm
bfitout<-fitEmaxB(y=predvals,dose=doselev,prior=prior,modType=4,
									vcest=vcpred,mcmc=mcmc,nproc=3)
#plot(bfitout)
bpred<-predict(bfitout,dosevec=doselev)
best<-bpred$predMed
bse<-bpred$se

test_that("z-stat dose estimated using vcest: 4-parm",{
expect_lt(max(abs((best-poppred)/bse)),
						2.5)
})

##############################################
### pbo adjusted
LL<-model.matrix(~ dosefac*visfac-1,preddat)

predpbo<-data.frame(dosefac=factor(rep(0,nd-1),levels=doselev),
										visfac=factor(rep(nv,nd-1),levels=c(1:4)))
predvalspbo<-predvals[2:nd]-predvals[1]
Lpbo<-model.matrix(~ dosefac*visfac-1,predpbo)

LL<-LL[-1,]

vcpredpbo<-(LL-Lpbo)%*%tcrossprod(vcov(modfit),(LL-Lpbo))

bfitoutpbo<-fitEmaxB(y=predvalspbo,dose=doselev[-1],prior=prior,modType=3,
									pboAdj=TRUE,vcest=vcpredpbo,mcmc=mcmc,nproc=3)
#plot(bfitoutpbo)
bpredpbo<-predict(bfitoutpbo,dosevec=doselev[-1])
bestpbo<-bpredpbo$predMed
bsepbo<-bpredpbo$se

test_that("z-stat dose estimated using vcest:pbo-adj 3 parm",{
expect_lt(max(abs((bestpbo-(poppred[-1]-poppred[1]))/bsepbo)),
						2.5)
})

### repeat with 4 parm
bfitoutpbo<-fitEmaxB(y=predvalspbo,dose=doselev[-1],prior=prior,modType=4,
									pboAdj=TRUE,vcest=vcpredpbo,mcmc=mcmc,nproc=3)
#plot(bfitoutpbo)
bpredpbo<-predict(bfitoutpbo,dosevec=doselev[-1])
bestpbo<-bpredpbo$predMed
bsepbo<-bpredpbo$se

test_that("z-stat dose estimated using vcest:pbo-adj 4 parm",{
expect_lt(max(abs((bestpbo-(poppred[-1]-poppred[1]))/bsepbo)),
						2.5)
})


############################
### sim check of vcest option
### continuous data subject to missingness
### same conditions as base case

### currently set for 20 processors
###

runsim<-function(j,seed,nsim){
	set.seed(seed[j])
	nrep<-10
	nd<-6
	nv<-4
	sig<-sqrt(10/4)
	
	doselev<-c(0,1,2,4,8,16)
	dose<-sort(rep(doselev,nrep*nv))
	
	id<-sort(rep(1:(nrep*nd),nv))
	vis<-rep(1:nv,nd*nrep)
	
	dosefac<-factor(dose)
	visfac<-factor(vis)
	
	led50<-log(3)
	emax<-5
	e0<-0
	popparm<-c(led50,emax,e0)
	misprop<-0.15  # mcar
	ntot<-nd*nrep*nv
	
	modmean<-emaxfun(dose,popparm)+vis
	poppred<-tapply(modmean,list(dose,vis),mean)[,nv]
	popcov<-matrix(rep(0.25,nv^2),ncol=nv)
	popcov[1,2]<-0.75
	popcov[2,1]<-0.75
	diag(popcov)<-1
	popcov<-popcov*sig^2
	
	### bayes prior
	prior<-emaxPrior.control(epmu=5,epsca=20,
								difTargetmu=0,difTargetsca = 20,
								dTarget = 16,p50=4,
								sigmalow=0.001,sigmaup = 1)
	mcmc=mcmc.control(chains=1,iter=10000)
		
	
	bdest<-matrix(numeric(nd*nsim),nrow=nsim)
	bdcov<-matrix(numeric(nd*nsim),nrow=nsim)
	for(i in 1:nsim){
		y<-modmean+as.vector(t(rmvnorm(nd*nrep,rep(0,nv),popcov)))
		
		### impose missing data
		misid<-sample(1:ntot,round(ntot*misprop),replace = FALSE)
		y<-y[-misid]
		dosefacm<-dosefac[-misid]
		visfacm<-visfac[-misid]
		vism<-vis[-misid]
		idm<-id[-misid]
		
		modfit<-gls(y ~ dosefacm*visfacm-1, 
								correlation = corSymm(form = ~ vism | idm),
								weights = varIdent(form = ~ 1 | vism))
		
		preddat<-data.frame(dosefacm=factor(doselev),
												visfacm=factor(rep(nv,nd),levels=c(1:4)))
		predvals<-predict(modfit,preddat)
		
		L<-model.matrix(~ dosefacm*visfacm-1,preddat)
		vcpred<-L%*%tcrossprod(vcov(modfit),L)
		
		## bayes model fit
	
		bfitout<-fitEmaxB(y=predvals,dose=doselev,prior=prior,modType=3,
											vcest=vcpred,mcmc=mcmc,nproc=1)
		plot(bfitout)
		bpred<-predict(bfitout,dosevec=doselev)
		bdest[i,]<-bpred$predMed
		
		bdcov[i,]<-(bpred$lb<poppred & bpred$ub>poppred)
	}
	return(list(bdest=bdest,bdcov=bdcov))
}	

nsim<-67
### set up independent stream of random numbers for
### each simulation iteration.
RNGkind("L'Ecuyer-CMRG")
set.seed(12357)
seed<-matrix(integer(nprocdef*7),ncol=7)
seed[1,]<-as.integer(.Random.seed)
for(i in 2:nprocdef){
 seed[i,]<-nextRNGStream(seed[i-1,])
}
 
cl<-makeCluster(nprocdef)
registerDoParallel(cl)	
outsim<-foreach(j=1:nprocdef, .packages=c('nlme','DoseFinding','clinDR')) %dopar%{
	runsim(j,seed,nsim)
}
stopCluster(cl)
RNGkind("default")
		
bdest<-NULL
bdcov<-NULL
for(i in 1:nprocdef){
	bdest<-rbind(bdest,outsim[[i]]$bdest)
	bdcov<-rbind(bdcov,outsim[[i]]$bdcov)
}

test_that("coverage using vcest simulation result",{
expect_lt(abs(max(apply(bdcov-0.9,2,mean))),0.05)
})


#########################################################################
#########################################################################
#########################################################################
##### binary


set.seed(12357)

modType<-4
dvec1<-c(0,.1,.3,.6,1)
dvec2<-c(0,.1,.2,.4,.6,1)
nd1<-length(dvec1)
nd2<-length(dvec2)
n1<-rep(10000,nd1)
n2<-rep(10000,nd2)

parms<-c(log(0.25),1,1.4,-0.5,-0.85)

mlev1<-plogis(emaxfun(dvec1,parms[1:4]))
mlev2<-plogis(emaxfun(dvec2,parms[c(1:3,5)]))

y1<-rbinom(nd1,n1,mlev1)
y2<-rbinom(nd2,n2,mlev2)

### fitEmaxB inputs
y<-c(rep(1,nd1),rep(0,nd1),rep(1,nd2),rep(0,nd2))
counts<-c(y1,n1-y1,y2,n2-y2)
prots<-c(rep(1,2*nd1),rep(2,2*nd2))
dvec<-c(dvec1,dvec1,dvec2,dvec2)

prior<-emaxPrior.control(0,4,0,4,1.0,.5,parmDF=5,binary=TRUE)
mcmc<-mcmc.control(chains=3,warmup=500,iter=3000,seed=53453,propInit=0.15,adapt_delta = .9)

suppressWarnings(testout<-fitEmaxB(y,dvec,modType=modType,
									prot=prots,
									count=counts,binary=TRUE,
									prior=prior,mcmc=mcmc,	
									diagnostics=FALSE,nproc=3))

pgen<-coef(testout)
estimate<-apply(pgen,2,mean)
se<-sqrt(diag(var(pgen)))
z<-(estimate-parms)/se

### check parameter estimates
test_that("model parameters agree within 2se",{
	expect_lt(as.numeric(max(abs(z))),2.0)
})

### check predictions
predout<-predict(testout,dosevec=c(.25,.75),int=2)

poppred<-plogis(emaxfun(c(.25,.75),parms[c(1:3,5)]))
popref<-plogis(emaxfun(0,parms[c(1:3,5)]))
z<-(predout$pred-poppred)/predout$se
zdif<-(predout$fitdif-(poppred-popref))/predout$sedif

### check predictions 
test_that("predictions agree within 2se",{
	expect_lt(as.numeric(max(abs(z))),2.0)
	expect_lt(as.numeric(max(abs(zdif))),2.0)
})

test_that("check absolute levels",{
	expect_that(as.numeric(predout$pred),
							equals(poppred,tol=2*.5/sqrt(10000),scale=1))
})

#########################################################################
#########################################################################
#### repeat for 3 parm use first prot
set.seed(12357)

modType<-3
dvec1<-c(0,.1,.3,.6,1)
dvec2<-c(0,.1,.2,.4,.6,1)
nd1<-length(dvec1)
nd2<-length(dvec2)
n1<-rep(10000,nd1)
n2<-rep(10000,nd2)

parms<-c(log(0.25),1,1.4,-0.5,-0.85)

mlev1<-plogis(emaxfun(dvec1,parms[1:4]))
mlev2<-plogis(emaxfun(dvec2,parms[c(1:3,5)]))

y1<-rbinom(nd1,n1,mlev1)
y2<-rbinom(nd2,n2,mlev2)

### fitEmaxB inputs
y<-c(rep(1,nd1),rep(0,nd1),rep(1,nd2),rep(0,nd2))
counts<-c(y1,n1-y1,y2,n2-y2)
prots<-c(rep(1,2*nd1),rep(2,2*nd2))
dvec<-c(dvec1,dvec1,dvec2,dvec2)

prior<-emaxPrior.control(0,4,0,4,1.0,0.5,parmDF=5,binary=TRUE)
mcmc<-mcmc.control(chains=3,warmup=500,iter=3000,seed=53453,propInit=0.15,adapt_delta = .9)


suppressWarnings(testout<-fitEmaxB(y,dvec,modType=modType,
									prot=prots,
									count=counts,binary=TRUE,
									prior=prior,mcmc=mcmc,	
									diagnostics=FALSE,nproc=3))

pgen<-coef(testout)
estimate<-apply(pgen,2,mean)
se<-sqrt(diag(var(pgen)))
z<-(estimate[1:3]-parms[c(1,3,4)])/se[1:3]

### check parameter estimates
test_that("model parameters agree within 2se",{
	expect_lt(as.numeric(max(abs(z))),2.0)
})

### check predictions
predout<-predict(testout,dosevec=c(.25,.75),int=1)

poppred<-plogis(emaxfun(c(.25,.75),parms[c(1,3,4)]))
popref<-plogis(emaxfun(0,parms[c(1,3,4)]))
z<-(predout$pred-poppred)/predout$se
zdif<-(predout$fitdif-(poppred-popref))/predout$sedif

### check predictions 
test_that("predictions agree within 2se",{
	expect_lt(as.numeric(max(abs(z))),2.0)
	expect_lt(as.numeric(max(abs(zdif))),2.0)
})

test_that("check absolute levels",{
	expect_that(as.numeric(predout$pred),
							equals(poppred,tol=2*.5/sqrt(10000),scale=1))
})



#######################################################
###### check CI and prediction intervals from plot.fitEmaxB
######

runsim<-function(j,seed,nsim){
	set.seed(seed[j])
	
	doselev<-c(0,5,25,50,100,350)
	nd<-length(doselev)
	n<-50*c(78,81,81,81,77,80)
	
	### population parameters for simulation
	e0<-2.465375 
	ed50<-67.481113 
	emax<-15.127726
	sdy<-8.0
	pop.parm<-c(log(ed50),emax,e0)    
	
	modType<-4
	prior<-emaxPrior.control(0,30,0,30,350,50,0.1,30,parmDF=5)
	mcmc<-mcmc.control(chains=1,warmup=500,iter=5000,seed=53453,propInit=0.15,adapt_delta = .9)
	estan<-selEstan('mrmodel')
	
	dose<-rep(doselev,n)
	meanlev<-emaxfun(doselev,pop.parm)  
	meanrep<-emaxfun(dose,pop.parm)  
	clev<-0.9
	covci<-matrix(logical(nsim*nd),ncol=nd)
	covpi<-matrix(logical(nsim*nd),ncol=nd)
	covdifci<-matrix(logical(nsim*nd),ncol=nd)
	covdifpi<-matrix(logical(nsim*nd),ncol=nd)
	for (i in 1:nsim){
		y<-rnorm(sum(n),meanrep,sdy)
		ymean<-tapply(y,dose,mean)
		msSat<-(summary(lm(y~factor(dose)))$sigma)^2
		suppressWarnings(testout<-fitEmaxB(ymean,doselev,prior=prior,modType=modType,count=n,
											mcmc=mcmc,msSat=msSat,
											estan=estan,diagnostics = FALSE,nproc = 1))
		if(is.null(testout)){
			covci[i,]<-NA
			covpi[i,]<-NA
			covdifci[i,]<-NA
			covdifpi[i,]<-NA   
		}else{
			intout<-plot(testout,clev=clev,plot=FALSE)$plotdata
			covci[i,]<-meanlev>=intout[,'cil'] & meanlev<=intout[,'cih']
			### coverage for independent sample means from the same design
			y<-rnorm(sum(n),meanrep,sdy)
			ym<-tapply(y,dose,mean)
			covpi[i,]<-ym>=intout[,'pil'] & ym<=intout[,'pih']  
			### difference with pbo
			intout<-plot(testout,clev=clev,plot=FALSE,plotDif=TRUE)$plotdata
			covdifci[i,]<-meanlev-meanlev[1]>=intout[,'cil'] & meanlev-meanlev[1]<=intout[,'cih']
			### coverage for independent sample means from the same design
			ym<-ym-ym[1]
			covdifpi[i,]<-ym>=intout[,'pil'] & ym<=intout[,'pih']   
		}
	}
	return(outsim=list(covci=covci,covpi=covpi,covdifci=covdifci,covdifpi=covdifpi))
}

clev<-0.9
nsim<-67
### set up independent stream of random numbers for
### each simulation iteration.
RNGkind("L'Ecuyer-CMRG")
set.seed(12357)
seed<-matrix(integer(nprocdef*7),ncol=7)
seed[1,]<-as.integer(.Random.seed)
for(i in 2:nprocdef){
 seed[i,]<-nextRNGStream(seed[i-1,])
}
 
cl<-makeCluster(nprocdef)
registerDoParallel(cl)	
outsim<-foreach(j=1:nprocdef, .packages=c('clinDR')) %dopar%{
	runsim(j,seed,nsim)
}
stopCluster(cl)
RNGkind("default")

covci<-NULL
covpi=NULL
covdifci<-NULL
covdifpi=NULL
for(i in 1:nprocdef){
	covci<-rbind(covci,outsim[[i]]$covci)
	covpi<-rbind(covpi,outsim[[i]]$covpi)
	covdifci<-rbind(covdifci,outsim[[i]]$covdifci)
	covdifpi<-rbind(covdifpi,outsim[[i]]$covdifpi)
}

nsim<-nrow(covci)
test_that("plot.fitEmaxB CI for continuous data agree within 3se",{
	expect_that(clev,
							equals(as.numeric(mean(apply(covci,2,mean,na.rm=TRUE))),
										 tolerance=0.05,scale=1))
})
test_that("plot.fitEmaxB PI for continuous data agree within 3se",{
	expect_that(clev,
							equals(as.numeric(mean(apply(covpi,2,mean,na.rm=TRUE))),
										 tolerance=0.05,scale=1))
})


test_that("plot.fitEmaxB CI DIF for continuous data agree within 3se",{
	expect_that(clev,
							equals(as.numeric(mean(apply(covdifci[,-1],2,mean,na.rm=TRUE))),
										 tolerance=0.05,scale=1))
})
test_that("plot.fitEmaxB PI DIF for continuous data agree within 3se",{
	expect_that(clev,
							equals(as.numeric(mean(apply(covdifpi[,-1],2,mean,na.rm=TRUE))),
										 tolerance=0.05,scale=1))
})

#####################################
### repeat with binary data

runsim<-function(j,seed,nsim){
	set.seed(seed[j])
	
	doselev<-c(0,5,25,50,100,350)
	nd<-length(doselev)
	n<-10*c(78,81,81,81,77,80)
	
	### population parameters for simulation
	e0<-qlogis(.2) 
	ed50<-67.481113 
	emax<-qlogis(.95)
	pop.parm<-c(log(ed50),emax,e0)    
	meanlev<-plogis(emaxfun(doselev,pop.parm))
	y01<-c(rep(1,length(meanlev)),rep(0,length(meanlev)))
	d01<-c(doselev,doselev)
	clev<-0.9
	covci<-matrix(logical(nsim*nd),ncol=nd)
	covpi<-matrix(logical(nsim*nd),ncol=nd)
	covdifci<-matrix(logical(nsim*nd),ncol=nd)
	covdifpi<-matrix(logical(nsim*nd),ncol=nd)
	
	modType<-4
	prior<-emaxPrior.control(0,4,0,4,350,50,parmDF=5,binary=TRUE)
	mcmc<-mcmc.control(chains=1,warmup=500,iter=5000,seed=53453,propInit=0.15,adapt_delta = .9)
	estan<-selEstan('mrmodel')
	
	for (i in 1:nsim){
		y<-rbinom(length(n),n,meanlev)
		n01<-c(y,n-y)
		suppressWarnings(testout<-fitEmaxB(y01,d01,prior,modType=4,
											count=n01,diagnostics = FALSE,binary=TRUE,nproc=1))
		if(is.null(testout)){
			covci[i,]<-NA
			covpi[i,]<-NA
			covdifci[i,]<-NA
			covdifpi[i,]<-NA   
		}else{
			intout<-plot(testout,clev=clev,plot=FALSE)$plotdata
			covci[i,]<-meanlev>=intout[,'cil'] & meanlev<=intout[,'cih']
			### coverage for independent sample means from the same design
			
			ypred<-rbinom(length(n),n,meanlev)
			ym<-ypred/n
			covpi[i,]<-ym>=intout[,'pil'] & ym<=intout[,'pih']  
			### difference with pbo
			intout<-plot(testout,clev=clev,plot=FALSE,plotDif=TRUE)$plotdata
			covdifci[i,]<-meanlev-meanlev[1]>=intout[,'cil'] & meanlev-meanlev[1]<=intout[,'cih']
			### coverage for independent sample means from the same design
			ym<-ym-ym[1]
			covdifpi[i,]<-ym>=intout[,'pil'] & ym<=intout[,'pih']   
		}
	}
	return(outsim=list(covci=covci,covpi=covpi,covdifci=covdifci,covdifpi=covdifpi))
}

clev<-0.9
nsim<-67
### set up independent stream of random numbers for
### each simulation iteration.
RNGkind("L'Ecuyer-CMRG")
set.seed(12357)
seed<-matrix(integer(nsim*7),ncol=7)
seed[1,]<-as.integer(.Random.seed)
for(i in 2:nsim){
 seed[i,]<-nextRNGStream(seed[i-1,])
}
 
cl<-makeCluster(nprocdef)
registerDoParallel(cl)	
outsim<-foreach(j=1:nprocdef, .packages=c('clinDR')) %dopar%{
	runsim(j,seed,nsim)
}
stopCluster(cl)
RNGkind("default")

covci<-NULL
covpi=NULL
covdifci<-NULL
covdifpi=NULL
for(i in 1:nprocdef){
	covci<-rbind(covci,outsim[[i]]$covci)
	covpi<-rbind(covpi,outsim[[i]]$covpi)
	covdifci<-rbind(covdifci,outsim[[i]]$covdifci)
	covdifpi<-rbind(covdifpi,outsim[[i]]$covdifpi)
}

test_that("plot.fitEmaxB CI for binary data agree within 3se",{
	expect_that(clev,
							equals(as.numeric(mean(apply(covci,2,mean,na.rm=TRUE))),
										 tolerance=0.04,scale=1))
})
test_that("plot.fitEmaxB CI for binary data agree within 3se",{
	expect_that(clev,
							equals(as.numeric(mean(apply(covci,2,mean,na.rm=TRUE))),
										 tolerance=0.04,scale=1))
})
test_that("plot.fitEmaxB PI for binary data agree within 3se",{
	expect_that(clev,
							equals(as.numeric(mean(apply(covpi,2,mean,na.rm=TRUE))),
										 tolerance=0.05,scale=1))
})
test_that("plot.fitEmaxB PI for binary data agree within 3se",{
	expect_that(clev,
							equals(as.numeric(mean(apply(covpi,2,mean,na.rm=TRUE))),
										 tolerance=0.05,scale=1))
})


test_that("plot.fitEmaxB CI DIF for binary data agree within 3se",{
	expect_that(clev,
							equals(as.numeric(mean(apply(covdifci[,-1],2,mean,na.rm=TRUE))),
										 tolerance=0.05,scale=1))
})
test_that("plot.fitEmaxB CI DIF for binary data agree within 3se",{
	expect_that(clev,
							equals(as.numeric(mean(apply(covdifci[,-1],2,mean,na.rm=TRUE))),
										 tolerance=0.04,scale=1))
})
test_that("plot.fitEmaxB PI DIF for binary data agree within 3se",{
	expect_that(clev,
							equals(as.numeric(mean(apply(covdifpi[,-1],2,mean,na.rm=TRUE))),
										 tolerance=0.05,scale=1))
})
test_that("plot.fitEmaxB PI DIF for binary data agree within 3se",{
	expect_that(clev,
							equals(as.numeric(mean(apply(covdifpi[,-1],2,mean,na.rm=TRUE))),
										 tolerance=0.04,scale=1))
})

#################################################################################
#################################################################################
#### include covariates
#################################################################################

## 1 covariate, 2 protocols
set.seed(12357)

doselev<-c(0,5,25,50,100,350)
n<-c(78,81,81,81,77,80)
n1<-sum(n)
n2<-sum(n[1:4])

doselev<-c(doselev,doselev[1:4])
n<-c(n,n[1:4])

### population parameters for simulation
e0<-2.465375 
ed50<-67.481113 
emax<-15.127726
sdy<-8.0
x1<-rnorm(n1)
x1<-x1-mean(x1)
x2<-rnorm(n2)
x2<-x2-mean(x2)
x<-matrix(c(x1,x2),ncol=1)
pop<-c(log(ed50),emax,e0)    
dose<-rep(doselev,n)
bparm<-1
meanlev<-emaxfun(dose,pop) + x%*%bparm 

y<-rnorm(n1+n2,meanlev,sdy)
prots<-c(rep(1,n1),rep(2,n2))

basemu<-0
basevar<-matrix((10*sdy)^2,nrow=1,ncol=1)
prior<-emaxPrior.control(0,30,0,30,350,50,0.1,30,parmDF=5,basemu=basemu,basevar=basevar)
mcmc<-mcmc.control(chains=3,warmup=500,iter=3000,seed=53453,propInit=0.15,adapt_delta = .9)

suppressWarnings(testout<-fitEmaxB(y,dose,prior=prior,modType=4,prot=prots,xbase=x,
									mcmc=mcmc,diagnostics=FALSE,nproc=3))

parms<-coef(testout)
estimate<-apply(parms,2,mean)
se<-sqrt(diag(var(parms)))
z<-(estimate-c(pop[1],1,pop[2:3],pop[3],bparm))/se

### check parameter estimates
test_that("model parameters agree within 2.5se",{
	expect_lt(as.numeric(max(abs(z))),2.5)
})

### check predictions
predout<-predict(testout,dosevec=c(20,80),int=2)

poppred<-emaxfun(c(20,80),pop[c(1:3)])
z<-(predout$pred-poppred)/predout$se
zdif<-(predout$fitdif-(poppred-e0))/predout$sedif


### check predictions 
test_that("predictions agree within 2.5se",{
	expect_lt(as.numeric(max(abs(z))),2.5)
	expect_lt(as.numeric(max(abs(zdif))),2.5)
})

test_that("check absolute levels",{
	expect_that(as.numeric(predout$pred),
							equals(poppred,tol=2*sdy/sqrt(70),scale=1))
	expect_that(as.numeric(predout$fitdif),
							equals((poppred-e0),tol=2*sdy/sqrt(70),scale=1))
})

#############################################################
## 3 covariates, 2 protocols
set.seed(12357)

doselev<-c(0,5,25,50,100,350)
n<-c(78,81,81,81,77,80)
n1<-sum(n)
n2<-sum(n[1:4])

doselev<-c(doselev,doselev[1:4])
n<-c(n,n[1:4])

### population parameters for simulation
e0<-2.465375 
ed50<-67.481113 
emax<-15.127726
sdy<-8.0
x1<-matrix(rnorm(3*n1),ncol=3)
x1<-scale(x1,center=TRUE,scale=FALSE)
x2<-matrix(rnorm(3*n2),ncol=3)
x2<-scale(x2,center=TRUE,scale=FALSE)
x<-rbind(x1,x2)
pop<-c(log(ed50),emax,e0)    
dose<-rep(doselev,n)
bparm<-c(2,-1,0.5)
meanlev<-emaxfun(dose,pop) + x%*%bparm 

y<-rnorm(n1+n2,meanlev,sdy)
prots<-c(rep(1,n1),rep(2,n2))

basemu<-numeric(3)
basevar<-diag(3)*(10*sdy)^2
prior<-emaxPrior.control(0,30,0,30,350,50,0.1,30,parmDF=5,basemu=basemu,basevar=basevar)
mcmc<-mcmc.control(chains=3,warmup=500,iter=3000,seed=53453,propInit=0.15,adapt_delta = .9)

suppressWarnings(testout2<-fitEmaxB(y,dose,prior=prior,modType=4,prot=prots,xbase=x,
									 mcmc=mcmc,diagnostics=FALSE,nproc=3))

parms<-coef(testout2)
estimate<-apply(parms,2,mean)
se<-sqrt(diag(var(parms)))
z<-(estimate-c(pop[1],1,pop[2:3],pop[3],bparm))/se

### check parameter estimates
test_that("model parameters agree within 2.5se",{
	expect_lt(as.numeric(max(abs(z))),2.5)
})

### check predictions
predout<-predict(testout2,dosevec=c(20,80),int=2)

poppred<-emaxfun(c(20,80),pop[c(1:3)])
z<-(predout$pred-poppred)/predout$se
zdif<-(predout$fitdif-(poppred-e0))/predout$sedif


### check predictions 
test_that("predictions agree within 2.5se",{
	expect_lt(as.numeric(max(abs(z))),2.5)
	expect_lt(as.numeric(max(abs(zdif))),2.5)
})

test_that("check absolute levels",{
	expect_that(as.numeric(predout$pred),
							equals(poppred,tol=2*sdy/sqrt(70),scale=1))
	expect_that(as.numeric(predout$fitdif),
							equals((poppred-e0),tol=2*sdy/sqrt(70),scale=1))
})

##########################################################
### check with larger n for better asymptotics
### 3-parm model, 2 covariates, 1 protocol

set.seed(12357)

doselev<-c(0,5,25,50,100,350)
n<-5*c(78,81,81,81,77,80)
ntot<-sum(n)

### population parameters for simulation
e0<-2.465375 
ed50<-67.481113 
emax<-15.127726
sdy<-8.0
pop<-c(log(ed50),emax,e0)    
dose<-rep(doselev,n)
meanlev<-emaxfun(dose,pop)  

x<-matrix(rnorm(2*ntot),ncol=2)
x<-scale(x,center=TRUE,scale=FALSE)
bparm<-c(2,-1)
meanlev<-meanlev + x%*%bparm 

y<-rnorm(ntot,meanlev,sdy)

basemu<-numeric(2)
basevar<-diag(2)*(10*sdy)^2
prior<-emaxPrior.control(0,30,0,30,350,50,0.1,30,parmDF=5,basemu=basemu,basevar=basevar)
mcmc<-mcmc.control(chains=1,warmup=500,iter=3000,seed=53453,propInit=0.15,adapt_delta = .9)

suppressWarnings(testout3<-fitEmaxB(y,dose,prior=prior,modType=3,xbase=x,
									mcmc=mcmc,diagnostics=FALSE,nproc=1))

parms<-coef(testout3)
estimate<-apply(parms,2,mean)
se<-sqrt(diag(var(parms)))
z<-(estimate-c(pop[1],pop[2:3],bparm))/se

### check parameter estimates
test_that("model parameters agree within 2.5se",{
	expect_lt(as.numeric(max(abs(z))),2.5)
})

### check predictions
predout<-predict(testout3,dosevec=c(20,80),int=1)

poppred<-emaxfun(c(20,80),pop[c(1:3)])
z<-(predout$pred-poppred)/predout$se
zdif<-(predout$fitdif-(poppred-e0))/predout$sedif


### check predictions 
test_that("predictions agree within 2.5se",{
	expect_lt(as.numeric(max(abs(z))),2.5)
	expect_lt(as.numeric(max(abs(zdif))),2.5)
})

test_that("check absolute levels",{
	expect_that(as.numeric(predout$pred),
							equals(poppred,tol=2*sdy/sqrt(70),scale=1))
	expect_that(as.numeric(predout$fitdif),
							equals((poppred-e0),tol=2*sdy/sqrt(70),scale=1))
})

##########################################################
### binary, covariates
### check with larger n for better asymptotics
### 4-parm model, 2 covariates, 1 protocol

set.seed(12357)

doselev<-c(0,5,25,50,100,350)
n<-5*c(78,81,81,81,77,80)
ntot<-sum(n)

### population parameters for simulation
e0<- -1.5 
ed50<-67.481113 
emax<-4.0
pop<-c(log(ed50),emax,e0)    
dose<-rep(doselev,n)
meanlev<-emaxfun(dose,pop)  

x<-matrix(rnorm(2*ntot),ncol=2)
x<-scale(x,center=TRUE,scale=FALSE)
bparm<-c(2,-1)
meanlev<-plogis(meanlev + x%*%bparm)

y<-rbinom(ntot,1,meanlev)

basemu<-numeric(2)
basevar<-diag(2)*(4)^2
prior<-emaxPrior.control(0,30,0,30,350,50,parmDF=5,basemu=basemu,basevar=basevar,binary=TRUE)
mcmc<-mcmc.control(chains=1,warmup=500,iter=3000,seed=53453,propInit=0.15,adapt_delta = .9)

suppressWarnings(testout4b<-fitEmaxB(y,dose,prior=prior,modType=4,xbase=x,
									mcmc=mcmc,diagnostics=FALSE,binary=TRUE,nproc=1))

parms<-coef(testout4b)
estimate<-apply(parms,2,mean)
se<-sqrt(diag(var(parms)))
z<-(estimate-c(pop[1],1,pop[2:3],bparm))/se

### check parameter estimates
test_that("model parameters agree within 2.5se",{
	expect_lt(as.numeric(max(abs(z))),2.5)
})

### check predictions
predout<-predict(testout4b,dosevec=c(20,80),int=1,xvec=c(0,0))

poppred<-plogis(emaxfun(c(20,80),pop[c(1:3)]))
z<-(predout$pred-poppred)/predout$se
zdif<-(predout$fitdif-(poppred-plogis(e0)))/predout$sedif


### check predictions 
test_that("predictions agree within 2.5se",{
	expect_lt(as.numeric(max(abs(z))),2.5)
	expect_lt(as.numeric(max(abs(zdif))),2.5)
})

test_that("check absolute levels",{
	expect_that(as.numeric(predout$pred),
							equals(poppred,tol=0.05,scale=1))
	expect_that(as.numeric(predout$fitdif),
							equals((poppred-plogis(e0)),tol=0.05,scale=1))
})

### binary, covariates
### check with larger n for better asymptotics
### 3-parm model, 2 covariates, 1 protocol

set.seed(20572)

doselev<-c(0,5,25,50,100,350)
n<-4*c(78,81,81,81,77,80)
ntot<-sum(n)

### population parameters for simulation
e0<- -1.5 
ed50<-67.481113 
emax<-4.0
pop<-c(log(ed50),emax,e0)    
dose<-c(rep(doselev,n/2),rep(doselev,n/2))
prot<-sort(rep(1:2,ntot/2))
meanlev<-emaxfun(dose,pop)+0.5*(prot==2)  

x1<-matrix(rnorm(ntot),ncol=2)
x1<-scale(x1,center=TRUE,scale=FALSE)
x2<-matrix(rnorm(ntot),ncol=2)
x2<-scale(x2,center=TRUE,scale=FALSE)
x<-rbind(x1,x2)
bparm<-c(1.,-0.5)
meanlev<-plogis(meanlev + x%*%bparm)

y<-rbinom(ntot,1,meanlev)

basemu<-numeric(2)
basevar<-diag(2); basevar[2,1]<-.25; basevar[1,2]<-.25
basevar<-basevar*(4)^2  ## off-diagonal elements

prior<-emaxPrior.control(0,30,0,30,350,50,parmDF=5,basemu=basemu,basevar=basevar,binary=TRUE)
mcmc<-mcmc.control(chains=1,warmup=500,iter=3000,seed=53453,propInit=0.15,adapt_delta = .9)

suppressWarnings(testout5b<-fitEmaxB(y,dose,prot=prot,prior=prior,modType=3,xbase=x,
									mcmc=mcmc,diagnostics=FALSE,binary=TRUE,nproc=1))

parms<-coef(testout5b)
estimate<-apply(parms,2,mean)
se<-sqrt(diag(var(parms)))
z<-(estimate-c(pop[1],pop[c(2:3)],pop[3]+0.5,bparm))/se

### check parameter estimates
test_that("model parameters agree within 2.5se",{
	expect_lt(as.numeric(max(abs(z))),2.5)
})

### check predictions
predout<-predict(testout5b,dosevec=c(20,80),int=1,xvec=c(0,0))

poppred<-plogis(emaxfun(c(20,80),c(pop[c(1:2)],e0)))
z<-(predout$pred-poppred)/predout$se
zdif<-(predout$fitdif-(poppred-plogis(e0)))/predout$sedif
### check predictions 
test_that("predictions agree within 2.5se",{
	expect_lt(as.numeric(max(abs(z))),2.5)
	expect_lt(as.numeric(max(abs(zdif))),2.5)
})

test_that("check absolute levels",{
	expect_that(as.numeric(predout$pred),
							equals(poppred,tol=0.05,scale=1))
	expect_that(as.numeric(predout$fitdif),
							equals((poppred-plogis(e0)),tol=0.05,scale=1))
})

### check predictions
predout<-predict(testout5b,dosevec=c(20,80),int=2,xvec=c(0,0))

poppred<-plogis(emaxfun(c(20,80),c(pop[c(1:2)],e0+0.5)))
z<-(predout$pred-poppred)/predout$se
zdif<-(predout$fitdif-(poppred-plogis(e0+0.5)))/predout$sedif


### check predictions 
test_that("predictions agree within 2.5se",{
	expect_lt(as.numeric(max(abs(z))),2.5)
	expect_lt(as.numeric(max(abs(zdif))),2.5)
})

test_that("check absolute levels",{
	expect_that(as.numeric(predout$pred),
							equals(poppred,tol=0.05,scale=1))
	expect_that(as.numeric(predout$fitdif),
							equals((poppred-plogis(e0+0.5)),tol=0.05,scale=1))
})



###################################
###### binary data (with missingness)
###### applied with vcest specified
###### for 2-stage model fit

set.seed(123578)
nrep<-100
nd<-6
nv<-4
ntot<-nd*nrep*nv

doselev<-c(0,1,2,4,8,16)
dose<-sort(rep(doselev,nrep*nv))
dindlev<-1:nd
dind<-sort(rep(dindlev,nrep*nv))

id<-sort(rep(1:(nrep*nd),nv))
vis<-rep(1:nv,nd*nrep)
k<-(dind-1)*nv+vis

led50<-log(3)
emax<-qlogis(0.75)-qlogis(0.25)
e0<-qlogis(0.25)
popparm<-c(led50,emax,e0)
tau<-1.0
misprop<-0.15  # mcar

modmeanL<-emaxfun(dose,popparm)+0.5*(vis-1)
poppredL<-tapply(modmeanL,list(vis,dose),mean) # reverse to match binary code
poppred<-matrix(numeric(nd*nv),nrow=nv)

fint<-function(theta,logit,sig){
	plogis(logit+theta)*dnorm(theta,mean=0,sd=sig)
}
fint<-Vectorize(fint,vectorize.args='theta')
for(j in 1:nd){
	for(i in 1:nv){
		poppred[i,j]<-integrate(fint,-Inf,Inf,logit=poppredL[i,j],sig=tau)$value
	}
}

theta<-rnorm(nd*nrep,0,tau)
theta<-rep(theta,rep(nv,nd*nrep))
y<-rbinom(ntot,1,plogis(modmeanL+theta))
dosefac<-factor(dose)
visfac<-factor(vis)

### impose missing data
misid<-sample(1:ntot,round(ntot*misprop),replace = FALSE)
y<-y[-misid]
dosefac<-dosefac[-misid]
dind<-dind[-misid]
visfac<-visfac[-misid]
vis<-vis[-misid]
id<-id[-misid]
k<-k[-misid]

checkid<-(length(unique(id))==nd*nrep) ## must be true or skip fixed


#### fit saturated model (replace by package eventually)
fitModel<-function(id,y,trt,visit,prmean0,prsd0,prmean,prsd,
									 gparm1=3,gparm2=1.5,
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

mcmc<-list(chains=3,thin=4,warmup=2000,iter=4*3333+2000,
					 adapt_delta=0.95,seed=12357)


modfit<-fitModel(id,y,dind,vis,
				 prmean0=rep(qlogis(0.25),nv),prsd0=rep(3.0,nv),
				 prmean=rep(0,nv),prsd=rep(3.0,nv),mcmc=mcmc)

outmod<-modfit$stanfit
beta<-modfit$beta
sigma<-modfit$sigma
theta<-modfit$theta
nmc<-length(sigma)

##########################################################
### convergence diagnostics
# traceplot(outmod,pars=c('beta'))
# traceplot(outmod,pars=c('sigma'))
# traceplot(outmod,pars=c('theta[1]','theta[2]','theta[3]','theta[4]'))
# 
# stan_rhat(outmod,pars=c('beta'))
# stan_rhat(outmod,pars='theta')
# Rhat(as.matrix(outmod,pars='sigma'))
# stan_ac(outmod,pars=c('beta'))
# stan_ac(outmod,pars=c('sigma'))
# 
# 
# stan_mcse(outmod)
# stan_ess(outmod)


j<-4
lpm<-matrix(numeric(nmc*nd),ncol=nd)
for(i in 1:nd){
	for(k in 1:nmc){
		lpm[k,i]<-qlogis(integrate(fint,-Inf,Inf,logit=beta[k,j,i],
															 sig=sigma[k])$value)
	}
}

predvals<-apply(lpm,2,mean)
vcpred<-var(lpm)

## bayes model fit
prior<-emaxPrior.control(epmu=qlogis(0.5),epsca=2,
							difTargetmu=0,difTargetsca = 2,
							dTarget = 16,p50=4,binary=TRUE)
mcmc=mcmc.control(chains=3)

bfitoutb<-fitEmaxB(y=predvals,dose=doselev,prior=prior,modType=3,
									vcest=vcpred,mcmc=mcmc,nproc=3,binary=TRUE)
plot(bfitoutb)
bpredb<-predict(bfitoutb,dosevec=doselev)

test_that("binary vcest-specified ci cov",{
expect_true(all(bpredb$lb<poppred[nv,] & bpredb$ub>poppred[nv,]))
})

test_that("binary vcest-specified pred check",{
expect_that(as.numeric(bpredb$predMed),
						equals(poppred[nv,],tol=0.04,scale=1))
})


