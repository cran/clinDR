context('emaxsimB evaluations updated prior')

####################################################################
###### check mcp testing
######
set.seed(12357)
nsim<-1000

doselev<-c(0,0.033,0.067,0.167,0.33,1)
tareff<-0.313
ndose<-length(doselev)
n<-rep(95,ndose)
sdy<-1

parm.mat<-matrix(c( .5,.5,  .05,.5, .95,.5),byrow=T,ncol=2)

parm.mat[,1]<-qgamma(parm.mat[,1],shape=1.05,rate=3.0)  ###ed50
parm.mat[,2]<-qgamma(parm.mat[,2],shape=2.5,rate=1.18) ###lambda


###########################################
#### use pop emax model 
e0<-0
ed50<-qgamma(0.5,shape=1.05,rate=3.0)
lambda<-qgamma(0.5,shape=2.5,rate=1.18)
emax<-(ed50^lambda+1)*(tareff)
pop<-c(log(ed50),lambda,emax,e0)
meanlev<-emaxfun(doselev,parm=pop)

gen.parm<-FixedMean(n,doselev,meanlev,sdy,parm=pop)  

### 3 contrasts and linear
emaxMods3<-Mods(sigEmax=cbind(parm.mat[,1],parm.mat[,2]),
								linear=NULL,doses=doselev,placEff=e0,maxEff=1)
emaxMat3<-optContr(emaxMods3,w=n)

#### power from DoseFinding
altmod<-Mods(sigEmax=cbind(parm.mat[1,1],parm.mat[1,2]), doses=doselev, 
						 placEff = e0, maxEff=max(meanlev),
						 direction ="increasing")

pow3<-as.numeric(powMCT(emaxMat3,alpha=0.05,altModels=altmod,n=n, 
												sigma=sdy,placAdj=FALSE,
												alternative="one.sided",critV=TRUE))

prior<-emaxPrior.control(0,3,0,3,0.5,1.0,0.1,5,parmDF=5)
mcmc<-mcmc.control(chains=1,warmup=500,iter=5000,seed=53453,propInit=0.15,adapt_delta = .95)


D3 <- emaxsimB(nsim,gen.parm,prior,modType=4,
							 mcmc=mcmc,testMods=emaxMods3,check=FALSE,nproc=nprocdef)


test_that("check mcp power calculation",{
	expect_that(mean(D3$pVal<0.05),equals(pow3,tol=2.0*sqrt(pow3*(1-pow3)/nsim)))
})
####################################################################

####################################################################
### large sample near-exact with negative emax  4-parm
set.seed(12357)
nsim<-1000

doselev<-c(0,0.033,0.067,0.167,0.33,1)
tareff<- -0.313
ndose<-length(doselev)
n<-rep(95,ndose)
sdy<-0.01

parm.mat<-matrix(c( .5,.5,  .05,.5, .95,.5),byrow=T,ncol=2)

parm.mat[,1]<-qgamma(parm.mat[,1],shape=1.05,rate=3.0)  ###ed50
parm.mat[,2]<-qgamma(parm.mat[,2],shape=2.5,rate=1.18) ###lambda

###########################################
#### use pop emax model 
e0<-0
ed50<-qgamma(0.5,shape=1.05,rate=3.0)
lambda<-qgamma(0.5,shape=2.5,rate=1.18)
emax<-(ed50^lambda+1)*(tareff)
pop<-c(log(ed50),lambda,emax,e0)
meanlev<-emaxfun(doselev,parm=pop)

gen.parm<-FixedMean(n,doselev,meanlev,sdy,parm=pop)  

### 1 contrast
emaxMods1<-Mods(sigEmax=cbind(parm.mat[1,1],parm.mat[1,2]),
								doses=doselev,placEff=e0,maxEff=-1)
emaxMat1<-optContr(emaxMods1,w=n)

prior<-emaxPrior.control(0,3,0,3,1.0,0.5,0.001,0.1,parmDF=5)
mcmc<-mcmc.control(chains=1,warmup=500,iter=5000,seed=53453,propInit=0.15,adapt_delta = .95)

D1 <- emaxsimB(nsim,gen.parm,prior,modType=4,
							mcmc=mcmc,testMods=emaxMods1,nproc=nprocdef)

parmmed<-apply(coef(D1),2,mean)

test_that("check asymptotic distributions for 4-parm model",{
	expect_that(as.numeric(parmmed),equals(pop,tol=0.01))
})

pop1<-c(pop,emaxfun(1.0,pop))
popmat<-matrix(rep(pop1,nsim),ncol=length(pop1),byrow=T)

cov95<-apply(t(D1$estlb[,,1])<popmat & t(D1$estub[,,1])>popmat,2,mean)
test_that("check posterior 95 bounds for model parameter",{
	expect_that(as.numeric(cov95),equals(rep(0.95,length(pop1)),tol=0.01))
})
cov90<-apply(t(D1$estlb[,,2])<popmat & t(D1$estub[,,2])>popmat,2,mean)
test_that("check posterior 90 bounds for model parameter",{
	expect_that(as.numeric(cov90),equals(rep(0.90,length(pop1)),tol=0.015))
})
cov80<-apply(t(D1$estlb[,,3])<popmat & t(D1$estub[,,3])>popmat,2,mean)
test_that("check posterior 80 bounds for model parameter",{
	expect_that(as.numeric(cov80),equals(rep(0.80,length(pop1)),tol=0.015))
})


test_that("emaxsim internal predicted values match population",{
	expect_that(max(apply(D1$fitpredv,2,mean)),
							equals(max(emaxfun(doselev,pop)),tol=0.01,scale=1))
})



####################################################################
##########  binary
####################################################################

####################################################################
###### check mcp testing
######
set.seed(12357)
nsim<-1000

doselev<-c(0,0.033,2*0.033,4*0.033,8*0.033,16*0.033)
ndose<-length(doselev)
n<-rep(40,ndose)

parm.mat<-matrix(c( .5,.5,  .05,.5, .95,.5),byrow=T,ncol=2)

parm.mat[,1]<-qgamma(parm.mat[,1],shape=1.05,rate=3.0)  ###ed50
parm.mat[,2]<-rep(1,3) ###lambda


###########################################
#### use pop emax model 
e0<-qlogis(.2)
ed50<-0.15
lambda<-1
emax<-qlogis(0.55)-e0
pop<-c(log(ed50),lambda,emax,e0)
meanlev<-plogis(emaxfun(doselev,parm=pop))

gen.parm<-FixedMean(n,doselev,meanlev,sdy,parm=pop,binary=TRUE)  

### variance matrix for logistic rates
V<-diag(1/(meanlev*(1-meanlev)*n))

### 3 contrasts and linear
emaxMods3<-Mods(sigEmax=cbind(parm.mat[,1],parm.mat[,2]),
								linear=NULL,doses=doselev,placEff=e0,maxEff=1)
emaxMat3<-optContr(emaxMods3,S=V)

#### power from DoseFinding
altmod<-Mods(sigEmax=cbind(parm.mat[1,1],parm.mat[1,2]), doses=doselev, 
						 placEff = e0, maxEff=max(qlogis(meanlev)-qlogis(meanlev[1])),
						 direction ="increasing")

pow3<-as.numeric(powMCT(emaxMat3,alpha=0.05,altModels=altmod, 
												S=V,df=Inf,placAdj=FALSE,
												alternative="one.sided",critV=TRUE))

prior<-emaxPrior.control(qlogis(0.2),4,0,4,16*0.033,0.05,parmDF=5,binary=TRUE)
mcmc<-mcmc.control(chains=1,warmup=500,iter=5000,seed=53453,propInit=0.15,adapt_delta = .95)


D3 <- emaxsimB(nsim,gen.parm,prior,modType=4,
							 mcmc=mcmc,testMods=emaxMods3,check=FALSE,binary=TRUE,
							 seed=12357,nproc=nprocdef)

test_that("check mcp power calculation",{
	expect_that(mean(D3$pVal<0.05),equals(pow3,tol=0.02,scale=1))
})

test_that("emaxsim internal predicted values match population for binary",{
	expect_that(max(apply(D3$fitpredv,2,mean)),
							equals(max(plogis(emaxfun(doselev,pop))),tol=0.03,scale=1))
})

###
### 3 parameter fit, different format for testing contrasts
D3m <- emaxsimB(nsim,gen.parm,prior,modType=3,ed50contr=parm.mat[,1],
							lambdacontr=parm.mat[,2],mcmc=mcmc,binary=TRUE,
							seed=12357,nproc=nprocdef)

test_that("check mcp power calculation",{
	expect_that(mean(D3m$pVal<0.05),equals(pow3,tol=0.02,scale=1))
})

test_that("emaxsim internal predicted values match population for binary",{
	expect_that(max(apply(D3m$fitpredv,2,mean)),
							equals(max(plogis(emaxfun(doselev,pop))),tol=0.02,scale=1))
})



####################################################################
### simulation to test prediction and parameters
set.seed(12357)
nsim<-1000

doselev<-c(0,0.033,2*0.033,4*0.033,8*0.033,16*0.033)
ndose<-length(doselev)
n<-10*rep(95,ndose)
###########################################
#### use pop emax model 
e0<-qlogis(.2)
ed50<-0.15
lambda<-1
emax<-qlogis(0.55)-e0
pop<-c(log(ed50),lambda,emax,e0)
meanlev<-plogis(emaxfun(doselev,parm=pop))

gen.parm<-FixedMean(n,doselev,meanlev,sdy,parm=pop,binary=TRUE)  

prior<-emaxPrior.control(qlogis(0.2),4,0,4,16*0.033,0.05,parmDF=5,binary=TRUE)
mcmc<-mcmc.control(chains=1,warmup=500,iter=5000,seed=53453,propInit=0.15,adapt_delta = .95)

### ci for a dose not in design
customCode<-function(parms,residSD,pVal,dose,y,customParms){
	dtar<-customParms[1]
	levl<-customParms[2]
	levh<-customParms[3]
	led50<-parms[,1]
	lambda<-parms[,2]
	emax<-parms[,3]
	e0<-parms[,4]
	predval<-emaxfun(dtar,cbind(led50,lambda,emax,e0))
	
	lb<-quantile(predval,levl)
	ub<-quantile(predval,levh)
	return(c(predval=median(predval),lb=lb,ub=ub))
}


D4 <- emaxsimB(nsim,gen.parm,prior,modType=4,
							 mcmc=mcmc,binary=TRUE,seed=12357,
							 customCode=customCode,
							 customParms = c(0.2,.05,.95),nproc=nprocdef)

###
test_that("check DR for coverage",{
	expect_that(as.numeric(apply(abs(D4$fitpredv-D4$predpop)/D4$sepredv<1.96,2,mean)),
							equals(rep(0.95,ndose),tol=2*sqrt(.95*.05/nsim)))
})

lb<-numeric(nsim)
ub<-numeric(nsim)
for(i in 1:nsim){
	lb[i]<-plogis(D4$customOut[[i]][2])
	ub[i]<-plogis(D4$customOut[[i]][3])
}

mean(lb<=plogis(emaxfun(.2,pop)) & plogis(emaxfun(.2,pop))<=ub)

test_that("check DR for coverage from custom code",{
	expect_that(mean(lb<=plogis(emaxfun(.2,pop)) & plogis(emaxfun(.2,pop))<=ub),
							equals(0.90,tol=4*sqrt(.90*.10/nsim)))
})


####################
### gof tests
###
set.seed(12357)
nsim<-25
idmax<-5
doselev<-c(0,5,25,50,100)
n<-c(78,81,81,81,77)
Ndose<-length(doselev)

### population parameters for simulation
e0<-2.465375 
ed50<-67.481113 
emax<-4.127726
sdy<-1.0
pop<-c(log(ed50),emax,e0)    
meanlev<-emaxfun(doselev,pop)  
meanlev[5]<-meanlev[5]+1.0

###FixedMean is specialized constructor function for emaxsim
gen<-FixedMean(n,doselev,meanlev,sdy)  

prior<-emaxPrior.control(epmu=0,epsca=30,difTargetmu=0,difTargetsca=30,dTarget=100,p50=50,sigmalow=0.1,
										 sigmaup=30,parmDF=5)
mcmc<-mcmc.control(chains=1,warmup=500,iter=5000,seed=53453,
									 propInit=0.15,adapt_delta = 0.95)

D1 <- emaxsimB(nsim,gen, prior, modType=3,seed=12357,
							 mcmc=mcmc,check=FALSE,nproc=nprocdef)

test_that("gofP is very high",{
	expect(all(D1$gofP>0.40),'gofP not large as expected')
})

set.seed(12357)
nsim<-25
idmax<-5
doselev<-c(0,5,25,50,100)
n<-c(78,81,81,81,77)
Ndose<-length(doselev)

### population parameters for simulation
e0<-2.465375 
ed50<-67.481113 
emax<-4.127726
sdy<-1.0
pop<-c(log(ed50),emax,e0)    
meanlev<-emaxfun(doselev,pop)  
meanlev[5]<-meanlev[1]   ### back to pbo

###FixedMean is specialized constructor function for emaxsim
gen<-FixedMean(n,doselev,meanlev,sdy)  

prior<-emaxPrior.control(epmu=0,epsca=30,difTargetmu=0,difTargetsca=30,dTarget=100,p50=50,sigmalow=0.1,
										 sigmaup=30,parmDF=5)
mcmc<-mcmc.control(chains=1,warmup=500,iter=5000,seed=53453,
									 propInit=0.15,adapt_delta = 0.95)

D1low <- emaxsimB(nsim,gen, prior, modType=3,seed=12357,
									mcmc=mcmc,check=FALSE,nproc=nprocdef)

test_that("gofP is very low",{
	expect(all(D1low$gofP<0.05),'gofP not low as expected')
})

### binary
set.seed(12357)
nsim<-25

doselev<-c(0,0.033,2*0.033,4*0.033,8*0.033,16*0.033)
ndose<-length(doselev)
n<-10*rep(95,ndose)
###########################################
#### use pop emax model 
e0<-qlogis(.2)
ed50<-0.15
lambda<-1
emax<-qlogis(0.55)-e0
pop<-c(log(ed50),lambda,emax,e0)
meanlev<-plogis(emaxfun(doselev,parm=pop))
meanlev[6]<-0.55    ## higher than expected

gen.parm<-FixedMean(n,doselev,meanlev,sdy,parm=pop,binary=TRUE)  

prior<-emaxPrior.control(qlogis(0.2),4,0,4,16*0.033,0.05,parmDF=5,binary=TRUE)
mcmc<-mcmc.control(chains=1,warmup=500,iter=5000,seed=53453,propInit=0.15,adapt_delta = .95)

D1b <- emaxsimB(nsim,gen.parm,prior,modType=4,
							 mcmc=mcmc,binary=TRUE,seed=12357,nproc=nprocdef)

test_that("gofP is very high: binary",{
	expect(all(D1b$gofP>0.50),'gofP not large as expected')
})

####
set.seed(12357)
nsim<-25

doselev<-c(0,0.033,2*0.033,4*0.033,8*0.033,16*0.033)
ndose<-length(doselev)
n<-10*rep(95,ndose)
###########################################
#### use pop emax model 
e0<-qlogis(.2)
ed50<-0.15
lambda<-1
emax<-qlogis(0.55)-e0
pop<-c(log(ed50),lambda,emax,e0)
meanlev<-plogis(emaxfun(doselev,parm=pop))
meanlev[6]<-meanlev[1]    ## descend to pbo

gen.parm<-FixedMean(n,doselev,meanlev,sdy,parm=pop,binary=TRUE)  

prior<-emaxPrior.control(qlogis(0.2),4,0,4,16*0.033,0.05,parmDF=5,binary=TRUE)
mcmc<-mcmc.control(chains=1,warmup=500,iter=5000,seed=53453,propInit=0.15,adapt_delta = .95)

D1blow <- emaxsimB(nsim,gen.parm,prior,modType=4,
							 mcmc=mcmc,binary=TRUE,seed=12357,nproc=nprocdef)

test_that("gofP is very high: binary",{
	expect(all(D1blow$gofP<0.05),'gofP not small as expected')
})

####
set.seed(1244444)
nsim<-8

doselev<-c(0,0.033,2*0.033,4*0.033,8*0.033,16*0.033)
ndose<-length(doselev)
n<-10*rep(95,ndose)
###########################################
#### use pop emax model 
e0<-qlogis(.2)
ed50<-0.15
lambda<-1
emax<-qlogis(0.55)-e0
pop<-c(log(ed50),lambda,emax,e0)
meanlev<-plogis(emaxfun(doselev,parm=pop))

gen.parm<-FixedMean(n,doselev,meanlev,sdy,parm=pop,binary=TRUE)  

prior<-emaxPrior.control(qlogis(0.2),4,0,4,16*0.033,0.05,parmDF=5,binary=TRUE)
mcmc<-mcmc.control(chains=1,warmup=500,iter=5000,seed=53453,propInit=0.15,adapt_delta = .95)

D2 <- emaxsimB(nsim,gen.parm,prior,modType=4,
							 mcmc=mcmc,binary=TRUE,seed=12357,nproc=1)

out8<-D2[8]

fit8<-fitEmaxB(out8$y,out8$dose,prior=prior, modType = 4, binary=TRUE, diagnostics = FALSE)
parms<-coef(fit8)
gofP<-bpchkMonoEmax(fit8)

test_that("emaxsimB gofP equals bpchkMonoEmax return value ",{
	expect_equal(gofP,D2$gofP[8],tol=0.025,scale=1)
})

################################
##### test with bayesian random parameters
#####

set.seed(12357)

n<-c(99,95,98,94,98,98)
doselev<-c(0,5,10,25,50,150)

dose<-rep(doselev,n)


mcmc<-mcmc.control(chains=1,warmup=500,iter=5000,seed=53453,propInit=0.15,adapt_delta = .95)


priorbay<-emaxPrior.control(epmu=0,epsca=1,difTargetmu=1,difTargetsca=0.5,
												dTarget=150,p50=25,effDF=100,parmDF=100,
												 sigmalow=1.2,sigmaup=1.8)

genparmbay<-randomEmax(priorbay,n,doselev, modType='4' )



prior<-emaxPrior.control(epmu=0,epsca=4,difTargetmu=0,difTargetsca=4,
												 dTarget=150,p50=7.5,
												 sigmalow=0.01,sigmaup=3)

Dbay <- emaxsimB(nsim=1000,genparmbay,prior,modType=4,
							 mcmc=mcmc,binary=FALSE,seed=12357,nproc=nprocdef)

sink("NUL")
sout<-summary(Dbay)
sink()
test_that("emaxsimB random generated parms cov ",{
	expect_equal(as.numeric(sout$covMod),rep(0.9,length(doselev)-1),
							 tol=0.02,scale=1)
})
test_that("emaxsimB random generated parms mse",{
	expect_true(all(sout$mseMod<sout$mseMean))	
})

#### repeat with effective emax parameterization
priorE<-emaxPrior.control(epmu=0,epsca=4,difTargetmu=0,difTargetsca=4,
												 dTarget=2000,p50=7.5,
												 sigmalow=0.01,sigmaup=3)

DbayE <- emaxsimB(nsim=1000,genparmbay,priorE,modType=4,
							 mcmc=mcmc,binary=FALSE,seed=12357,nproc=nprocdef)
sink("NUL")
soutE<-summary(DbayE)
sink()
test_that("emaxsimB random generated parms cov 2 ",{
	expect_equal(as.numeric(soutE$covMod),rep(0.9,length(doselev)-1),
							 tol=0.02,scale=1)
})
test_that("emaxsimB random generated parms mse",{
	expect_true(all(soutE$mseMod<soutE$mseMean))	
})

################################
##### test with bayesian random parameters
##### binary data

set.seed(12357)

n<-c(99,95,98,94,98,98)
doselev<-c(0,5,10,25,50,150)

mcmc<-mcmc.control(chains=1,warmup=500,iter=5000,seed=53453,propInit=0.15,adapt_delta = .95)


priorbay<-emaxPrior.control(epmu=qlogis(.25),epsca=1,
												difTargetmu=qlogis(.75)-qlogis(.25),difTargetsca=1,
												dTarget=150,p50=50,effDF=100,parmDF=100,binary=TRUE)

genparmbayb<-randomEmax(priorbay,n,doselev, modType='4' )

priorb<-emaxPrior.control(epmu=0,epsca=4,difTargetmu=0,difTargetsca=4,
												 dTarget=150,p50=7.5,binary=TRUE)

Dbayb <- emaxsimB(nsim=1000,genparmbayb,priorb,modType=4,
							 mcmc=mcmc,seed=12357,binary=TRUE,nproc=nprocdef)

sink("NUL")
soutb<-summary(Dbayb)
sink()
test_that("emaxsimB random generated parms cov 2 ",{
	expect_equal(as.numeric(soutb$covMod),rep(0.9,length(doselev)-1),
							 tol=0.02,scale=1)
})
test_that("emaxsimB random generated parms mse",{
	expect_true(all(soutb$mseMod<soutb$mseMean))	
})
