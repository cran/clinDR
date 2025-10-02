context('checkMonoEmax local parameters')


### continuous with aggregate data

data("metaData")
exdat<-metaData[metaData$taid==6 & metaData$poptype==1,]

prior<-emaxPrior.control(epmu=0,epsca=10,difTargetmu=0,difTargetsca=10,dTarget=80.0,
        p50=3.75,sigmalow=0.01,sigmaup=20)
mcmc<-mcmc.control(chains=3)

msSat<-sum((exdat$sampsize-1)*(exdat$sd)^2)/(sum(exdat$sampsize)-length(exdat$sampsize))
fitout<-fitEmaxB(exdat$rslt,exdat$dose,prior,modType=4,
				count=exdat$sampsize,msSat=msSat,mcmc=mcmc,nproc=3)
parm<-coef(fitout)[,1:4]  #use first intercept
sigsim<-sigma(fitout)

p1<-suppressWarnings(checkMonoEmax(exdat$rslt,exdat$dose,parm,(sigsim)^2,
									nvec=exdat$sampsize,trend='negative'))

test_that("continuous model fit is good",{
  expect_gt(p1,0.5)
})

### make high group decline (do not rerun mcmc for altered fit)
ynm<-ifelse(exdat$dose==80,exdat$rslt+40,exdat$rslt)
set.seed(12357)
p1nm<-suppressWarnings(checkMonoEmax(ynm,exdat$dose,parm,
													(sigsim)^2,trend='negative'))

test_that("coninuous model fit is bad",{
  expect_lt(p1nm,0.05)
})

##########
#### repeat with individual patient data
set.seed(12357)
popparm<-apply(parm[,1:4],2,median)
popmean<-emaxfun(exdat$dose,popparm)
popmean[6]<-popmean[6]
dose<-rep(exdat$dose,exdat$sampsize)
popmean<-rep(popmean,exdat$sampsize)
y<-rnorm(length(popmean),popmean,exdat$sd)

prior<-emaxPrior.control(epmu=0,epsca=10,difTargetmu=0,difTargetsca=10,
												 dTarget=80,p50=3.75,
										 sigmalow=0.01,sigmaup=3)
mcmc<-mcmc.control(chains=3)

fitouti<-fitEmaxB(y,dose,prior,modType=4, mcmc=mcmc,nproc=3)

parmi<-coef(fitouti)
sigsimi<-sigma(fitouti)

p1v<-suppressWarnings(checkMonoEmax(y,dose,parmi,(sigsimi)^2,
									trend='negative'))

test_that("individual data continuous model fit is good",{
  expect_gt(p1v,0.5)
})


ynm<-ifelse(dose==80,y+40,y)

p1nv<-suppressWarnings(checkMonoEmax(ynm,dose,parm,(sigsim)^2,
									trend='negative'))

test_that("individual data coninuous model fit is bad",{
  expect_lt(p1nv,0.05)
})


##########################################################################
### convert continuous to binary for testing
set.seed(12357)
data(metaData)
exdat<-metaData[metaData$taid==8,]

cy<-round(exdat$sampsize*exdat$rslt)
y<-c(rep(1,length(cy)),rep(0,length(cy)))
cy<-c(cy,exdat$sampsize-cy)
drep<-c(exdat$dose,exdat$dose)


prior<-suppressWarnings(prior.control(epmu=0,epsd=4,emaxmu=0,emaxsd=4,p50=.15,binary=TRUE))
mcmc<-mcmc.control(chains=3)

fitoutb<-fitEmaxB(y,drep,prior,modType=4,count=cy,
								 mcmc=mcmc,binary=TRUE,nproc=3)

parms<-coef(fitoutb)
popparms<-apply(parms,2,median)

dose<-rep(drep,cy)
modp<-plogis(emaxfun(dose,popparms))
modp[dose==4]<-modp[dose==4]+.10
modp[dose==2.5]<-modp[dose==2.5]+.05

yb<-rbinom(length(modp),1,modp)

fitoutb2<-fitEmaxB(yb,dose,prior,modType=4,
								 mcmc=mcmc,binary=TRUE,nproc=3)
parmb<-coef(fitoutb2)

set.seed(12357)
p1b<-suppressWarnings(checkMonoEmax(yb,dose,parmb,
				rep(1,nrow(parmb)),trend='negative',binary=TRUE))

test_that("binary model fit is good",{
  expect_gt(p1b,0.2)
})

## non-monotone 
ybnm<-ifelse(dose==10,rbinom(sum(dose==10),1,.3),yb)

set.seed(12357)
p1bnm<-suppressWarnings(checkMonoEmax(ybnm,dose,parmb,
						rep(1,nrow(parmb)),trend='negative',binary=TRUE))

test_that("binary model fit is bad",{
  expect_lt(p1bnm,0.05)
})

###########
#### repeat with aggregate data
#### first with weighted averages of means
dvec<-sort(unique(dose))

ybv<-tapply(yb,dose,mean)
ybnmv<-tapply(ybnm,dose,mean)
nvec<-as.vector(table(dose))

set.seed(12357)
p1bv<-suppressWarnings(checkMonoEmax(ybv,dvec,parmb,
			rep(1,nrow(parmb)),nvec=nvec,trend='negative',binary=TRUE))

test_that("aggregate good binary model check",{
  expect_equal(p1b,p1bv)
})

set.seed(12357)
p1bnmv<-suppressWarnings(checkMonoEmax(ybnmv,dvec,parmb,
					rep(1,nrow(parmb)),nvec=nvec,trend='negative',binary=TRUE))

test_that("aggregate bad coninuous model fit",{
  expect_equal(p1bnm,p1bnmv)
})

### repeat with y as 0/1 counts

yc<-c(rep(1,7),rep(0,7))
dvec<-c(dvec,dvec)
nvec<-rep(0,14)
for(i in 1:7){
	nvec[i]<-sum(yb==1 & dose==dvec[i])	
	nvec[i+7]<-sum(yb==0 & dose==dvec[i])	
}
set.seed(12357)
p1bv0<-suppressWarnings(checkMonoEmax(yc,dvec,parmb,
			rep(1,nrow(parmb)),nvec=nvec,trend='negative',binary=TRUE))

test_that("aggregate good coninuous model check alternative form",{
  expect_equal(p1bv,p1bv0)
})

############
### with covariates

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

testout2<-fitEmaxB(y,dose,prior=prior,modType=4,prot=prots,xbase=x,
									 mcmc=mcmc,diagnostics=FALSE,nproc=3)

parms<-coef(testout2)
parms<-parms[,c(1:4,6:8)]
sigsim<-sigma(testout2)

y1<-y[prots==1]
dose1<-dose[prots==1]

p1<-suppressWarnings(checkMonoEmax(y1,dose1,parms,(sigsim)^2,xbase=x1, 
															trend='positive'))

test_that("continuous model fit with covariate is good",{
  expect_gt(p1,0.5)
})

### check lack of fit
y1nm<-y1
y1nm[dose1==350]<-rnorm(sum(dose1==350),5,sdy)

p1nm<-suppressWarnings(checkMonoEmax(y1nm,dose1,parms,(sigsim)^2,
										xbase=x1, trend='positive'))


test_that("cont model fit with covariate is bad",{
  expect_lt(p1nm,0.05)
})

### binary, covariates
### check with larger n for better asymptotics
### 4-parm model, 2 covariates, 1 protocol

set.seed(20572)

doselev<-c(0,5,25,50,100,350)
n<-2*c(78,81,81,81,77,80)
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

testout5b<-fitEmaxB(y,dose,prot=prot,prior=prior,modType=4,xbase=x,
										mcmc=mcmc,diagnostics=FALSE,binary=TRUE,nproc=1)

parms<-coef(testout5b)
parms<-parms[,c(1:4,6:7)]


y1<-y[prot==1]
dose1<-dose[prot==1]
p1b<-suppressWarnings(checkMonoEmax(y1,dose1,parms,sigma2=
									numeric(length(y1)),xbase=x1,binary=TRUE))

test_that("binary model fit with covariate is good",{
  expect_gt(p1,0.25)
})

### check lack of fit
y1nm<-y1
y1nm[dose1==350]<-rbinom(sum(dose1==350),1,0.2)

p1bnm<-suppressWarnings(checkMonoEmax(y1nm,dose1,parms,(sigsim)^2,
													xbase=x1, trend='positive',binary=TRUE))

test_that("binary model fit with covariate is bad",{
  expect_lt(p1bnm,0.05)
})


###################################################################
###################################################################
#### update changed to testing bpchkMonEmax
###
### continuous with aggregate data

data("metaData")
exdat<-metaData[metaData$taid==6 & metaData$poptype==1,]

prior<-emaxPrior.control(epmu=0,epsca=10,difTargetmu=0,difTargetsca=10,dTarget=80.0,
        p50=3.75,sigmalow=0.01,sigmaup=20)
mcmc<-mcmc.control(chains=3)

msSat<-sum((exdat$sampsize-1)*(exdat$sd)^2)/(sum(exdat$sampsize)-length(exdat$sampsize))
fitout<-fitEmaxB(exdat$rslt,exdat$dose,prior,modType=4,
				count=exdat$sampsize,msSat=msSat,mcmc=mcmc,nproc=3)

parm<-coef(fitout)[,1:4]  #use first intercept
sigsim<-sigma(fitout)

p1<-bpchkMonoEmax(fitout,trend='negative')

test_that("continuous model fit is good",{
  expect_gt(p1,0.5)
})

### make high group decline (do not rerun mcmc for altered fit)
fitoutmod<-fitout
fitoutmod$y<-ifelse(exdat$dose==80,exdat$rslt+40,exdat$rslt)
set.seed(12357)
p1nm<-bpchkMonoEmax(fitoutmod,trend='negative')

test_that("coninuous model fit is bad",{
  expect_lt(p1nm,0.05)
})

##########
#### repeat with individual patient data
set.seed(12357)
popparm<-apply(parm[,1:4],2,median)
popmean<-emaxfun(exdat$dose,popparm)
popmean[6]<-popmean[6]
dose<-rep(exdat$dose,exdat$sampsize)
popmean<-rep(popmean,exdat$sampsize)
y<-rnorm(length(popmean),popmean,exdat$sd)

prior<-emaxPrior.control(epmu=0,epsca=10,difTargetmu=0,difTargetsca=10,
												 dTarget=80,p50=3.75,
										 sigmalow=0.01,sigmaup=3)
mcmc<-mcmc.control(chains=3)

fitouti<-fitEmaxB(y,dose,prior,modType=4, mcmc=mcmc,nproc=3)

p1v<-bpchkMonoEmax(fitouti,trend='negative')

test_that("individual data continuous model fit is good",{
  expect_gt(p1v,0.5)
})


fitoutimod<-fitouti
fitoutimod$y<-ifelse(dose==80,y+40,y)

p1nv<-bpchkMonoEmax(fitoutimod,trend='negative')

test_that("individual data coninuous model fit is bad",{
  expect_lt(p1nv,0.05)
})


##########################################################################
### convert continuous to binary for testing
set.seed(12357)
data(metaData)
exdat<-metaData[metaData$taid==8,]

cy<-round(exdat$sampsize*exdat$rslt)
y<-c(rep(1,length(cy)),rep(0,length(cy)))
cy<-c(cy,exdat$sampsize-cy)
drep<-c(exdat$dose,exdat$dose)


prior<-suppressWarnings(prior.control(epmu=0,epsd=4,emaxmu=0,emaxsd=4,p50=.15,binary=TRUE))
mcmc<-mcmc.control(chains=3)

fitoutb<-fitEmaxB(y,drep,prior,modType=4,count=cy,
								 mcmc=mcmc,binary=TRUE,nproc=3)

parms<-coef(fitoutb)
popparms<-apply(parms,2,median)

dose<-rep(drep,cy)
modp<-plogis(emaxfun(dose,popparms))
modp[dose==4]<-modp[dose==4]+.10
modp[dose==2.5]<-modp[dose==2.5]+.05

yb<-rbinom(length(modp),1,modp)

fitoutb2<-fitEmaxB(yb,dose,prior,modType=4,
								 mcmc=mcmc,binary=TRUE,nproc=3)
parmb<-coef(fitoutb2)

set.seed(12357)
p1b<-bpchkMonoEmax(fitoutb2,trend='negative')

test_that("binary model fit is good",{
  expect_gt(p1b,0.2)
})

## non-monotone 
fitoutb2mod<-fitoutb2
fitoutb2mod$y<-ifelse(dose==10,rbinom(sum(dose==10),1,.3),yb)

set.seed(12357)
p1bnm<-bpchkMonoEmax(fitoutb2mod,trend='negative')

test_that("binary model fit is bad",{
  expect_lt(p1bnm,0.05)
})

###########
#### repeat with aggregate data
#### weighted averages of means
dvec<-sort(unique(dose))
yvec<-c(rep(1,length(dvec)),rep(0,length(dvec)))
dvec<-c(dvec,dvec)

n1<-tapply(yb,dose,sum)
ntot<-as.vector(table(dose))
nvec<-c(n1,ntot-n1)

fitoutb2<-fitEmaxB(yvec,dvec,prior,modType=4,count=nvec,
								 mcmc=mcmc,binary=TRUE,nproc=3)
set.seed(12357)
p1bv<-bpchkMonoEmax(fitoutb2,trend='negative')

test_that("aggregate good binary model check",{
  expect_equal(p1b,p1bv)
})

############
### with covariates
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

testout2<-fitEmaxB(y,dose,prior=prior,modType=4,prot=prots,xbase=x,
									 mcmc=mcmc,diagnostics=FALSE,nproc=3)

parms<-coef(testout2)

y1<-y[prots==1]
dose1<-dose[prots==1]

p1<-bpchkMonoEmax(testout2,trend='positive',protSel=1)
p2<-bpchkMonoEmax(testout2,trend='positive',protSel=2)

test_that("continuous model fit with covariate is good",{
  expect_gt(p1,0.5)
})

### check lack of fit
testoutmod2<-testout2
testoutmod2$y[testout2$dose==350] <-rnorm(sum(testout2$dose==350),5,sdy)
testoutmod2$y[testout2$dose==100] <-rnorm(sum(testout2$dose==100),5,sdy)
testoutmod2$y[testout2$dose==50] <-rnorm(sum(testout2$dose==50),3,sdy)

p1nm<-bpchkMonoEmax(testoutmod2, trend='positive',protSel=1)
p2nm<-bpchkMonoEmax(testoutmod2, trend='positive',protSel=2)

test_that("cont model fit with covariate is bad",{
  expect_lt(p1nm,0.05)
})
test_that("cont model fit with covariate is bad",{
  expect_lt(p2nm,0.05)
})

#### repeat with 3-parm fit
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

testout2<-fitEmaxB(y,dose,prior=prior,modType=3,prot=prots,xbase=x,
									 mcmc=mcmc,diagnostics=FALSE,nproc=3)

parms<-coef(testout2)

y1<-y[prots==1]
dose1<-dose[prots==1]

p1<-bpchkMonoEmax(testout2,trend='positive',protSel=1)
p2<-bpchkMonoEmax(testout2,trend='positive',protSel=2)

test_that("continuous model fit with covariate is good: 3-parm",{
  expect_gt(p1,0.5)
})
test_that("continuous model fit with covariate is good: 3-parm",{
  expect_gt(p2,0.5)
})

### check lack of fit
testoutmod2<-testout2
testoutmod2$y[testout2$dose==350] <-rnorm(sum(testout2$dose==350),5,sdy)
testoutmod2$y[testout2$dose==100] <-rnorm(sum(testout2$dose==100),5,sdy)
testoutmod2$y[testout2$dose==50] <-rnorm(sum(testout2$dose==50),3,sdy)

p1nm<-bpchkMonoEmax(testoutmod2, trend='positive',protSel=1)
p2nm<-bpchkMonoEmax(testoutmod2, trend='positive',protSel=2)

test_that("cont model fit with covariate is bad:3-parm",{
  expect_lt(p1nm,0.05)
})
test_that("cont model fit with covariate is bad:3-parm",{
  expect_lt(p2nm,0.05)
})

### no covariates
prior$basemu<-NULL
prior$basevar<-NULL
testout2<-fitEmaxB(y,dose,prior=prior,modType=3,prot=prots,
									 mcmc=mcmc,diagnostics=FALSE,nproc=3)

parms<-coef(testout2)

y1<-y[prots==1]
dose1<-dose[prots==1]

p1<-bpchkMonoEmax(testout2,trend='positive',protSel=1)
p2<-bpchkMonoEmax(testout2,trend='positive',protSel=2)

test_that("continuous model fit without covariate is good: 3-parm",{
  expect_gt(p1,0.5)
})
test_that("continuous model fit without covariate is good: 3-parm",{
  expect_gt(p2,0.5)
})

### check lack of fit
testoutmod2<-testout2
testoutmod2$y[testout2$dose==350] <-rnorm(sum(testout2$dose==350),5,sdy)
testoutmod2$y[testout2$dose==100] <-rnorm(sum(testout2$dose==100),5,sdy)
testoutmod2$y[testout2$dose==50] <-rnorm(sum(testout2$dose==50),3,sdy)

p1nm<-bpchkMonoEmax(testoutmod2, trend='positive',protSel=1)
p2nm<-bpchkMonoEmax(testoutmod2, trend='positive',protSel=2)

test_that("cont model fit without covariate is bad:3-parm",{
  expect_lt(p1nm,0.05)
})
test_that("cont model fit without covariate is bad:3-parm",{
  expect_lt(p2nm,0.05)
})

### repeat without covariates, 2 prot and 4-parm
prior$basemu<-NULL
prior$basevar<-NULL
testout2<-fitEmaxB(y,dose,prior=prior,modType=4,prot=prots,
									 mcmc=mcmc,diagnostics=FALSE,nproc=3)

parms<-coef(testout2)

y1<-y[prots==1]
dose1<-dose[prots==1]

p1<-bpchkMonoEmax(testout2,trend='positive',protSel=1)
p2<-bpchkMonoEmax(testout2,trend='positive',protSel=2)

test_that("continuous model fit without covariate, 2-prot is good: 4-parm",{
  expect_gt(p1,0.5)
})
test_that("continuous model fit without covariate, 2-prot is good: 4-parm",{
  expect_gt(p2,0.5)
})

### check lack of fit
testoutmod2<-testout2
testoutmod2$y[testout2$dose==350] <-rnorm(sum(testout2$dose==350),5,sdy)
testoutmod2$y[testout2$dose==100] <-rnorm(sum(testout2$dose==100),5,sdy)
testoutmod2$y[testout2$dose==50] <-rnorm(sum(testout2$dose==50),3,sdy)

p1nm<-bpchkMonoEmax(testoutmod2, trend='positive',protSel=1)
p2nm<-bpchkMonoEmax(testoutmod2, trend='positive',protSel=2)

test_that("cont model fit without covariate,2-prot is bad:4-parm",{
  expect_lt(p1nm,0.05)
})
test_that("cont model fit without covariate,2-prot is bad:4-parm",{
  expect_lt(p2nm,0.05)
})

### binary, covariates
### check with larger n for better asymptotics
### 4-parm model, 2 covariates, 1 protocol

set.seed(20572)

doselev<-c(0,5,25,50,100,350)
n<-2*c(78,81,81,81,77,80)
ntot<-sum(n)

### population parameters for simulation
e0<- -1.5 
ed50<-67.481113 
emax<-4.0
pop<-c(log(ed50),emax,e0)    
dose<-c(rep(doselev,n/2),rep(doselev,n/2))
prot<-sort(rep(1,ntot))   #1 prot for bp check
meanlev<-emaxfun(dose,pop)  

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

testout5b<-fitEmaxB(y,dose,prot=prot,prior=prior,modType=4,xbase=x,
										mcmc=mcmc,diagnostics=FALSE,binary=TRUE,nproc=1)

p1b<-bpchkMonoEmax(testout5b)

test_that("binary model fit with covariate is good",{
  expect_gt(p1,0.25)
})

### check lack of fit
testout5bmod<-testout5b
testout5bmod$y[testout5b$dose==350]<-rbinom(sum(testout5b$dose==350),1,0.2)

p1bnm<-bpchkMonoEmax(testout5bmod, trend='positive')

test_that("binary model fit with covariate is bad",{
  expect_lt(p1bnm,0.05)
})
