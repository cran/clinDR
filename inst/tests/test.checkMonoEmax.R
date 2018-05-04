context('checkMonoEmax')


### continuous with aggregate data

data("examples14")
exdat<-examples14[[6]]

prior<-prior.control(epmu=0,epsd=10,emaxmu=0,emaxsd=10,p50=0.25,
										 sigmalow=0.01,sigmaup=3)
mcmc<-mcmc.control(chains=3)

fitout<-fitEmaxB(exdat$y,exdat$dose,prior,modType=4,prot=exdat$prot,
								 count=exdat$nsize,msSat=(exdat$sd)^2,mcmc=mcmc)

parm<-coef(fitout)
sigsim<-sigma(fitout)

p1<-checkMonoEmax(exdat$y,exdat$dose,parm,(sigsim)^2,
									nvec=exdat$nsize,trend='negative')

test_that("continuous model fit is good",{
  expect_gt(p1,0.5)
})

### make high group decline (do not rerun mcmc for altered fit)
ynm<-ifelse(exdat$dose==1,exdat$y+1.5,exdat$y)
set.seed(12357)
p1nm<-checkMonoEmax(ynm,exdat$dose,parm,(sigsim)^2,trend='negative')

test_that("coninuous model fit is bad",{
  expect_lt(p1nm,0.05)
})

##########
#### repeat with individual patient data
set.seed(12357)
popparm<-apply(parm[,1:4],2,median)
popmean<-emaxfun(exdat$dose,popparm)
popmean[6]<-popmean[6]+1.5*exdat$sd   ### ensure positive
dose<-rep(exdat$dose,exdat$nsize)
popmean<-rep(popmean,exdat$nsize)
y<-rnorm(length(popmean),popmean,exdat$sd)


prior<-prior.control(epmu=0,epsd=10,emaxmu=0,emaxsd=10,p50=0.25,
										 sigmalow=0.01,sigmaup=3)
mcmc<-mcmc.control(chains=3)

fitouti<-fitEmaxB(y,dose,prior,modType=4, mcmc=mcmc)

parmi<-coef(fitouti)
sigsimi<-sigma(fitouti)

p1v<-checkMonoEmax(y,dose,parmi,(sigsimi)^2,
									trend='negative')

test_that("individual data continuous model fit is good",{
  expect_gt(p1v,0.5)
})


ynm<-ifelse(dose==1,y+0.75,y)

p1nv<-checkMonoEmax(ynm,dose,parm,(sigsim)^2,
									trend='negative')

test_that("individual data coninuous model fit is bad",{
  expect_lt(p1nv,0.05)
})


##########################################################################
### convert continuous to binary for testing
set.seed(12357)
exdat<-examples14[[8]]
prior<-prior.control(epmu=0,epsd=4,emaxmu=0,emaxsd=4,p50=.15,binary=TRUE)
mcmc<-mcmc.control(chains=3)

fitoutb<-fitEmaxB(exdat$y,exdat$dose,prior,modType=4,count=exdat$nsize,
								 mcmc=mcmc,binary=TRUE)

parms<-coef(fitoutb)
popparms<-apply(parms,2,median)

dose<-rep(exdat$dose,exdat$nsize)
modp<-plogis(emaxfun(dose,popparms))
modp[dose==0.4]<-modp[dose==0.4]+.05

yb<-rbinom(length(modp),1,modp)

fitoutb2<-fitEmaxB(yb,dose,prior,modType=4,
								 mcmc=mcmc,binary=TRUE)
parmb<-coef(fitoutb2)

set.seed(12357)
p1b<-checkMonoEmax(yb,dose,parmb,rep(1,nrow(parmb)),trend='negative',logit=TRUE)

test_that("binary model fit is good",{
  expect_gt(p1b,0.2)
})

## non-monotone 
ybnm<-ifelse(dose==1,rbinom(sum(dose==1),1,.3),yb)

set.seed(12357)
p1bnm<-checkMonoEmax(ybnm,dose,parmb,rep(1,nrow(parmb)),trend='negative',logit=TRUE)

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
p1bv<-checkMonoEmax(ybv,dvec,parmb,rep(1,nrow(parmb)),nvec=nvec,trend='negative',logit=TRUE)

test_that("aggregate good binary model check",{
  expect_equal(p1b,p1bv)
})

set.seed(12357)
p1bnmv<-checkMonoEmax(ybnmv,dvec,parmb,rep(1,nrow(parmb)),nvec=nvec,trend='negative',logit=TRUE)

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
p1bv0<-checkMonoEmax(yc,dvec,parmb,rep(1,nrow(parmb)),nvec=nvec,trend='negative',logit=TRUE)

test_that("aggregate good coninuous model check alternative form",{
  expect_equal(p1bv,p1bv0)
})


