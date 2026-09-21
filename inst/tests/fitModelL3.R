setwd("C:/OneDrive/OneDrive - Pfizer/consulting/alopecia/missingData")
library(rstan)
### import model fit and MI functions
source('sourcev2.R')

### read simulated data
testnum<-2
fname<-paste0('sim1015data',testnum,'_longm.csv')
datin<-read.csv(fname)

dat<-datin
dat$trt<-dat$trt+1
nvisit<-max(dat$time)


##########################################################
####  fit bayes longitudinal model
####
mcmc<-list(chains=3,thin=4,warmup=2000,iter=4*3333+2000,
					 adapt_delta=0.95,seed=12357)

outres<-fitModel(id=dat$id,y=dat$resp,trt=dat$trt,visit=dat$time,
				 prmean0=rep(-1.7,nvisit),prsd0=rep(3.0,nvisit),
				 prmean=rep(0,nvisit),prsd=rep(3.0,nvisit),
				 gparm1=2,gparm2=1,
				 mcmc=mcmc)

### extract model output
outmod<-outres$stanfit
beta<-outres$beta
sigma<-outres$sigma
theta<-outres$theta

##########################################################
### convergence diagnostics
traceplot(outmod,pars=c('beta'))
traceplot(outmod,pars=c('sigma'))
traceplot(outmod,pars=c('theta[1]','theta[2]','theta[3]','theta[4]'))
pset<-paste0('theta[',(max(dat$id)-3):max(dat$id),']')
traceplot(outmod,pars=pset)

stan_rhat(outmod,pars=c('beta'))
stan_rhat(outmod,pars='theta')
Rhat(as.matrix(outmod,pars='sigma'))
stan_ac(outmod,pars=c('beta'))
stan_ac(outmod,pars=c('sigma'))


stan_mcse(outmod)
stan_ess(outmod)

save.image(paste0('simevalL2',testnum,'.RData'))
##############################################################


##########################################################
##########################################################
### MI analysis with sensitivity parameter

### format data to selected visit and treatments with NA included
mdat<-inputmi(id=dat$id,trt=dat$trt,y=dat$resp,visit=dat$time,trtsel=1:2,vsel=max(dat$time))


### create matrix with imputation probabilities
mprobs<-miprobs(mdat=mdat,vsel=5,beta=beta,sigma=sigma,nimp=100)

### MAR
set.seed(12357)
miout<-midat(mprobs=mprobs,trt=mdat$trt,m=mdat$m,deltat=0,deltac=0)
miout

### MAR (return imputed responses)
set.seed(12357)
miout<-midat(mprobs=mprobs,trt=mdat$trt,m=mdat$m,deltat=0,deltac=0,returnYimp=TRUE)
miout$miest
miout$mise
miout$midf
dim(miout$yimp)

### MAR (use normal approximation to binomial)
set.seed(12357)
miout<-midat(mprobs=mprobs,trt=mdat$trt,m=mdat$m,deltat=0,deltac=0,f=binNorm)
miout



### tipping point examples
set.seed(12357)
miout<-midat(mprobs=mprobs,trt=mdat$trt,m=mdat$m,
						 deltac=0,deltat=0.2)
miout
	
miout<-midat(mprobs=mprobs,trt=mdat$trt,m=mdat$m,
						 deltac=-0.1,deltat=0.2)
miout
	



###################################################################
####  assign non-covid to NR rather than impute
####
#### create a simulated covid variable
####
misID<-unique(dat$id[is.na(dat$resp)])
nmis<-length(misID)
covidmisID<-misID[c(1:4,(nmis-3):nmis)]
covid19<-rep(0,nrow(dat))
covid19[dat$id %in% covidmisID]<-1
datcovid<-dat
datcovid$covid19<-covid19



### assign NR to NA not covid19
datcovid$resp[is.na(datcovid$resp)& datcovid$covid19==0]<-0
datcovid[dat$id%in%misID,]

mcovdat<-inputmi(id=datcovid$id,trt=datcovid$trt,y=datcovid$resp,
								 visit=datcovid$time,trtsel=1:2,vsel=max(datcovid$time))
mcovprobs<-miprobs(mdat=mcovdat,vsel=5,beta=beta,sigma=sigma,nimp=100)

### MAR
set.seed(12357)
micovout<-midat(mprobs=mcovprobs,trt=mcovdat$trt,m=mcovdat$m,deltat=0,deltac=0)
micovout

### tipping point examples
set.seed(12357)
micovout<-midat(mprobs=mcovprobs,trt=mcovdat$trt,m=mcovdat$m,
						 deltac=0,deltat=0.20)
micovout
	
set.seed(12357)
micovout<-midat(mprobs=mcovprobs,trt=mcovdat$trt,m=mcovdat$m,
						 deltac=-0.1,deltat=0.2)
micovout

set.seed(12357)
micovout<-midat(mprobs=mcovprobs,trt=mcovdat$trt,m=mcovdat$m,
						 deltac=-0.999,deltat=0.999)
micovout
	

save.image(paste0('simevalL',testnum,'.RData'))
	
