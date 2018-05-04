context('fitEmax')

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


testout<-fitEmax(y,dose,modType=4,prot=prots,diagnostics=FALSE)

se<-sqrt(diag(vcov(testout)))
z<-(coef(testout)-c(pop[1],1,pop[2:3],pop[3]))/se

### check parameter estimates
test_that("model parameters agree within 2se",{
	expect_lt(as.numeric(max(abs(z))),2.0)
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
#### repeat with 3 parm, int=2, starting value given
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


testout<-fitEmax(ysub,dsub,iparm=pop,modType=3,prot=protsub,pboAdj=TRUE,
								 diagnostics=FALSE)

se<-sqrt(diag(vcov(testout)))
z<-(coef(testout)-pop[1:2])/se

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
#### no intercept with 4 parm, int=1, starting value given
#### non-zero reference dose

testout<-fitEmax(ysub,dsub,iparm=c(pop[1],1,pop[2],0),
								 modType=4,prot=protsub,pboAdj=TRUE,diagnostics=FALSE)

se<-sqrt(diag(vcov(testout)))
z<-(coef(testout)-c(pop[1],1,pop[2]))/se

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



###  simulation to test gof test (primarily)
###  also test predictions

set.seed(12357)
nsim<-1000
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
sdy<-1
pop<-c(log(ed50),emax,e0)    
dose<-rep(doselev,n)
meanlev<-emaxfun(dose,pop)  
poppred<-emaxfun(c(20,80),pop)
popref<-emaxfun(50,pop)

modtype<-4
if(modtype==4){pparm<-c(pop[1],1,pop[2:3],pop[3])
}else pparm<-pop
z<-matrix(numeric((modtype+1)*nsim),ncol=modtype+1)
gof<-numeric(nsim)
zabs<-matrix(numeric(nsim*2),ncol=2)
zdif<-matrix(numeric(nsim*2),ncol=2)
for(i in 1:nsim){
	y<-rnorm(n1+n2,meanlev,sdy)
	prots<-c(rep('d',n1),rep('c',n2))
	testout<-fitEmax(y,dose,modType=modtype,prot=prots,diagnostics=FALSE)
	se<-sqrt(diag(vcov(testout)))
	z[i,]<-(coef(testout)-pparm)/se
	gof[i]<-testout$gofTest
	predout<-predict(testout,dosevec=c(20,80),int=1,dref=50)
	zabs[i,]<-(predout$pred-poppred)/predout$se
	zdif[i,]<-(predout$fitdif-(poppred-popref))/predout$sedif
}
### check parameter estimates
test_that("model parameters agree within 2se",{
	expect_that(0.05,
							equals(as.numeric(max(apply(abs(z)>1.96,2,mean))),
										 tolerance=2*sqrt(.05*.95/nsim),scale=1))
})
test_that("GOF test statistic",{
	expect_that(0.05,
							equals(as.numeric(mean(gof<0.05)),tolerance=2*sqrt(.05*.95/nsim),scale=1))
})
test_that("predictions agree within 2se",{
	expect_that(0.05,
							equals(as.numeric(max(apply(abs(zabs)>1.96,2,mean))),
										 tolerance=2*sqrt(.05*.95/nsim),scale=1))
})
test_that("dif predictions agree within 2se",{
	expect_that(0.05,
							equals(as.numeric(max(apply(abs(zdif)>1.96,2,mean))),
										 tolerance=2*sqrt(.05*.95/nsim),scale=1))
})




### repeat with 3 parm model and grouped data, one protocol
### and replicated measurement per dose/protocol condition
set.seed(12357)
modtype<-3
if(modtype==4){pparm<-c(pop[1],1,pop[2:3])
}else pparm<-pop
z<-matrix(numeric(modtype*nsim),ncol=modtype)
gof<-numeric(nsim)
zabs<-matrix(numeric(nsim*2),ncol=2)
zdif<-matrix(numeric(nsim*2),ncol=2)
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
	testout<-fitEmax(ysum,doselev,count=n,modType=modtype,msSat=msSat,diagnostics=FALSE)
	se<-sqrt(diag(vcov(testout)))
	z[i,]<-(coef(testout)-pparm)/se
	gof[i]<-testout$gofTest
	predout<-predict(testout,dosevec=c(20,80),int=1,dref=50)
	zabs[i,]<-(predout$pred-poppred)/predout$se
	zdif[i,]<-(predout$fitdif-(poppred-popref))/predout$sedif
}
### check parameter estimates
test_that("grouped data model parameters agree within 2se",{
	expect_that(0.05,
							equals(as.numeric(max(apply(abs(z)>1.96,2,mean))),
										 tolerance=2*sqrt(.05*.95/nsim),scale=1))
})
test_that("grouped data GOF test statistic",{
	expect_that(0.05,
							equals(as.numeric(mean(gof<0.05)),tolerance=2*sqrt(.05*.95/nsim),scale=1))
})

test_that("predictions agree within 3se",{
	expect_that(0.05,
							equals(as.numeric(max(apply(abs(zabs)>1.96,2,mean))),
										 tolerance=3*sqrt(.05*.95/nsim),scale=1))
})
test_that("dif predictions agree within 2se",{
	expect_that(0.05,
							equals(as.numeric(max(apply(abs(zdif)>1.96,2,mean))),
										 tolerance=2*sqrt(.05*.95/nsim),scale=1))
})





#########################################################################
##### binary
set.seed(12357)

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

### fitEmax inputs
y<-c(rep(1,nd1),rep(0,nd1),rep(1,nd2),rep(0,nd2))
counts<-c(y1,n1-y1,y2,n2-y2)
prots<-c(rep(1,2*nd1),rep(2,2*nd2))
dvec<-c(dvec1,dvec1,dvec2,dvec2)


testout<-fitEmax(y,dvec,modType=4,parms[1:4],
								 prot=prots,
								 count=counts,binary=TRUE,diagnostics=FALSE)

se<-sqrt(diag(vcov(testout)))
z<-(coef(testout)-parms)/se

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

##############################################################
### repeat without starting values, 3 parameter model
parms<-c(log(0.25),1.4,-0.5,-0.85)

testout<-fitEmax(y,dvec,modType=3,
								 prot=prots,
								 count=counts,binary=TRUE,diagnostics=FALSE)

se<-sqrt(diag(vcov(testout)))
z<-(coef(testout)-parms)/se

### check parameter estimates
test_that("model parameters agree within 2se: no parms",{
	expect_lt(as.numeric(max(abs(z))),2.0)
})

### check predictions
predout<-predict(testout,dosevec=c(.25,.75),int=2)

poppred<-plogis(emaxfun(c(.25,.75),parms[c(1:2,4)]))
popref<-plogis(emaxfun(0,parms[c(1:2,4)]))
z<-(predout$pred-poppred)/predout$se
zdif<-(predout$fitdif-(poppred-popref))/predout$sedif

### check predictions 
test_that("predictions agree within 2se: no parms",{
	expect_lt(as.numeric(max(abs(z))),2.0)
	expect_lt(as.numeric(max(abs(zdif))),2.0)
})

test_that("check absolute levels: no parms",{
	expect_that(as.numeric(predout$pred),
							equals(poppred,tol=2*.5/sqrt(10000),scale=1))
})


#### simulations to check gof (primarily)
set.seed(12357)
nsim<-1000
dvec1<-c(0,.1,.3,.6,1)
dvec2<-c(0,.1,.2,.4,.6,1)
nd1<-length(dvec1)
nd2<-length(dvec2)
n1<-rep(10000,nd1)
n2<-rep(10000,nd2)
prots<-c(rep('a',2*nd1),rep('b',2*nd2))  # non-integer protocol identifiers
dvec<-c(dvec1,dvec1,dvec2,dvec2)

parms<-c(log(0.25),1,1.4,-0.5,-0.85)
modtype<-4
if(modtype==4){pparm<-parms
}else pparm<-parms[c(1,3:5)]
mlev1<-plogis(emaxfun(dvec1,pparm[1:4]))
mlev2<-plogis(emaxfun(dvec2,pparm[c(1:3,5)]))
poppred<-plogis(emaxfun(c(.2,.5,.8),parms[c(1:3,4)]))
popref<-plogis(emaxfun(0,parms[c(1:3,4)]))

z<-matrix(numeric((modtype+1)*nsim),ncol=modtype+1)
gof<-numeric(nsim)
zabs<-matrix(numeric(nsim*3),ncol=3)
zdif<-matrix(numeric(nsim*3),ncol=3)
for(i in 1:nsim){
	y1<-rbinom(nd1,n1,mlev1)
	y2<-rbinom(nd2,n2,mlev2)
	
	### fitEmax inputs
	y<-c(rep(1,nd1),rep(0,nd1),rep(1,nd2),rep(0,nd2))
	counts<-c(y1,n1-y1,y2,n2-y2)
	testout<-fitEmax(y,dvec,modType=modtype,parms[1:4],
									 prot=prots,
									 count=counts,binary=TRUE,diagnostics=FALSE)
	se<-sqrt(diag(vcov(testout)))
	z[i,]<-(coef(testout)-pparm)/se
	gof[i]<-testout$gofTest
	predout<-predict(testout,dosevec=c(.2,.5,.8),int=1,dref=0)
	zabs[i,]<-(predout$pred-poppred)/predout$se
	zdif[i,]<-(predout$fitdif-(poppred-popref))/predout$sedif
}
### check parameter estimates
test_that("binary model parameters agree within 2se",{
	expect_that(0.05,
							equals(as.numeric(max(apply(abs(z)>1.96,2,mean))),
										 tolerance=2*sqrt(.05*.95/nsim),scale=1))
})
test_that("binary GOF test statistic",{
	expect_that(0.05,
							equals(as.numeric(mean(gof<0.05)),tolerance=2*sqrt(.05*.95/nsim),scale=1))
})

test_that("binary predictions agree within 2se",{
	expect_that(0.05,
							equals(as.numeric(max(apply(abs(zabs)>1.96,2,mean))),
										 tolerance=2*sqrt(.05*.95/nsim),scale=1))
})
test_that("binary dif predictions agree within 2se",{
	expect_that(0.05,
							equals(as.numeric(max(apply(abs(zdif)>1.96,2,mean))),
										 tolerance=2*sqrt(.05*.95/nsim),scale=1))
})

#######################################################
###### check CI and prediction intervals from plot.fitEmax
######

set.seed(12357)

doselev<-c(0,5,25,50,100,350)
nd<-length(doselev)
n<-10*c(78,81,81,81,77,80)

### population parameters for simulation
e0<-2.465375 
ed50<-67.481113 
emax<-15.127726
sdy<-8.0
pop.parm<-c(log(ed50),emax,e0)    
dose<-rep(doselev,n)
meanlev<-emaxfun(doselev,pop.parm)  
meanrep<-emaxfun(dose,pop.parm)  
clev<-0.9
nsim<-1000
covci<-matrix(logical(nsim*nd),ncol=nd)
covpi<-matrix(logical(nsim*nd),ncol=nd)
covdifci<-matrix(logical(nsim*nd),ncol=nd)
covdifpi<-matrix(logical(nsim*nd),ncol=nd)
for (i in 1:nsim){
	y<-rnorm(sum(n),meanrep,sdy)
	testout<-fitEmax(y,dose,modType=4,diagnostics = FALSE)
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

test_that("plot.fitEmax CI for continuous data agree within 2.5se",{
	expect_that(clev,
							equals(as.numeric(max(apply(covci,2,mean,na.rm=TRUE))),
										 tolerance=2.5*sqrt(.1*.9/nsim),scale=1))
})
test_that("plot.fitEmax PI for continuous data agree within 2.5se",{
	expect_that(clev,
							equals(as.numeric(max(apply(covpi,2,mean,na.rm=TRUE))),
										 tolerance=2.5*sqrt(.1*.9/nsim),scale=1))
})


test_that("plot.fitEmax CI DIF for continuous data agree within 2.5se",{
	expect_that(clev,
							equals(as.numeric(max(apply(covdifci[,-1],2,mean,na.rm=TRUE))),
										 tolerance=2.5*sqrt(.1*.9/nsim),scale=1))
})
test_that("plot.fitEmax PI DIF for continuous data agree within 2.5se",{
	expect_that(clev,
							equals(as.numeric(max(apply(covdifpi[,-1],2,mean,na.rm=TRUE))),
										 tolerance=2.5*sqrt(.1*.9/nsim),scale=1))
})


#####################################
### repeat with binary data
set.seed(12357)

doselev<-c(0,5,25,50,100,350)
nd<-length(doselev)
n<-5*c(78,81,81,81,77,80)

### population parameters for simulation
e0<-qlogis(.2) 
ed50<-67.481113 
emax<-qlogis(.95)
pop.parm<-c(log(ed50),emax,e0)    
dose<-rep(doselev,n)
meanlev<-plogis(emaxfun(doselev,pop.parm))
meanrep<-plogis(emaxfun(dose,pop.parm))  
clev<-0.9
nsim<-1000
covci<-matrix(logical(nsim*nd),ncol=nd)
covpi<-matrix(logical(nsim*nd),ncol=nd)
covdifci<-matrix(logical(nsim*nd),ncol=nd)
covdifpi<-matrix(logical(nsim*nd),ncol=nd)
for (i in 1:nsim){
	y<-rbinom(sum(n),1,meanrep)
	testout<-fitEmax(y,dose,modType=4,diagnostics = FALSE,binary=TRUE)
	if(is.null(testout)){
		covci[i,]<-NA
		covpi[i,]<-NA
		covdifci[i,]<-NA
		covdifpi[i,]<-NA   
	}else{
		intout<-plot(testout,clev=clev,plot=FALSE)$plotdata
		covci[i,]<-meanlev>=intout[,'cil'] & meanlev<=intout[,'cih']
		### coverage for independent sample means from the same design
		y<-rbinom(sum(n),1,meanrep)   
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

test_that("plot.fitEmax CI for binary data agree within 2.5se",{
	expect_that(clev,
							equals(as.numeric(max(apply(covci,2,mean,na.rm=TRUE))),
										 tolerance=2.5*sqrt(.1*.9/nsim),scale=1))
})
test_that("plot.fitEmax PI for binary data agree within 2.5se",{
	expect_that(clev,
							equals(as.numeric(max(apply(covpi,2,mean,na.rm=TRUE))),
										 tolerance=2.5*sqrt(.1*.9/nsim),scale=1))
})


test_that("plot.fitEmax CI DIF for binary data agree within 2.5se",{
	expect_that(clev,
							equals(as.numeric(max(apply(covdifci[,-1],2,mean,na.rm=TRUE))),
										 tolerance=2.5*sqrt(.1*.9/nsim),scale=1))
})
test_that("plot.fitEmax PI DIF for binary data agree within 2.5se",{
	expect_that(clev,
							equals(as.numeric(max(apply(covdifpi[,-1],2,mean,na.rm=TRUE))),
										 tolerance=2.5*sqrt(.1*.9/nsim),scale=1))
})


