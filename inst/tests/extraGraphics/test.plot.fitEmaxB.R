library(clinDR)


if(file.exists("./clinDR/inst/tests/extraGraphics/pdfoutput")) setwd("./clinDR/inst/tests/extraGraphics/pdfoutputB")

if(file.exists("output.fitEmaxB_new.pdf")) file.rename("output.fitEmaxB_new.pdf", "output.fitEmaxB_old.pdf")

pdf(file="output.fitEmaxB_new.pdf")


set.seed(12357)

doselev<-c(0,5,25,50,100,350)
n<-c(78,81,81,81,77,80)

### population parameters for simulation
e0<-2.465375 
ed50<-67.481113 
emax<-15.127726
sdy<-8.0
pop.parm<-c(log(ed50),emax,e0)    
dose<-rep(doselev,n)
meanlev<-emaxfun(dose,pop.parm)  

y<-rnorm(sum(n),meanlev,sdy)


prior<-prior.control(0,30,0,30,50,0.1,30,edDF=5)
mcmc<-mcmc.control(chains=3,warmup=500,iter=3000,seed=53453,propInit=0.15,adapt_delta = .9)

testout<-fitEmaxB(y,dose,prior=prior,modType=4,
								 mcmc=mcmc)

### basic test
plot(testout)
#plot(testout, xat=c(-100, 0, 5, 100, 350)) ## fail planned
#plot(testout, xat=c(0, 5, 100, 350, 400))  ## fail planned
plot(testout, xat=c(0, 100, 350))

plot(testout, log=TRUE)
#plot(testout, log=TRUE, xat=c(-100, 0, 5, 100, 350)) ## fail planned
#plot(testout, log=TRUE, xat=c(0, 5, 100, 350, 400))  ## fail planned
plot(testout, log=TRUE, xat=c(0, 5, 100, 350))

plot(testout,plotDif=TRUE)
plot(testout,plotDif=TRUE, log=TRUE)

### with symbol
symbol<-1+1*(dose==100)
plot(testout,symbol=symbol)
plot(testout,symbol=symbol, log=TRUE)

plot(testout,symbol=symbol,
	 symbolLabel='TESTGROUP',
	 symbolShape=c(8,10),symbolColor=c('blue','red'))
plot(testout,symbol=symbol, log=TRUE,
     symbolLabel='TESTGROUP',
     symbolShape=c(8,10),symbolColor=c('blue','red'))

plot(testout,symbol=symbol,
	 symbolLabel='TESTGROUP',
	 symbolShape=c(8,10),symbolColor=c('blue','red'),plotDif=TRUE)
plot(testout,symbol=symbol, log=TRUE,
     symbolLabel='TESTGROUP',
     symbolShape=c(8,10),symbolColor=c('blue','red'),plotDif=TRUE)

### residual plot

plot(testout,plotResid=TRUE)
plot(testout,plotResid=TRUE, log=TRUE)
#####################################################################
##### tests with multiple protocols and 3-parm emax
##### one protocol much larger than other
set.seed(12357)
doselev<-c(0,5,25,50,100,350)
n<-80

### population parameters for simulation
e0<-2.465375 
ed50<-67.481113 
emax<-15.127726
sdy<-8.0
pop.parm<-c(log(ed50),emax,e0)    
dose<-rep(doselev,n)
meanlev<-emaxfun(dose,pop.parm)  

y<-rnorm(n*length(doselev),meanlev,sdy)

prot<-rep('a',length(y))
prot[1:(5*length(doselev))]<-'b'

symbol=1*(dose<300)+2*(dose>=300)
symbolLabel='TESTNAME'
symbolColor=c('red','blue','green')
symbolShape=c(8,10)

prior<-prior.control(0,30,0,30,50,0.1,30,edDF=5)
mcmc<-mcmc.control(chains=3,warmup=500,iter=3000,seed=53453,propInit=0.15,adapt_delta = .9)

testout2<-fitEmaxB(y,dose,prior=prior,modType=3,prot=prot,mcmc=mcmc)

plot(testout2,bwidth=10,symbol=symbol,symbolLabel=symbolLabel,
	 symbolShape=symbolShape,symbolColor=symbolColor)

plot(testout2,bwidth=10,symbol=symbol,symbolLabel=symbolLabel,
     symbolShape=symbolShape,symbolColor=symbolColor, log=TRUE)

plot(testout2,bwidth=10,symbol=symbol,symbolLabel=symbolLabel,
	 symbolShape=symbolShape,symbolColor=symbolColor,plotDif=TRUE)

plot(testout2,bwidth=10,symbol=symbol,symbolLabel=symbolLabel,
     symbolShape=symbolShape,symbolColor=symbolColor,plotDif=TRUE, log=TRUE)

testout3<-fitEmaxB(y,dose,prior=prior,modType=3,mcmc=mcmc)

symbmod<-symbol
symbmod[prot==1 & dose==350]<-1
plot(testout3,bwidth=10,symbol=symbmod,symbolLabel=symbolLabel,
	 symbolShape=symbolShape,symbolColor=symbolColor)
plot(testout3,bwidth=10,symbol=symbmod,symbolLabel=symbolLabel,
     symbolShape=symbolShape,symbolColor=symbolColor, log=TRUE)


###############################################################
#### 3 parm, int=2, starting value given
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
pop.parm<-c(log(ed50),emax,e0)    
dose<-rep(doselev,n)
meanlev<-emaxfun(dose,pop.parm)  


y<-rnorm(n1+n2,meanlev,sdy)
prots<-c(rep('d',n1),rep('ee',n2))

ysub<-y[dose!=0]
dsub<-dose[dose!=0]
protsub<-prots[dose!=0]

prior<-prior.control(0,30,0,30,50,0.1,30,edDF=5)
mcmc<-mcmc.control(chains=3,warmup=500,iter=3000,seed=53453,propInit=0.15,adapt_delta = .9)

testout<-fitEmaxB(ysub,dsub,prior=prior,modType=3,prot=protsub,pboAdj=TRUE,mcmc=mcmc)

plot(testout)
plot(testout, log=TRUE)

testout4<-fitEmaxB(ysub,dsub,prior=prior,
				  modType=4,prot=protsub,pboAdj=TRUE,mcmc=mcmc)

plot(testout4)
plot(testout4, log=TRUE)

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
prots<-c(rep('eg',2*nd1),rep('fg',2*nd2))
dvec<-c(dvec1,dvec1,dvec2,dvec2)

prior<-prior.control(0,4,0,4,.5,edDF=5,binary=TRUE)
mcmc<-mcmc.control(chains=3,warmup=500,iter=3000,seed=53453,propInit=0.15,adapt_delta = .9)


testout<-fitEmaxB(y,dvec,modType=4,
									prot=prots,
									count=counts,binary=TRUE,
									prior=prior,mcmc=mcmc,	
									diagnostics=FALSE)


plot(testout)
plot(testout, log=TRUE)
plot(testout, log=TRUE, xat=c(0.1, 0.2, 0.4, 0.6, 1.0))


plot(testout,plotDif=TRUE)
plot(testout,plotDif=TRUE, log=TRUE)
plot(testout,plotDif=TRUE, log=TRUE, xat=c(0.1, 0.2, 0.4, 0.6, 1.0))


testout2<-fitEmaxB(y[prots=='fg'],dvec[prots=='fg'],prior=prior,modType=4,
          count=counts[prots=='fg'],binary=TRUE,diagnostics=FALSE,
          mcmc=mcmc)

plot(testout2)
plot(testout2, log=TRUE)

plot(testout2,plotDif=TRUE)
plot(testout2,plotDif=TRUE, log=TRUE)


dev.off()


