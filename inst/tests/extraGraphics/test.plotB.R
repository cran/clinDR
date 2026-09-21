library(clinDR)


### note change in random number generator and set.seed
### produced some changes with earlier versions of the
### plots


if(file.exists("./clinDR/inst/tests/extraGraphics/pdfoutput")){
	pvar<-"./clinDR/inst/tests/extraGraphics/pdfoutput"
} else pvar<-NULL

if(file.exists(file.path(pvar,"output.plotB_new.pdf"))) file.rename(file.path(pvar,"output.plotB_new.pdf"), 
																																				 file.path(pvar,"output.plotB_old.pdf"))
	
pdf(file=paste(file.path(pvar,"output.plotB_new.pdf")))


set.seed(12357)
data(metaData) 
dat<-metaData[metaData$taid==6 & metaData$poptype==1,]
attach(dat)

msSat<-sum((dat$sampsize-1)*(dat$sd)^2)/(sum(dat$sampsize)-length(dat$sampsize))
fit<-fitEmax(rslt,dose,modType=3,count=sampsize,msSat=msSat)

prior<-emaxPrior.control(epmu=0,epsca=10,difTargetmu=0,difTargetsca=10,dTarget=80.0,
												 p50=3.75,sigmalow=0.01,sigmaup=20)

fitb<-fitEmaxB(rslt,dose,prior,modType=4,count=sampsize,msSat=msSat,
							 mcmc=mcmc.control(iter=20000))
parms<-coef(fitb)
sigsim<-sigma(fitb)


##### basic plot
outB<-plotB(rslt,dose,parms[,1:4],(sigsim)^2,count=sampsize,
      ylab="Y",xat=c(0, 40, 80))

plot(outB)
plot(outB, xat=c(0, 40, 80))

plot(outB, log=TRUE)
plot(outB, log=TRUE, xat=c(0, 5, 40, 80))


##### resid plot
plot(outB,plotResid=TRUE,predict=FALSE)
plot(outB,plotResid=TRUE,predict=FALSE, log=TRUE)

##### test symbol plotting
symbol<-rep(2,length(dose))
symbol[1:floor(length(dose)/2)]<-1
outB<-plotB(rslt,dose,parms[,1:4],(sigsim)^2,count=sampsize,symbol=symbol,
      ylab="Cholesterol")

plot(outB,symbolShape=8:9,symbolColor=c('red','blue'),
	 symbolLabel='TESTNAME')
plot(outB,symbolShape=8:9,symbolColor=c('red','blue'),
     symbolLabel='TESTNAME', log=TRUE)


#### use medians and no prediction
plot(outB,symbolShape=8:9,symbolColor=c('red','blue'),
	 symbolLabel='TESTNAME',plotMed=TRUE,predict=FALSE)

plot(outB,symbolShape=8:9,symbolColor=c('red','blue'),
     symbolLabel='TESTNAME',plotMed=TRUE,predict=FALSE, log=TRUE)

### residual plot
plot(outB,symbolShape=8:9,symbolColor=c('red','blue'),
	 symbolLabel='TESTNAME',plotResid=TRUE)

plot(outB,symbolShape=8:9,symbolColor=c('red','blue'),
     symbolLabel='TESTNAME',plotResid=TRUE, log=TRUE)

##### test as responder variable (Note: with small sample sizes
##### create patient level data   the predictive intervals and
#####                             points are discrete with big
#####                             jumps
ymeans<-predict(fitb,dose)$pred
ymeans<-rep(ymeans,sampsize)
yvec<-rnorm(length(ymeans),ymeans,sd)
dvec<-rep(dose,sampsize)
outBB<-plotB(yvec,dvec,parms[,1:4],(sigsim)^2,
      ylab="Change in LDL", log=TRUE)
plot(outBB,log=TRUE)
outBB<-plotB(yvec,dvec,parms[,1:4],(sigsim)^2,
      ylab="Change in LDL",binary='BinRes',BinResLev=-40,BinResDir='<')
plot(outBB,log=TRUE)
plot(outBB)

##### test active comparator
set.seed(12357)
nac<-10
msd<-median(sigsim)
yac<-rnorm(nac,mean(rslt),msd)
ac<-rnorm(nrow(parms),mean(yac),msd/sqrt(nac))

outac<-plotB(rslt,dose,parms[,1:4],(sigsim)^2,count=sampsize,
      ylab="Change in LDL vs Active Comparator",activeControl=TRUE,yac=yac,ac=ac,
	  labac='Active',plotDif=TRUE)

plot(outac,labac='Active',xat=c(0,.25,.5,.75,1))
plot(outac,labac='Active', log=TRUE, xat=c(0,.1,.5,1))

symbol<-rep(2,length(dose))
symbol[1:floor(length(dose)/2)]<-1
outac<-plotB(rslt,dose,parms[,1:4],(sigsim)^2,count=sampsize,symbol=symbol,
		symbolShape=9:10,symbolColor=c('red','blue'),symbolLabel='TRT GROUP',
      ylab="LDL with Active Comparator",activeControl=TRUE,yac=yac,ac=ac,
	  labac='Prednisone',shapeac=8,colac='green',xat=c(0,.25,.5,.75,1))


plot(outac,xat=c(0,.25,.5,.75,1),
             symbolShape=9:10,symbolColor=c('red','blue'),symbolLabel='TESTNAME',
             ylab="LDL with Active Comparator",activeControl=TRUE,yac=yac,ac=ac,
             labac='Prednisone',shapeac=8,colac='green')

#### with responder outcome
symvec<-rep(symbol,sampsize)
outac<-plotB(yvec,dvec,parms[,1:4],(sigsim)^2,
      ylab="Change in LDL vs Active Comparator",
	  activeControl=TRUE,yac=yac,ac=ac,
	  binary='BinRes',BinResLev=-40,BinResDir = '<')

outac<-plotB(yvec,dvec,parms[,1:4],(sigsim)^2,symbol=symvec,
		symbolShape=9:10,symbolColor=c('red','blue'),symbolLabel='TESTNAME',
      ylab="Change in LDL vs Active Comparator",
	  activeControl=TRUE,yac=yac,ac=ac,
	  binary='BinRes',BinResLev=-40,shapeac=8,colac='orange')

outac<-plotB(yvec,dvec,parms[,1:4],(sigsim)^2,
             ylab="Change in EDD vs Active Comparator",
             activeControl=TRUE,yac=yac,ac=ac,
             binary='BinRes',BinResLev=-40)

plot(outac,symbolShape=9:10,symbolColor=c('red','blue'),
     xat=c(0, 40, 80))


dev.off()
detach(dat)
