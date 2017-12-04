library(clinDR)

if(file.exists("./clinDR/inst/tests/extraGraphics/pdfoutput")) setwd("./clinDR/inst/tests/extraGraphics/pdfoutput")

if(file.exists("output.plotB_new.pdf")) file.rename("output.plotB_new.pdf", "output.plotB_old.pdf")

pdf(file=paste("output.plotB_new.pdf"))


set.seed(12357)
data(examples14)
dat<-examples14[[6]]
attach(dat)

fit<-fitEmax(y,dose,modType=3,count=nsize,msSat=sd^2)

prior<-prior.control(epmu=1.78,epsd=1,emaxmu=-2.82,emaxsd=1,
										 p50=0.04,sigmalow=0.01,1)

fitb<-fitEmaxB(y,dose,prior,modType=4,count=nsize,msSat=sd^2)

parms<-as.matrix(fitb$estanfit)


##### basic plot
outB<-plotB(y,dose,parms[,1:4],(parms[,5])^2,
      ylab="Y")

plot(outB)
plot(outB, xat=c(0, 25, 100)/100)

plot(outB, log=TRUE)
plot(outB, log=TRUE, xat=c(0, 5, 25, 100)/100)


##### test as responder variable (fit is poor due to changing variances)
outBB<-plotB(y,dose,parms[,1:4],(parms[,5])^2,
      ylab="Change in LDL",binary='BinRes',BinResLev=-0.5,BinResDir='<')
plot(outBB)
plot(outBB, log=TRUE)

##### resid plot
plot(outB,plotResid=TRUE,predict=FALSE)
plot(outB,plotResid=TRUE,predict=FALSE, log=TRUE)

##### test symbol plotting
symbol<-rep(2,length(dose))
symbol[1:floor(length(dose)/2)]<-1
outB<-plotB(y,dose,parms[,1:4],(parms[,5])^2,symbol=symbol,
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

##### test active comparator
set.seed(12357)
nac<-50
msd<-median(parms[,5])
yac<-rnorm(nac,y,msd)
ac<-rnorm(nrow(parms),y,msd/sqrt(nac))

outac<-plotB(y,dose,parms[,1:4],(parms[,5])^2,
      ylab="Change in LDL vs Active Comparator",activeControl=TRUE,yac=yac,ac=ac,
	  dac=10,labac='Active',plotDif=TRUE)

plot(outac,labac='Active')
plot(outac,labac='Active', log=TRUE)

symbol<-rep(2,length(dose))
symbol[1:floor(length(dose)/2)]<-1
outac<-plotB(y,dose,parms[,1:4],(parms[,5])^2,symbol=symbol,
		symbolShape=9:10,symbolColor=c('red','blue'),symbolLabel='TESTNAME',
      ylab="LDL with Active Comparator",activeControl=TRUE,yac=yac,ac=ac,
	  dac=10,labac='Prednisone',shapeac=9,colac='red')


hold<-plotB(y,dose,parms[,1:4],(parms[,5])^2,symbol=symbol,
             symbolShape=9:10,symbolColor=c('red','blue'),symbolLabel='TESTNAME',
             ylab="LDL with Active Comparator",activeControl=TRUE,yac=yac,ac=ac,
             dac=10,labac='Prednisone',shapeac=9,colac='red')
plot(hold)

#### with responder outcome
outac<-plotB(y,dose,parms[,1:4],(parms[,5])^2,
      ylab="Change in LDL vs Active Comparator",
	  activeControl=TRUE,yac=yac,ac=ac,dac=2,
	  binary='BinRes',BinResLev=0,BinResDir = '<')

outac<-plotB(y,dose,parms[,1:4],(parms[,5])^2,symbol=symbol,
		symbolShape=9:10,symbolColor=c('red','blue'),symbolLabel='TESTNAME',
      ylab="Change in LDL vs Active Comparator",
	  activeControl=TRUE,yac=yac,ac=ac,dac=2,
	  binary='BinRes',BinResLev=0,shapeac=9,colac='red')

outac<-plotB(y,dose,parms[,1:4],(parms[,5])^2,
             ylab="Change in EDD vs Active Comparator",
             activeControl=TRUE,yac=yac,ac=ac,dac=2,
             binary='BinRes',BinResLev=0)

plot(outac,symbolShape=9:10,symbolColor=c('red','blue'),
     xat=c(0, 50, 100)/100)


dev.off()

