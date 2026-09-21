######################################################
###### DRDensity plots
library(clinDR)

if(file.exists("./clinDR/inst/tests/extraGraphics/pdfoutput")){
	pvar<-"./clinDR/inst/tests/extraGraphics/pdfoutput"
} else pvar<-NULL

if(file.exists(file.path(pvar,"output.densityplot_new.pdf"))) file.rename(file.path(pvar,"output.densityplot_new.pdf"), 
																																				 file.path(pvar,"output.densityplot_old.pdf"))
	
pdf(file=paste(file.path(pvar,"output.densityplot_new.pdf")))




set.seed(12357)
data(metaData) 
dat<-metaData[metaData$taid==6 & metaData$poptype==1,]
attach(dat)

msSat<-sum((dat$sampsize-1)*(dat$sd)^2)/(sum(dat$sampsize)-length(dat$sampsize))

fit<-fitEmax(rslt,dose,modType=3,count=sampsize,msSat=msSat)


dgrid<-0:100
fitout<-predict(fit, dosevec=dgrid,clev=0.95)
qL<-fitout$lbdif
qH<-fitout$ubdif
fitout<-predict(fit, dosevec=dgrid,clev=0.90)
qL<-cbind(qL,fitout$lbdif)
qH<-cbind(qH,fitout$ubdif)
fitout<-predict(fit, dosevec=dgrid,clev=0.80)
qL<-cbind(qL,fitout$lbdif)
qH<-cbind(qH,fitout$ubdif)
fitout<-predict(fit, dosevec=dgrid,clev=0.50)
qL<-cbind(qL,fitout$lbdif)
qH<-cbind(qH,fitout$ubdif)

qlev<-c(0.025,0.05,0.10,0.25)

rownames(qL)<-NULL
colnames(qL)<-NULL
rownames(qH)<-NULL
colnames(qH)<-NULL

DRDensityPlot(dgrid,qL,qH,qlevL=qlev,xlab='Dose',ylab='Diff with PBO')


#####################################################
### plotBDensity

dgrid<-seq(0,100,0.5)

prior<-emaxPrior.control(epmu=0,epsca=10,difTargetmu=0,difTargetsca=10,dTarget=80.0,
												 p50=3.75,sigmalow=0.01,sigmaup=20)
mcmc<-mcmc.control(chains=3)


fitb<-fitEmaxB(rslt,dose,prior,modType=4,count=sampsize,msSat=msSat)

parms<-coef(fitb)
pout<-plotBdensity(dgrid,parms[,1:4])

### difference with pbo
poutdif<-plotBdensity(dgrid,parms[,1:4],plotDif=TRUE,
       xlab='Dose',ylab='Dif with PBO')

dev.off()
detach(dat)