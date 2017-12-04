######################################################
###### DRDensity plots
library(clinDR)

if(file.exists("./clinDR/inst/tests/extraGraphics/pdfoutput")) setwd("./clinDR/inst/tests/extraGraphics/pdfoutput")

if(file.exists("output.densityplot_new.pdf")) file.rename("output.densityplot_new.pdf", "output.densityplot_old.pdf")

pdf(file=paste("output.densityplot_new.pdf"))



set.seed(12357)
data(examples14) 
dat<-examples14[[6]]
attach(dat)

fit<-fitEmax(y,dose,modType=3,count=nsize,msSat=sd^2)


dgrid<-1:100
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

qlev<-c(0.025,0.10,0.25)
fitout<-predict(fit, dosevec=dgrid,clev=0.95)
qL<-fitout$lbdif
qH<-fitout$ubdif
fitout<-predict(fit, dosevec=dgrid,clev=0.80)
qL<-cbind(qL,fitout$lbdif)
qH<-cbind(qH,fitout$ubdif)
fitout<-predict(fit, dosevec=dgrid,clev=0.50)
qL<-cbind(qL,fitout$lbdif)
qH<-cbind(qH,fitout$ubdif)

rownames(qL)<-NULL
colnames(qL)<-NULL
rownames(qH)<-NULL
colnames(qH)<-NULL


DRDensityPlot(dgrid,qL,qH,qlevL=qlev,xlab='Dose',ylab='Diff with PBO')



#####################################################
### plotBDensity

dgrid<-seq(1,100,0.5)

prior<-prior.control(epmu=1.78,epsd=1,emaxmu=-2.82,emaxsd=1,
										 p50=0.04,sigmalow=0.01,1)

fitb<-fitEmaxB(y,dose,prior,modType=4,count=nsize,msSat=sd^2)

parms<-as.matrix(fitb$estanfit)
plotBdensity(dgrid,cbind(parms[,1],parms[,2],
       parms[,3],parms[,4]))

### difference with pbo
plotBdensity(dgrid,cbind(parms[,1],parms[,2],
       parms[,3],parms[,4]),plotDif=TRUE,
       xlab='Dose',ylab='Dif with PBO')

dev.off()