library(clinDR)


if(file.exists("./clinDR/inst/tests/extraGraphics/pdfoutput")){
	pvar<-"./clinDR/inst/tests/extraGraphics/pdfoutput"
} else pvar<-NULL

if(file.exists(file.path(pvar,"output.plotD_new.pdf"))) file.rename(file.path(pvar,"output.plotD_new.pdf"), 
																																				 file.path(pvar,"output.plotD_old.pdf"))
	
pdf(file=paste(file.path(pvar,"output.plotD_new.pdf")))



data(metaData) 
dat<-metaData[metaData$taid==6 & metaData$poptype==1,]

with(dat,plotD(rslt,dose,sem=se,meansOnly=TRUE,ylab=
                 "Change from baseline LDL",xlab="Dose (mg)"))


with(dat,plotD(rslt,dose,sem=se,meansOnly=TRUE,ylab=
                 "Change from baseline LDL",xlab="Log Dose (mg)",log=TRUE))

dev.off()
