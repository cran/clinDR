library(clinDR)
if(file.exists("./clinDR/inst/tests/extraGraphics/pdfoutput")) setwd("./clinDR/inst/tests/extraGraphics/pdfoutput")

if(file.exists("output.plotD_new.pdf")) file.rename("output.plotD_new.pdf", "output.plotD_old.pdf")

pdf(file=paste("output.plotD_new.pdf"))



data(examples14) 
dat<-examples14[[6]]



with(dat,plotD(y,dose,sem=sem,meansOnly=TRUE,ylab=
                 "Change from baseline LDL",xlab="Dose (mg)"))


with(dat,plotD(y,dose,sem=sem,meansOnly=TRUE,ylab=
                 "Change from baseline LDL",xlab="Log Dose (mg)",log=TRUE))

dev.off()
