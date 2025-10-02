###Evaluations for print.emaxPrior
### not included in run files, separate execution required

rm(list=objects())
library(clinDR)
set.seed(12357)

priorC<-emaxPrior.control(epmu=0,epsca=4,difTargetmu=0,
												 difTargetsca=4,dTarget=20,
												 p50=(2+5)/2,
												 sigmalow=0.01,sigmaup=3)


print.emaxPrior(priorC,diffuse = TRUE,modType='4',
								docType='sap',file='vtmp1.txt')

print.emaxPrior(priorC,diffuse = FALSE,modType='4',
								docType='sap',file='vtmp2.txt')




priorB<-emaxPrior.control(epmu=0,epsca=4,difTargetmu=0,
												 difTargetsca=4,dTarget=20,
												 p50=(2+5)/2,binary=TRUE)

print.emaxPrior(priorB,diffuse = TRUE,modType='3',
								docType='sap',file='vtmp3.txt')

print.emaxPrior(priorB,diffuse = FALSE,modType='3',
								docType='sap',file='vtmp4.txt')
#################################

print.emaxPrior(priorC,diffuse = TRUE,modType='4',
								docType='prot',file='vtmp5.txt')

print.emaxPrior(priorC,diffuse = FALSE,modType='4',
								docType='protocol',file='vtmp6.txt')

print.emaxPrior(priorB,diffuse = TRUE,modType='3',
								docType='protocol',file='vtmp7.txt')

print.emaxPrior(priorB,diffuse = FALSE,modType='3',
								docType='prot',file='vtmp8.txt')

#################################
### no print with non-default meta parms
priorD<-emaxPrior.control(epmu=0,epsca=4,difTargetmu=0,
												 difTargetsca=4,dTarget=20,
												 p50=(2+5)/2,loged50mu=1.0,binary=TRUE)
print.emaxPrior(priorD,diffuse = TRUE,modType='3',
								docType='sap')

