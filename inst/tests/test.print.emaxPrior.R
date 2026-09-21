context('print.emaxPrior evaluations')
################################################################################
######             Test Code for print.emaxPrior                          ######
################################################################################

### old results for comparison
old<-new.env()
load('currentPriors.RData',env=old)

## binary

# pbo diffuse, eff informative (unusual)
prior.m.diffuse <- emaxPrior.control(mixP=3, mu_ep=c(-1.74, -1,87, -1.73),
                                     sd_ep= c(2, 2, 2), w_ep=c(0.34, 0.16, 0.5),
                                     difTargetmu = qlogis(0.7) - qlogis(0.15),
                                     difTargetsca= 1,
                                     dTarget = 100,
                                     p50 = 30,
                                     binary=TRUE)
#Diffuse t prior: Both epsca and difTargetsca > 2
prior.t.diffuse <- emaxPrior.control(epmu=0, epsca=4,
                                     difTargetmu = qlogis(0.7) - qlogis(0.15),
                                     difTargetsca= 3,
                                     dTarget = 100,
                                     p50 = 30,
                                     binary=TRUE)

# Informative mixture prior
prior.m.informative <- emaxPrior.control(mixP=3, mu_ep=c(-1.74, -1,87, -1.73),
                                         sd_ep= c(0.77, 1.62, 7.84), w_ep=c(0.34, 0.16, 0.5),
                                         difTargetmu = qlogis(0.7) - qlogis(0.15),
                                         difTargetsca= 1,
                                         dTarget = 100,
                                         p50 = 30,
                                         binary=TRUE)
# Informative t prior: At least one of the scale parameters epsca / difTargetsca < 2
prior.t.informative<- emaxPrior.control(epmu=qlogis(0.15), epsca=1,
                                        difTargetmu = qlogis(0.7) - qlogis(0.15),
                                        difTargetsca= 2,
                                        dTarget = 100,
                                        p50 = 30,
                                        binary=TRUE)

suppressWarnings(bmds<-print.emaxPrior(x=prior.m.diffuse, file=NULL, doc=TRUE, 
                      diffuse=c(TRUE,FALSE), docType = "sap"))
test_that('bmds agree',{expect_identical(bmds,old$bmds)})
# check creation of word docx
#bmds<-print.emaxPrior(x=prior.m.diffuse, doc=TRUE, diffuse=c(TRUE,FALSE),
#                      docType = "sap")

btds<-print.emaxPrior(x=prior.t.diffuse, file=NULL, doc=TRUE, 
                      diffuse=c(TRUE,TRUE), docType = "sap")
test_that('btds agree',{expect_identical(btds,old$btds)})


bmis<-print.emaxPrior(x=prior.m.informative, file=NULL, doc=TRUE, 
                      diffuse=c(pbo=FALSE,eff=FALSE), docType = "sap")
test_that('bmis agree',{expect_identical(bmis,old$bmis)})


btis<-print.emaxPrior(x=prior.t.informative, file=NULL, doc=TRUE, 
                      diffuse=c(pbo=FALSE,eff=TRUE), docType = "sap")
test_that('btis agree',{expect_identical(btis,old$btis)})


suppressWarnings(bmdp<-print.emaxPrior(x=prior.m.diffuse, file=NULL, doc=TRUE, 
                      diffuse=c(TRUE,FALSE), docType = "protocol"))
test_that('bmdp agree',{expect_identical(bmdp,old$bmdp)})

btdp<-print.emaxPrior(x=prior.t.diffuse, file=NULL, doc=TRUE, 
                      diffuse=c(TRUE,TRUE), docType = "protocol")
test_that('btdp agree',{expect_identical(btdp,old$btdp)})

bmip<-print.emaxPrior(x=prior.m.informative, file=NULL, doc=TRUE, 
                      diffuse=c(pbo=FALSE,eff=FALSE), docType = "protocol")
test_that('bmip agree',{expect_identical(bmip,old$bmip)})

btip<-print.emaxPrior(x=prior.t.informative, file=NULL, doc=TRUE, 
                      diffuse=c(pbo=FALSE,eff=TRUE), docType = "protocol")
test_that('btip agree',{expect_identical(btip,old$btip)})



# continuous

# Mixture prior so no diffuse check.  Can call with any diffuse setting for test
prior.m.cont <- emaxPrior.control(mixP=3, mu_ep=c(-48.6, -52.2, -50),
                                  sd_ep= c(11.8, 38.4, 88.0), w_ep=c(0.37, 0.13, 0.5),
                                  difTargetmu = -150,
                                  difTargetsca= 30,
                                  dTarget = 100,
                                  p50 = 30,
                                  sigmalow=10,sigmaup=30)

prior.t.cont <- emaxPrior.control(epmu=-50, epsca=100,
                                  difTargetmu = -150,
                                  difTargetsca= 30,
                                  dTarget = 100,
                                  p50 = 30,
                                  sigmalow=10,
                                  sigmaup=30
)


suppressWarnings(cmds<-print.emaxPrior(x=prior.m.cont, doc=TRUE, diffuse=c(TRUE,FALSE), 
                      file = NULL, docType = "sap"))
test_that('cmds agree',{expect_identical(cmds,old$cmds)})

ctds<-print.emaxPrior(x=prior.t.cont, file=NULL, doc=TRUE, 
                      diffuse=c(TRUE,TRUE), docType = "sap")
test_that('ctds agree',{expect_identical(ctds,old$ctds)})

cmis<-print.emaxPrior(x=prior.m.cont, doc=TRUE, diffuse=c(FALSE,TRUE), 
                      file = NULL, docType = "sap")
test_that('cmis agree',{expect_identical(cmis,old$cmis)})

ctis<-print.emaxPrior(x=prior.t.cont, doc=TRUE, diffuse=c(FALSE,FALSE), 
                      file = NULL, docType = "sap")
test_that('ctis agree',{expect_identical(ctis,old$ctis)})

suppressWarnings(cmdp<-print.emaxPrior(x=prior.m.cont, doc=TRUE, diffuse=c(TRUE,FALSE), 
                      file = NULL, docType = "protocol"))
test_that('cmdp agree',{expect_identical(cmdp,old$cmdp)})


ctdp<-print.emaxPrior(x=prior.t.cont, file=NULL, doc=TRUE, 
                      diffuse=c(TRUE,TRUE), docType = "protocol")
test_that('ctdp agree',{expect_identical(ctdp,old$ctdp)})

cmip<-print.emaxPrior(x=prior.m.cont, doc=TRUE, diffuse=c(FALSE,TRUE), 
                      file = NULL, docType = "protocol")
test_that('cmip agree',{expect_identical(cmip,old$cmip)})

ctip<-print.emaxPrior(x=prior.t.cont, doc=TRUE, diffuse=c(FALSE,FALSE), 
                      file = NULL, docType = "protocol")
test_that('ctip agree',{expect_identical(ctip,old$ctip)})



## code executed to create regression tests for future testing
#save(bmds,btds,bmis,btis,bmdp,btdp,bmip,btip,cmds,ctds,cmis,ctis,
#     cmdp,ctdp,cmip,ctip,
#     file='currentPriors.RData')
