context('emaxfun')

set.seed(12357)

e0<-rnorm(50,0,1)
emax<-rnorm(50,10,1)
led50<-log(rnorm(50,10,.1))
lambda<-rnorm(50,1,.01)
parm<-cbind(led50,lambda,emax,e0)
ed50<-exp(led50)

dose<-c(1:5)

### 4-parm tests

# 1 set of parameters

outnew<-emaxfun(dose,c(led50[1],lambda[1],emax[1],e0[1]))
outold<-e0[1]+(emax[1]*dose^lambda[1])/(dose^lambda[1]+ed50[1]^lambda[1])

test_that("4 parm multi-dose, 1 parm set",{
    expect_that(outnew,equals(outold))
})

outnew<-emaxfun(dose,parm)
outold1<-e0+(emax*dose[1]^lambda)/(dose[1]^lambda+ed50^lambda)
outold4<-e0+(emax*dose[4]^lambda)/(dose[4]^lambda+ed50^lambda)

test_that("4 parm multi-dose, multi-parm set",{
    expect_that(outnew[,1],equals(outold1))
    expect_that(outnew[,4],equals(outold4))
})

### 3-parm tests
# 1 set of parameters

outnew<-emaxfun(dose,c(led50[1],emax[1],e0[1]))
outold<-e0[1]+(emax[1]*dose)/(dose+ed50[1])

test_that("3 parm multi-dose, 1 parm set",{
    expect_that(outnew,equals(outold))
})

parm<-cbind(led50,emax,e0)
outnew<-emaxfun(dose,parm)
outold1<-e0+(emax*dose[1])/(dose[1]+ed50)
outold4<-e0+(emax*dose[4])/(dose[4]+ed50)

test_that("4 parm multi-dose, multi-parm set",{
    expect_that(outnew[,1],equals(outold1))
    expect_that(outnew[,4],equals(outold4))
})


e0<-0.75
dose<-1
led50<-log(0.5)
lambda<-2
target<-1.5
emax<-solveEmax(target,dose,led50,lambda,e0)

dose1<-solveDose(target,led50,lambda,emax,e0)

eout<-emaxfun(dose=dose1,parm=c(led50,lambda,emax,e0))-e0

test_that("solveDose,solveEmax positive effect",{
    expect_that(dose1,equals(1))
    expect_that(eout,equals(target))
})

e0<-10
dose<-1
led50<-log(0.5)
lambda<-2
target<- -1.5
emax<-solveEmax(target,dose,led50,lambda,e0)

dose1<-solveDose(target,led50,lambda,emax,e0)

eout<-emaxfun(dose=dose1,parm=c(led50,lambda,emax,e0))-e0

test_that("solveDose,solveEmax negative effect",{
    expect_that(dose,equals(1))
    expect_that(eout,equals(target))
})

###
e0<-0.75
dose<-1
led50<-log(0.5)
lambda<-2
target<-1.5
emax<-solveEmax(target,dose,led50,lambda,e0,pboadj=FALSE)

dose1<-solveDose(target,led50,lambda,emax,e0,pboadj=FALSE)

eout<-emaxfun(dose=dose1,parm=c(led50,lambda,emax,e0))

test_that("solveDose,solveEmax positive effect",{
    expect_that(dose1,equals(1))
    expect_that(eout,equals(target))
})

e0<-10
dose<-1
led50<-log(0.5)
lambda<-2
target<- -1.5
emax<-solveEmax(target,dose,led50,lambda,e0,pboadj=FALSE)

dose1<-solveDose(target,led50,lambda,emax,e0,pboadj=FALSE)

eout<-emaxfun(dose=dose1,parm=c(led50,lambda,emax,e0))

test_that("solveDose,solveEmax negative effect",{
    expect_that(dose,equals(1))
    expect_that(eout,equals(target))
})

