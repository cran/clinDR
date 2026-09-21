library(clinDR)
library(testthat)
library(DoseFinding)
library(nlme)
library(parallel)
library(doParallel)
library(RBesT)
library(MASS)
library(mvtnorm)

if(file.exists("./clinDR/inst/tests")) setwd("./clinDR/inst/tests")

nprocdef<-15

RNGkind("default")
test_file('test.emaxsimBlocalParmMix.R')

RNGkind("default")
test_file('test.fitEmaxBlocalParmMix.R')

#
RNGkind("default")
test_file('test.fitEmaxB.R')

#
RNGkind("default")
test_file('test.compareTvsMix.R')

#
RNGkind("default")
test_file('test.print.emaxPrior.R')

RNGkind("default")
test_file('test.genfunctions.R')

RNGkind("default")
test_file('test.checkMonoEmaxlocalParm.R')

RNGkind("default")
test_file('test.checkMonoEmax.R')

RNGkind("default")

