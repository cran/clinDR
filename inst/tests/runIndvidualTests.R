library(clinDR)
library(testthat)
library(DoseFinding)


if(file.exists("./clinDR/inst/tests")) setwd("./clinDR/inst/tests")

test_file('test.emaxsim.R')

test_file('test.emaxsimB.R')

test_file('test.SeEmax.R')

test_file('test.predict.emaxsim.R')

test_file('test.update.emaxsimobj.R')

test_file('test.targetcodes.R')

test_file('test.emaxfun.R')

test_file('test.fitEmax.R')

test_file('test.fitEmaxB.R')

test_file('test.genfunctions.R')

test_file('test.checkMonoEmax.R')

