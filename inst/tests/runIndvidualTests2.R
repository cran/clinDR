library(clinDR)
library(testthat)
library(DoseFinding)
library(parallel)
library(doParallel)


if(file.exists("./clinDR/inst/tests")) setwd("./clinDR/inst/tests")

nprocdef<-15

# RNGkind("default")
# test_file('test.emaxsim.R')
# 
# RNGkind("default")
# test_file('test.emaxsimBlocalParm.R')
# 
# RNGkind("default")
# test_file('test.emaxsimB.R')
# 
# RNGkind("default")
# test_file('test.SeEmax.R')
# 
# RNGkind("default")
# test_file('test.predict.emaxsim.R')
# 
# RNGkind("default")
# test_file('test.update.emaxsimobj.R')
# 
# RNGkind("default")
# test_file('test.targetcodes.R')
# 
# RNGkind("default")
# test_file('test.emaxfun.R')
#
#RNGkind("default")
#test_file('test.fitEmax.R')
#
#RNGkind("default")
#test_file('test.fitEmaxBlocalParm.R')
#
RNGkind("default")
test_file('test.fitEmaxB.R')

RNGkind("default")
test_file('test.genfunctions.R')

RNGkind("default")
test_file('test.checkMonoEmaxlocalParm.R')

RNGkind("default")
test_file('test.checkMonoEmax.R')

RNGkind("default")

