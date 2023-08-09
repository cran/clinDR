#######################################
######## contB
########
library(tidyverse)
library(clinDR)


# Design aspects
nsim <- 50
doselev <- c(0,5,25,50,150)
n <- c(100, 50, 50, 50, 100)
Ndose <- length(doselev)
 
# Parameters
led50 <- log(15)
lambda <- 1
e0 <- 10
target <- 6
targetDose <- 175
emax <- solveEmax(target = target, 
                        dose = targetDose,
                        led50 = led50,
                        lambda = lambda,
                        e0 = e0,
                        pboadj = TRUE)

pop <- c(led50 = led50, lambda = lambda, emax = emax, e0 = e0)
resSD <- 8

meanlev <- emaxfun(doselev, pop)
gen <- FixedMean(n, doselev, meanlev, resSD, parm = pop)
 
# Priors
prior <- emaxPrior.control(epmu = 5,
	epsca = 80,
	difTargetmu = 0,
	difTargetsca = 80,
	dTarget = 175,
	p50 = 15,
	sigmalow = 0.8,
	sigmaup = 70,
	parmDF = 5)
# Stan and MCMC settings
rstan::rstan_options(auto_write = TRUE)
mcmc<-mcmc.control(chains=1,warmup=500,iter=5000,
	propInit=0.15,adapt_delta = 0.95)
 
# Simulate outcomes
D1 <- emaxsimB(nsim, gen, prior,
	modType=4,mcmc=mcmc,
	check=FALSE, nproc=10,
	binary = FALSE, seed=12357)
summary(D1)
plot(D1)
 
# Additional plots of simulation results
pairv<-select(mutate(as_tibble(D1$est),ED50 = exp(led50)),-led50)
pairv<-rename(pairv,E0 = e0, 
	Emax = emax) 
pairv<-mutate(pairv,'Residual SD' = D1$residSD)
pairs(pairv)

#######################################
######## cont
########
library(tidyverse)
library(clinDR)


# Design aspects
nsim <- 30
doselev <- c(0,5,25,50,150)
n <- c(100, 50, 50, 50, 150)
Ndose <- length(doselev)
 
# mean values
meanlev <- c(2.0, 3.72, 5.31, 7.31, 8.0)
resSD <- 6

gen <- FixedMean(n, doselev, meanlev, resSD)
 

 
# Simulate outcomes
D1 <- emaxsim(nsim, gen, modType=4,
	nproc=5, binary = FALSE,
	seed=12357)

summary(D1, testalph=0.05)
plot(D1)
 
# Additional plots of simulation results
pairv<-select(mutate(as_tibble(D1$est4),ED50 = exp(led50)),-led50)
pairv<-rename(pairv,E0 = e0, 
	Emax = emax) 
pairv<-mutate(pairv,'Residual SD' = D1$residSD)
pairs(pairv)


#######################################
######## binB
########
library(tidyverse)
library(clinDR)


# Design aspects
nsim <- 50
doselev <- c(0,5,25,50,150,300)
n <- c(100, 50, 50, 50, 150,50)
Ndose <- length(doselev)
 
# Parameters
led50 <- log(15)
lambda <- 1
e0 <- -1
target <- 1.5
targetDose <- 200
emax <- solveEmax(target = target, 
                        dose = targetDose,
                        led50 = led50,
                        lambda = lambda,
                        e0 = e0,
                        pboadj = TRUE)

pop <- c(led50 = led50, lambda = lambda, emax = emax, e0 = e0)
meanlev <- plogis(emaxfun(doselev, pop))
gen <- FixedMean(n, doselev, meanlev, parm = pop, binary = TRUE)
 
# Priors
prior <- emaxPrior.control(epmu = 0,
	epsca = 4,
	difTargetmu = 0,
	difTargetsca = 4,
	dTarget = 200,
	p50 = 25,
	parmDF = 5,
	binary = TRUE)
# Stan and MCMC settings
rstan::rstan_options(auto_write = TRUE)
mcmc<-mcmc.control(chains=1,warmup=500,iter=5000,
	propInit=0.15,adapt_delta = 0.95)
 
# Simulate outcomes
D1 <- emaxsimB(nsim, gen, prior,
	modType=4,mcmc=mcmc,
	check=FALSE, nproc=10,
	binary = TRUE, seed=12357)
summary(D1)
plot(D1)
 
# Additional plots of simulation results
pairv<-select(mutate(as_tibble(D1$est),ED50 = exp(led50)),-led50)
pairs(pairv)


#######################################
######## bin
########
library(tidyverse)
library(clinDR)


# Design aspects
nsim <- 50
doselev <- c(0,5,25,50,150,300)
n <- c(100, 50, 50, 50, 150,50)
Ndose <- length(doselev)
 
# Parameters
led50 <- log(15)
lambda <- 1.5
e0 <- -1
target <- 1.5
targetDose <- 200
emax <- solveEmax(target = target, 
                        dose = targetDose,
                        led50 = led50,
                        lambda = lambda,
                        e0 = e0,
                        pboadj = TRUE)

pop <- c(led50 = led50, lambda = lambda, emax = emax, e0 = e0)
meanlev <- plogis(emaxfun(doselev, pop))
gen <- FixedMean(n, doselev, meanlev, parm = pop, binary = TRUE)
 

 
# Simulate outcomes
D1 <- emaxsim(nsim, gen, modType=4,
	nproc=10, binary = TRUE,
	seed=12357)

summary(D1, testalph=0.05)
plot(D1)
 
# Additional plots of simulation results
pairv<-select(mutate(as_tibble(D1$est4),ED50 = exp(led50)),-led50)
pairs(pairv)



