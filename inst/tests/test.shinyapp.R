###
### first continuous example
### bayes, 4-parm, emax mod parms input
library(clinDR)
library(tidyverse)

set.seed(12357)

# Design aspects
nsim <- 48
doselev <- c(0,5,25,50,100)
n <- c(78,81,81,81,77)
Ndose <- length(doselev)

# Parameters
led50 <- log(67.481113)
lambda <- 1.25
e0 <- 2.465375
target <- 2.464592
targetDose <- 100
emax <- solveEmax(target = target, 
									dose = targetDose,
									led50 = led50,
									lambda = lambda,
									e0 = e0,
									pboadj = TRUE)

pop <- c(led50 = led50, lambda = lambda, emax = emax, e0 = e0)
resSD <- 7.967897

meanlev <- emaxfun(doselev, pop)
gen <- FixedMean(n, doselev, meanlev, resSD, parm = pop)

# Priors
prior <- emaxPrior.control(epmu = 0,
													 epsca = 30,
													 difTargetmu = 0,
													 difTargetsca = 30,
													 dTarget = 100,
													 p50 = 50,
													 sigmalow = 0.1,
													 sigmaup = 30,
													 parmDF = 5)
# Stan and MCMC settings
rstan::rstan_options(auto_write = TRUE)
mcmc<-mcmc.control(chains=1,warmup=500,iter=5000,
									 propInit=0.15,adapt_delta = 0.95)

# Simulate outcomes
D1 <- emaxsimB(nsim, gen, prior,
							 modType=4,mcmc=mcmc,
							 check=FALSE, nproc=8, binary = FALSE)
summary(D1)
plot(D1)

# Additional plots of simulation results
library(tidyverse,quietly = TRUE)

D1$est %>%
	as_tibble() %>%
	mutate(ED50 = exp(led50)) %>%
	select(-led50) %>%
	rename(E0 = e0, 
				 Emax = emax) %>%
	mutate('Residual SD' = D1$residSD) %>%
	pairs(.)

rm(list=objects())

###################################################
#### second continuous example
#### ml, 3 parm, meanlevs input
####
library(clinDR)
library(tidyverse)

set.seed(12357)

# Design aspects
nsim <- 48
doselev <- c(0,5,25,50,100)
n <- c(78,81,81,81,77)
Ndose <- length(doselev)

# mean values
meanlev <- c(2.465375,2.750120,3.581204,4.222137,4.929967)
resSD <- 7.967897

gen <- FixedMean(n, doselev, meanlev, resSD)



# Simulate outcomes
D1 <- emaxsim(nsim, gen, modType=3, nproc=8, binary = FALSE)

summary(D1, testalph=0.05)
plot(D1)

# Additional plots of simulation results
library(tidyverse,quietly = TRUE)

D1$est3 %>%
	as_tibble() %>%
	mutate(ED50 = exp(led50)) %>%
	select(-led50) %>%
	rename(E0 = e0, 
				 Emax = emax) %>%
	mutate('Residual SD' = D1$residSD) %>%
	pairs(.)

rm(list=objects())

##################################################
#### first binary test
#### bayes 4-parm, input parms
####
library(clinDR)
library(tidyverse)

set.seed(12357)

# Design aspects
nsim <- 56
doselev <- c(0,.1,.3,.6,1)
n <- c(100,100,100,100,100)
Ndose <- length(doselev)

# Parameters
led50 <- log(0.25)
lambda <- 1
e0 <- -0.5
target <- 1.12
targetDose <- 1
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
													 dTarget = 1,
													 p50 = 0.5,
													 parmDF = 5,
													 binary = TRUE)
# Stan and MCMC settings
rstan::rstan_options(auto_write = TRUE)
mcmc<-mcmc.control(chains=1,warmup=500,iter=5000,
									 propInit=0.15,adapt_delta = 0.95)

# Simulate outcomes
D1 <- emaxsimB(nsim, gen, prior,
							 modType=4,mcmc=mcmc,
							 check=FALSE, nproc=8, binary = TRUE)
summary(D1)
plot(D1)

# Additional plots of simulation results
library(tidyverse,quietly = TRUE)

D1$est %>%
	as_tibble() %>%
	mutate(ED50 = exp(led50)) %>%
	select(-led50) %>%
	mutate(e0 = plogis(e0),
				 emax = plogis(e0+emax)) %>%
	pairs(.)


rm(list=objects())


###
###  second binary example
###  input proportions, ml, 4-parms
###
library(clinDR)
library(tidyverse)

set.seed(12357)

# Design aspects
nsim <- 56
doselev <- c(0,5,25,50,100)
n <- c(100,100,100,100,100)
Ndose <- length(doselev)

# proportions
meanlev <- c(0.3775407,0.4750208,0.5655300,0.6196906,0.6502185)
gen <- FixedMean(n, doselev, meanlev, binary = TRUE)



# Simulate outcomes
D1 <- emaxsim(nsim, gen, modType=4, nproc=8, binary = TRUE)

summary(D1, testalph=0.05)
plot(D1)

# Additional plots of simulation results
library(tidyverse,quietly = TRUE)

D1$est4 %>%
	as_tibble() %>%
	mutate(ED50 = exp(led50)) %>%
	select(-led50) %>%
	pairs(.)
