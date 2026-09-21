context('compare t and mixture pbo prior results')

# Functions for comparison test

mixnormplot <- function(x, probs, mu, sd) {
  
  # weighted sum of normal densities
  
  rowSums(sapply(seq_along(probs), function(k)
    
    probs[k] * dnorm(x, mean = mu[k], sd = sd[k])
    
  ))
  
}

##########
#### continuous data example
####
# Generate data

simdata <- function(e0,sdy){
  # Dose
  
  doselev <- c(0, 5, 25, 50, 100)
  
  # Sample Size
  
  ss <- c(35, 70, 70, 70, 70)
  
  # SD
  sdy <- sdy
  
  # Other Emax parameters , lambda=1
  
  ed50 <- 15                
  dtarget <- max(doselev)   
  diftarget <- -40
  e0 <- e0
  emax <- solveEmax(target =diftarget, dose=dtarget, led50 = log(ed50), lambda=1, e0)
  pop <- c(log(ed50), emax, e0)
  
  # Generate Data
  
  meanlev <- emaxfun(dose=doselev, parm=pop)
  meanvec <- rep(meanlev, ss)
  y <- rnorm(sum(ss), meanvec, sdy)
  dose <- rep(doselev, times = ss)
  
  simdata <- data.frame(Dose = dose, Response = y)
  
  return(simdata)
  
}

# Pooled variance from summary data

poolvar<- function(summ.data){
  msSat <- sum((summ.data$n-1)*(summ.data$sd)^2)/(sum(summ.data$n)-length(summ.data$n))
  return(msSat)
}

##### example from RBesT
# Use Crohn's disease example to build an informative MAP prior

h.data <- crohn
crohn_sigma <- 88
h.data$y.se <- crohn_sigma / sqrt(h.data$n)

# Run gMAP to generate MCMC samples for MAP prior for new study

map_mcmc <- suppressWarnings(gMAP(cbind(y, y.se) ~ 1 | study,
                 weights = n, data = h.data,
                 family = gaussian,
                 beta.prior = cbind(0, crohn_sigma),
                 tau.dist = "HalfNormal", tau.prior = cbind(0, crohn_sigma / 2)
                 ))


# Approximate MCMC samples with as a two component normal mixtures

map <- mixfit(map_mcmc, Nc=2)



# Generate a t-prior for map samples
# Extract the posterior MCMC sample as a matrix

map_samples <- posterior::as_draws_matrix(map_mcmc)
prior.sample <- map_samples[,"theta_pred"]

# Fit t-distribution

df.t <- 5
suppressWarnings(fit.prior.t <- fitdistr(prior.sample, "t", df=df.t))

#t-prior for placebo

#fit.prior.t$estimate

#Parameters of t-distribution prior

m.t <- fit.prior.t$estimate["m"]
s.t <- fit.prior.t$estimate["s"]

# Add robustified mixture normal distribution from map_robust

w.map <- map[1,]
mu.map <- map[2,]
sd.map <- map[3,]

##########################################################


# Generate simulated test data

e0 <- summary(map_mcmc)$theta.pred[[1]]
sdy <- summary(map_mcmc)$theta.pred[[3]]

# No conflict

no.conflict <- simdata(e0, sdy)
Dose<-no.conflict$Dose
y<-no.conflict$Response
nc.summ<-data.frame(Dose=sort(unique(Dose)),n=as.vector(table(Dose)),
                    mn=tapply(y,Dose,mean),
                    sd=tapply(y,Dose,sd),row.names=NULL)



# Prior Emax Control for t, map and robust
prior.emax.t <- clinDR::emaxPrior.control(epmu=m.t,epsca=s.t, 
                                  difTargetmu = -40, difTargetsca =20,
                                  dTarget =100, p50=15, 
                                  sigmalow =10, sigmaup = 100)
prior.emax.map <- clinDR::emaxPrior.control(mixP=3, w_ep=w.map,mu_ep=mu.map, sd_ep=sd.map,
                                            difTargetmu = -40, difTargetsca =20,
                                            dTarget =100, p50=15, 
                                            sigmalow =10, sigmaup = 100)


# MCMC set up

mcmc <- mcmc.control(chains=3,warmup =500,iter=1500)

# Pool variances
ms.nc <- poolvar(nc.summ)

# No conflict

fit.nc.t <- fitEmaxB(y=nc.summ$mn, dose = nc.summ$Dose, prior=prior.emax.t,
                     modType =4, count= nc.summ$n, msSat=ms.nc, mcmc=mcmc )

fit.nc.map <- fitEmaxB(y=nc.summ$mn, dose = nc.summ$Dose, prior=prior.emax.map,
                       modType =4, count= nc.summ$n, msSat=ms.nc, mcmc=mcmc )


doselev<-nc.summ$Dose
tres<-predict(fit.nc.t,doselev)
mapres<- predict(fit.nc.map, doselev)

### test t vs map
compEst<- (tres$pred-mapres$pred)/tres$se
test_that("t and mixture prior DR estimates are similar: no conflict", {
  expect_true(all(abs(compEst)<0.1))
})
compSe<-mapres$se/tres$se
test_that("t and mixture prior DR SE are similar:  no conflict", {
  expect_true(all(compSe>0.9 & compSe<1.1))
})
