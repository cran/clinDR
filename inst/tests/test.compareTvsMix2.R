context('compare t and mixture pbo prior results')


library(RBesT)
library(MASS)
#library(ggplot2)


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

map_mcmc <- gMAP(cbind(y, y.se) ~ 1 | study,
                 weights = n, data = h.data,
                 family = gaussian,
                 beta.prior = cbind(0, crohn_sigma),
                 tau.dist = "HalfNormal", tau.prior = cbind(0, crohn_sigma / 2)
                 
)
#print(map_mcmc)

# Approximate MCMC samples with as a two component normal mixtures

map <- mixfit(map_mcmc, Nc=2)
# print(map)
# round(ess(map, method = 'elir'))

# Create a robustified MAP as a three component normal mixtures

map_robust <- robustify(map, weight = 0.5, mean = -50)
#print(map_robust)

#round(ess(map_robust, method = 'elir'))

# Generate a t-prior for map samples
# Extract the posterior MCMC sample as a matrix

map_samples <- as.matrix(map_mcmc)
prior.sample <- map_samples[,"theta_pred"]

# Fit t-distribution

df.t <- 5
fit.prior.t <- fitdistr(prior.sample, "t", df=df.t)

#t-prior for placebo

#fit.prior.t$estimate

#Parameters of t-distribution prior

m.t <- fit.prior.t$estimate["m"]
s.t <- fit.prior.t$estimate["s"]

# Add robustified mixture normal distribution from map_robust

w.map <- map[1,]
mu.map <- map[2,]
sd.map <- map[3,]

w.map_robust <- map_robust[1,]
mu.map_robust <- map_robust[2,]
sd.map_robust <- map_robust[3,]

##########################################################
# Compare densities graphically
#
# hist(prior.sample , density = 20, breaks=1000, probability = TRUE)
# 
# # Add t distribution
# 
# curve(dt((x - m.t) / s.t, df = df.t) / s.t, add = TRUE, col = "red", lwd = 2)
# 
# # Add mixture normal distribution from map
# 
# curve(mixnormplot(x, w.map, mu.map, sd.map), add = TRUE, col = "blue", lwd = 2)
# 
# curve(mixnormplot(x, w.map_robust, mu.map_robust,sd.map_robust), add = TRUE, col = "green", lwd = 2)
# 
# legend("topright", legend = c("t-prior", "MAP", "Robust"), col = c("red", "blue", "green"), lwd = 2)

# Generate simulated test data

e0 <- summary(map_mcmc)$theta.pred[1]
sdy <- summary(map_mcmc)$theta.pred[2]

# No conflict

no.conflict <- simdata(e0, sdy)
Dose<-no.conflict$Dose
y<-no.conflict$Response
nc.summ<-data.frame(Dose=sort(unique(Dose)),n=as.vector(table(Dose)),
                    mn=tapply(y,Dose,mean),
                    sd=tapply(y,Dose,sd),row.names=NULL)

# Moderate conflict
mod.conflict <- simdata(e0 - 2*sdy, sdy)
Dose<-mod.conflict$Dose
y<-mod.conflict$Response
mod.summ<-data.frame(Dose=sort(unique(Dose)),n=as.vector(table(Dose)),
                    mn=tapply(y,Dose,mean),
                    sd=tapply(y,Dose,sd),row.names=NULL)

# Severe conflict
sev.conflict <- simdata(e0 - 6* sdy, sdy)

Dose<-sev.conflict$Dose
y<-sev.conflict$Response
sev.summ<-data.frame(Dose=sort(unique(Dose)),n=as.vector(table(Dose)),
                    mn=tapply(y,Dose,mean),
                    sd=tapply(y,Dose,sd),row.names=NULL)

# Prior Emax Control for t, map and robust
prior.emax.t <- emaxPrior.control(epmu=m.t,epsca=s.t, 
                                  difTargetmu = -40, difTargetsca =20,
                                  dTarget =100, p50=15, 
                                  sigmalow =10, sigmaup = 100)

prior.emax.map <- emaxPrior.control(mixP=2, w_ep=w.map,mu_ep=mu.map, sd_ep=sd.map,
                                    difTargetmu = -40, difTargetsca =20,
                                    dTarget =100, p50=15, 
                                    sigmalow =10, sigmaup = 100)

prior.emax.rob <- emaxPrior.control(mixP=3, w_ep=w.map_robust,mu_ep=mu.map_robust, sd_ep=sd.map_robust,
                                    difTargetmu = -40, difTargetsca =20,
                                    dTarget =100, p50=15, 
                                    sigmalow =10, sigmaup = 100)


# MCMC set up

mcmc <- mcmc.control(chains=3,warmup =500,iter=1500)

# Pool variances
ms.nc <- poolvar(nc.summ)
ms.mod <- poolvar(mod.summ)
ms.sev <- poolvar(sev.summ)

fit.nc.t <- fitEmaxB(y=nc.summ$mn, dose = nc.summ$Dose, prior=prior.emax.t,
                     modType =4, count= nc.summ$n, msSat=ms.nc, mcmc=mcmc )
e0.nc.t <- coef(fit.nc.t)[,'e0[1]']

fit.nc.map <- fitEmaxB(y=nc.summ$mn, dose = nc.summ$Dose, prior=prior.emax.map,
                       modType =4, count= nc.summ$n, msSat=ms.nc, mcmc=mcmc )
e0.nc.map <-coef(fit.nc.map)[,'e0[1]'] 

fit.nc.rob <- fitEmaxB(y=nc.summ$mn, dose = nc.summ$Dose, prior=prior.emax.rob,
                       modType =4, count= nc.summ$n, msSat=ms.nc, mcmc=mcmc )
e0.nc.rob <- coef(fit.nc.rob)[,'e0[1]']

### test t vs mixture
doselev<-nc.summ$Dose
tres<-predict(fit.nc.t,doselev)
mixres <-predict(fit.nc.map,doselev)


compEst<- (tres$pred-mixres$pred)/tres$se
test_that("t and mixture prior DR estimates are similar: no conflict", {
  expect_true(all(abs(compEst)<0.1))
})
compSe<-mixres$se/tres$se
test_that("t and mixture prior DR SE are similar: no conflict", {
  expect_true(all(compSe>0.9 & compSe<1.1))
})

### compare pbo posteriors graphically
# e0.mn <- nc.summ$mn[1]
# e0.se <- nc.summ$sd[1]/sqrt(nc.summ$n[1])
# x.plc <- seq(e0.mn - 2*e0.se , e0.mn +2*e0.se,length.out=1000)
# y.plc <- dnorm(x.plc, mean=e0.mn, sd=e0.se)
# plot(x.plc, y.plc, type = "l", lwd = 2, col = "grey", ylim=c(0, 0.2))
# 
# # Fill the area under the curve
# polygon(c(x.plc, rev(x.plc)), c(y.plc, rep(0, length(y.plc))),
#         col = rgb(0,0,0,0.2), border = NA)
# 
# lines(density(e0.nc.t), col = "red", lwd = 2)
# lines(density(e0.nc.map), col = "blue", lwd = 2)
# lines(density(e0.nc.rob), col = "green", lwd = 2)
# legend("topright",
#        legend = c("True","t", "map", "rob"),
#        col = c("grey","red","blue", "green"),
#        lwd = 2)

# Moderate Conflict

fit.mod.t <- fitEmaxB(y=nc.summ$mn, dose = nc.summ$Dose, prior=prior.emax.t,
                     modType =4, count= nc.summ$n, msSat=ms.nc, mcmc=mcmc )
e0.mod.t <- coef(fit.mod.t)[,'e0[1]']

fit.mod.map <- fitEmaxB(y=nc.summ$mn, dose = nc.summ$Dose, prior=prior.emax.map,
                       modType =4, count= nc.summ$n, msSat=ms.nc, mcmc=mcmc )
e0.mod.map <-coef(fit.mod.map)[,'e0[1]'] 

fit.mod.rob <- fitEmaxB(y=nc.summ$mn, dose = nc.summ$Dose, prior=prior.emax.rob,
                       modType =4, count= nc.summ$n, msSat=ms.nc, mcmc=mcmc )
e0.mod.rob <- coef(fit.mod.rob)[,'e0[1]']

### test t vs mixture
tres<-predict(fit.mod.t,doselev)
mixres <-predict(fit.mod.map,doselev)


compEst<- (tres$pred-mixres$pred)/tres$se
test_that("t and mixture prior DR estimates are similar: moderate conflict", {
  expect_true(all(abs(compEst)<0.1))
})
compSe<-mixres$se/tres$se
test_that("t and mixture prior DR SE are similar:  moderate conflict", {
  expect_true(all(compSe>0.9 & compSe<1.1))
})


### compare pbo posteriors graphically
# e0.mn <- mod.summ$mn[1]
# e0.se <- mod.summ$sd[1]/sqrt(mod.summ$n[1])
# x.plc <- seq(e0.mn - 2*e0.se , e0.mn +2*e0.se,length.out=1000)
# y.plc <- dnorm(x.plc, mean=e0.mn, sd=e0.se)
# plot(x.plc, y.plc, type = "l", lwd = 2, col = "grey", ylim=c(0, 0.2))
# 
# # Fill the area under the curve
# polygon(c(x.plc, rev(x.plc)), c(y.plc, rep(0, length(y.plc))),
#         col = rgb(0,0,0,0.2), border = NA)
# 
# lines(density(e0.mod.t), col = "red", lwd = 2)
# lines(density(e0.mod.map), col = "blue", lwd = 2)
# lines(density(e0.mod.rob), col = "green", lwd = 2)
# legend("topright",
#        legend = c("True","t", "map", "rob"),
#        col = c("grey","red","blue", "green"),
#        lwd = 2)


###!! Satrajit test with severe conflict (rob prior checks)
