\name{checkMonoEmax}
\alias{checkMonoEmax}

\title{Bayes posterior predictive test for Emax (monotone) model fit}
\description{Bayes posterior predictive test for an Emax (monotone) model  
fit comparing the best response from lower doses to
the response from the highest dose.}
  
\usage{
checkMonoEmax(y, 
      dose, 
      parm, 
      sigma2,
      nvec=rep(1,length(dose)),
      xbase=NULL,
      modelFun=emaxfun,
      trend='positive',
      binary= FALSE,logit=binary)
}

\arguments{
  \item{y}{Outcomes.  Continuous \code{y} can be individual data or group means.
  Binary \code{y} can be individual data, group proportions, or 0/1 data with
  correspondng counts, as is required by \code{fitEmaxB}.}
  \item{dose}{Doses corresponding to outcomes }
  \item{parm}{Matrix of simultated parameter values (each row is a
simulated parameter vector).  The \code{parm} values
must be constructed for use in the model function \code{modFun}.  The
default is a 4-parameter Emax model with parameters
(log(ED50),lambda,Emax,E0).  For a 3-parameter model, set lambda=1 for each
simulated parameter vector.}
  \item{sigma2}{Simulated draws from the residual variance (assumed
additive, homogeneous).  The length of \code{sigma2} must be the same
as the number of rows of \code{parm}.  \code{sigma2} is ignored when
\code{binary=TRUE}}
  \item{nvec}{The number of observations contributing to each \code{y}.  The default 
  is \code{1} for patient-level data.}
  \item{xbase}{Optional covariates matching \code{y}.  \code{nvec} must be 1 (patient-level)
  data.  The coeficients for \code{xbase} are the final columns of \code{parm}.}
  \item{modelFun}{The mean model function.  The first argument is a
scalar dose, and the second argument is a matrix of parameter values.
The rows of the matrix are random draws of parameter vectors for the
model.  The default function is the 4-parameter Emax function \code{emaxfun}.}
  \item{trend}{The default is 'positive', so high values for lower doses
  yield small Bayesian predictive probabilities.  Set \code{trend} to 'negative'
  for dose response curves with negative trends.}
  \item{binary}{If TRUE, the inverse logit 
  transform is applied to the
  (Emax) function output for comparison to dose group sample proportions, and the
  predictive data are sampled from a binomial distribution.}
  \item{logit}{\code{logit} is deprecated, use \code{binary}}
}
\details{
  A sample of parameters from the joint posterior distribution must be 
  supplied (typically produced by an MCMC program).  The Bayesian predictive p-value 
  is the posterior probability that a dose group
  sample mean in a new study with the same sample sizes would yield a 
  higher (or lower for negative trend) difference for one of the lower doses versus 
  the highest dose than was actually obtained from the real sample.  There must
  be at least two non-placebo dose groups (NA returned otherwise).  Placebo
  response is excluded from the comparisons.
  
  The function generates random numbers, so the random number generator/seed must
  be set before the function is called for exact reproducibility.
}
\value{
Returns a scalar Bayesian predictive p-value.
}

\author{ Neal Thomas }

\seealso{\code{\link{plot.plotB}}, \code{\link{plotD}},  \code{\link{plot.fitEmax}} }

\examples{
\dontrun{

data("metaData")
exdat<-metaData[metaData$taid==6 & metaData$poptype==1,]

prior<-emaxPrior.control(epmu=0,epsca=10,difTargetmu=0,difTargetsca=10,dTarget=80.0,
        p50=3.75,sigmalow=0.01,sigmaup=20)
mcmc<-mcmc.control(chains=3)

msSat<-sum((exdat$sampsize-1)*(exdat$sd)^2)/(sum(exdat$sampsize)-length(exdat$sampsize))
fitout<-fitEmaxB(exdat$rslt,exdat$dose,prior,modType=4,
				count=exdat$sampsize,msSat=msSat,mcmc=mcmc)
parms<-coef(fitout)[,1:4]  #use first intercept

checkMonoEmax(y=exdat$rslt, dose=exdat$dose, parm=parms, sigma2=(sigma(fitout))^2,
      nvec=exdat$sampsize, trend='negative')
      
}
\dontshow{

data("metaData")
exdat<-metaData[metaData$taid==6 & metaData$poptype==1,]


parms<-cbind(rnorm(5,1.85,1.9),rnorm(5,.5,.17),rnorm(5,-70.19,11),rnorm(5,2.3,2.7))
sig2<-rnorm(5,9.6,0.9)

checkMonoEmax(y=exdat$rslt, dose=exdat$dose, parm=parms, sigma2=sig2^2,
      nvec=exdat$sampsize, trend='negative')
}
}
\keyword{nonlinear}
