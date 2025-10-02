\name{print.emaxsimB}
\alias{print.emaxsimB}
\title{Print simulation output from emaxsimB}
\description{
  Prints key summary variables of Emax estimation peformance for each
simulation.  Can be used to identify simulated data sets yielding
unusual estimates.
}
\usage{
\method{print}{emaxsimB}(x, 
      nprint = min(nsim, 20), 
      id = x$idmax, 
      digits = 3, ...)
}

\arguments{
  \item{x}{Output of \code{\link{emaxsimB}} }
  \item{nprint}{ Number of simulations to print.  If a vector of
length 2, \code{nprint} is the range of simulations to print. }
  \item{id}{Output includes the stdBias for the dose with index \code{id}  vs
placebo }
  \item{digits}{ Number of decimal digits to print for Z and p-values }
  \item{\dots}{ Other print parameters (none currently implemented) }
}

\value{
Printed output returned as invisible matrix.
}
\author{ Neal Thomas }
\note{ 
The stdBias printed is the difference between the posterior mean of the 
dose response at the dose with index \code{id} and its population
value.  The difference is divided by the SE (posterior SD).
}
\seealso{ \code{\link{emaxsimB}}, \code{\link{summary.emaxsimB}},
\code{\link{plot.emaxsimB}} }
\examples{

\dontrun{
## emaxsimB changes the random number seed
nsim<-50
idmax<-5
doselev<-c(0,5,25,50,100)
n<-c(78,81,81,81,77)
Ndose<-length(doselev)

### population parameters for simulation
e0<-2.465375 
ed50<-67.481113 
dtarget<-100
diftarget<-2.464592
emax<-solveEmax(diftarget,dtarget,log(ed50),1,e0)

sdy<-7.967897
pop<-c(log(ed50),emax,e0)    
meanlev<-emaxfun(doselev,pop)  

###FixedMean is specialized constructor function for emaxsim
gen<-FixedMean(n,doselev,meanlev,sdy)  

prior<-emaxPrior.control(epmu=0,epsca=30,difTargetmu=0,
		difTargetsca=30,dTarget=100,p50=50,sigmalow=0.1,
		sigmaup=30,parmDF=5)
mcmc<-mcmc.control(chains=1,warmup=500,iter=5000,seed=53453,propInit=0.15,adapt_delta = 0.95)

D1 <- emaxsimB(nsim,gen, prior, modType=3,mcmc=mcmc,check=FALSE)

print(D1)
}
}
\keyword{nonlinear}
