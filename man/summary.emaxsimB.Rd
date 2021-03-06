\name{summary.emaxsimB}
\alias{summary.emaxsimB}

\title{ Summary of output of emaxsimB }
\description{
Detailed summary of repeated sampling properties of Bayesian Emax estimation
and comparison with simple pairwise comparisons.
}
\usage{
\method{summary}{emaxsimB}(object, testalpha = 0.05, 
	clev = c('0.9','0.95','0.8'),
	seSim = FALSE, ...)
}

\arguments{
  \item{object}{ Output  of \code{\link{emaxsimB}} }
  \item{testalpha}{ Alpha level for a one-sided MCP-MOD trend test.}
  \item{clev}{ Posterior proabilities for reported intervals}
  \item{seSim}{ If \code{TRUE}, then simulation standard errors are reported in
parentheses.  These should be distinguished from posterior SD
in the simulations.}
  \item{\dots}{Other unspecified parameters (none currently utilized) }
}

\value{
The function produces annotated output summarizing the
properties of the estimation procedures.  The summaries are also
returned as an invisible list for extracting results.
}
\author{ Neal Thomas }

\details{For pairwise comparisons, the 'most favorable pairwise comparison' means the dose with the best difference versus placebo is compared to the population
mean response for the selected dose, thus the target value for coverage,
bias, and RMSE changes depending on the selected dose.}

\seealso{  \code{\link{emaxsim}}, \code{\link{print.emaxsim}},
\code{\link{plot.emaxsim}}  }

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


summary(D1,testalph=0.05,clev='0.95')
}
}
\keyword{nonlinear}
