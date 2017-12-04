\name{summary.emaxsimB}
\alias{summary.emaxsimB}

\title{ Summary of output of emaxsimB }
\description{
Detailed summary of repeated sampling properties of Bayesian Emax estimation
and comparison with simple pairwise comparisons.
}
\usage{
\method{summary}{emaxsimB}(object, testalpha = 0.05, 
	clev = c('0.95','0.9','0.8'),
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
repeated sampling properties of the estimation procedures.  It does not return any output.
}
\author{ Neal Thomas }

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
emax<-4.127726
sdy<-7.967897
pop<-c(log(ed50),emax,e0)    
meanlev<-emaxfun(doselev,pop)  

###FixedMean is specialized constructor function for emaxsim
gen<-FixedMean(n,doselev,meanlev,sdy)  

prior<-prior.control(epmu=0,epsd=30,emaxmu=0,emaxsd=30,p50=50,sigmalow=0.1,
		sigmaup=30,edDF=5)
mcmc<-mcmc.control(chains=1,warmup=500,iter=5000,seed=53453,propInit=0.15,adapt_delta = 0.95)

D1 <- emaxsimB(nsim,gen, prior, modType=3,mcmc=mcmc,check=FALSE)


summary(D1,testalph=0.05,clev='0.95',seSim=TRUE)
}
}
\keyword{nonlinear}
