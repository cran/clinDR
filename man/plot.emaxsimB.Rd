\name{plot.emaxsimB}
\alias{plot.emaxsimB}

\title{  Plot the output of emaxsimB}
\description{
  A Q-Q plot of the posterior mean of the mean dose response at a specified dose minus
the population value divided by the posterior SD of the mean difference.
}
\usage{
\method{plot}{emaxsimB}(x, id = x$idmax, difPlacebo= TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Output  of \code{\link{emaxsimB}} }
  \item{id}{ Index of the dose to be assessed (placebo index=1). }
  \item{difPlacebo}{ If true (default), the estimates and population
values are differences with placebo.  IF false, absolute dose response
values are used.}
  \item{\dots}{ Optional parameters passed to the plotting function }
}

\value{
ggplot object is returned
}

\author{Neal Thomas }

\seealso{  \code{\link{emaxsimB}}, \code{\link{print.emaxsimB}},
\code{\link{summary.emaxsimB}}  }

\examples{
\dontrun{
## emaxsimB changes the random number seeds
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

plot(D1,id=3)
}
}
\keyword{Bayes}
\keyword{Emax}
